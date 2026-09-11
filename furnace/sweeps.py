"""Parametric sweep driver: mutates param copies, never edits params.yaml. Step 6."""

from __future__ import annotations

import csv
import json
import math
import sys
import time
import types
from pathlib import Path
from typing import Any

import openmc

from furnace.triso import particle_atom_counts, _stage_outer_radius, _shell_vol
from furnace.geometry import frustum_volume, _particle_effective_density

_RESULTS_DIR = Path(__file__).parent.parent / 'results'
_CONVERGENCE_CSV = _RESULTS_DIR / 'convergence_n_slabs.csv'

# Fixed discretisation bias from n_slabs=8 convergence study (|k(8)−k(16)|).
# Direction: staircase overestimates k (particles overflow from cone to cylinder);
# bias is conservative — adding it to k+2σ makes the upper bound more conservative.
# Loaded once at import from convergence_n_slabs.csv; 0.0 if the file is absent.
def _load_dk_disc(n_slabs: int) -> float:
    if not _CONVERGENCE_CSV.exists():
        return 0.0
    with _CONVERGENCE_CSV.open() as f:
        for row in csv.DictReader(f):
            try:
                if int(row['n_slabs']) == n_slabs and row.get('dk_disc', '').strip():
                    return float(row['dk_disc'])
            except (ValueError, KeyError):
                pass
    return 0.0


_DK_DISC_N8 = _load_dk_disc(8)

_CSV_FIELDNAMES = [
    'tag', 'overrides_json',
    'k_eff', 'sigma', 'k_plus_2sigma', 'k_plus_3sigma',
    'dk_disc', 'k_plus_2sigma_plus_dk_disc',
    'u235_mass_g', 'state', 'stage',
    'z_flood_cm', 'water_density_gcc',
    'n_slabs', 'bed_volume_cm3', 'pf_achieved', 'bed_height_cm',
    'h_per_u235', 'c_per_u235', 'thermal_flux_fraction',
    'wall_time_s', 'seed', 'openmc_version',
]

# Quick-mode overrides: small charge mass (5 g ≈ 12 000 particles, matches smoke test),
# 1 000 transport particles, 2 inactive + 10 active = 12 total batches.
# Geometry + transport finishes in ~30 s per case. Results have no NCS value —
# quick mode is only for confirming the pipeline runs end-to-end.
_QUICK = {
    'model': {'particles': 1000, 'batches': 12, 'inactive': 2},
    'dimensions': {'bed': {'charge_mass_g': 5.0}},
}


# ---------------------------------------------------------------------------
# Params helpers
# ---------------------------------------------------------------------------

def _unfreeze(obj: Any) -> Any:
    """Recursively convert MappingProxyType / tuples to mutable dicts / lists."""
    if isinstance(obj, types.MappingProxyType):
        return {k: _unfreeze(v) for k, v in obj.items()}
    if isinstance(obj, dict):
        return {k: _unfreeze(v) for k, v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [_unfreeze(v) for v in obj]
    return obj


def _freeze(obj: Any) -> Any:
    """Recursively wrap dicts in MappingProxyType and lists in tuples."""
    if isinstance(obj, dict):
        return types.MappingProxyType({k: _freeze(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return tuple(_freeze(item) for item in obj)
    return obj


def _merge(base: dict, patch: dict) -> None:
    """Deep-merge patch into base in-place; underscore-prefixed keys are skipped."""
    for k, v in patch.items():
        if k.startswith('_'):
            continue
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _merge(base[k], v)
        else:
            base[k] = v


# ---------------------------------------------------------------------------
# Vessel capacity
# ---------------------------------------------------------------------------

def vessel_capacity_g(params, state: str, stage: str) -> float:
    """Maximum charge mass (g) that fits inside the retort at the given packing.

    Available volume = cone frustum (cone base to retort ID) + retort cylinder
    (from cone top to vacuum boundary at retort height).  Everything above that
    is outside the vacuum boundary, so any charge that would put the bed above
    z_rt is not tested.
    """
    dim = params['dimensions']
    mdl = params['model']

    r_thr   = dim['nozzle']['throat_diameter_cm'] / 2.0
    r_ret   = dim['retort']['id_cm'] / 2.0
    z_ct    = dim['cone']['vertical_drop_cm']
    z_rt    = z_ct + dim['retort']['height_cm']
    half_ang = dim['cone']['included_angle_deg'] / 2.0

    v_frustum = frustum_volume(r_ret, r_thr, half_angle_deg=half_ang)
    v_cyl     = math.pi * r_ret**2 * (z_rt - z_ct)
    v_vessel  = v_frustum + v_cyl

    pf = float(mdl['packing_fraction_static'] if state == 'collapsed'
               else mdl['packing_fraction_static'] / mdl['bed_expansion_ratio'])

    rho_eff = _particle_effective_density(stage, params)
    return v_vessel * pf * rho_eff


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------

def _compute_hc_ratios(
    stats,
    params,
    stage: str,
    fill_mat: openmc.Material,
) -> tuple[float, float]:
    """H/²³⁵U and C/²³⁵U atom ratios from actual model materials and volumes.

    Solid-phase C and U-235 come from particle_atom_counts() applied to the
    modified params; gas-phase H and C come from the OpenMC material atom
    densities so that any density or composition override is automatically
    reflected.
    """
    n_p = stats.n_particles
    N_C_solid, _, N_U235 = particle_atom_counts(stage, params)
    N_C_solid *= n_p
    N_U235 *= n_p

    # Void volume from actual achieved packing fraction
    V_gas = stats.V_bulk_cm3 * (1.0 - stats.pf_achieved)

    N_H = 0.0
    N_C_gas = 0.0
    for nuc, nd in fill_mat.get_nuclide_atom_densities().items():
        n_cm3 = nd * 1e24   # atoms/b-cm → atoms/cm³
        if nuc[0] == 'H':
            N_H += n_cm3 * V_gas
        elif nuc[0] == 'C':
            N_C_gas += n_cm3 * V_gas

    if N_U235 <= 0.0:
        return 0.0, 0.0
    return N_H / N_U235, (N_C_solid + N_C_gas) / N_U235


def _thermal_flux_fraction(sp: openmc.StatePoint) -> float:
    """Thermal / total flux fraction from the flux_spectrum_bed tally (3 groups)."""
    try:
        tally = sp.get_tally(name='flux_spectrum_bed')
        flux = tally.get_values(scores=['flux']).flatten()
        total = float(flux.sum())
        return float(flux[0]) / total if total > 0.0 else float('nan')
    except Exception:
        return float('nan')


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def _tag_in_csv(csv_path: Path, tag: str) -> bool:
    if not csv_path.exists():
        return False
    with csv_path.open() as f:
        return any(row.get('tag') == tag for row in csv.DictReader(f))


def _read_row(csv_path: Path, tag: str) -> dict:
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            if row.get('tag') == tag:
                return dict(row)
    return {}


def _append_to_csv(csv_path: Path, row: dict) -> None:
    """Append one row; write header if the file is new. Uses fcntl for crash safety."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not csv_path.exists()
    with csv_path.open('a', newline='') as f:
        try:
            import fcntl
            fcntl.flock(f, fcntl.LOCK_EX)
        except ImportError:
            pass  # Windows: no fcntl; single-process use is safe without it
        w = csv.DictWriter(f, fieldnames=_CSV_FIELDNAMES, extrasaction='ignore')
        if write_header:
            w.writeheader()
        w.writerow(row)
        f.flush()
        try:
            import fcntl
            fcntl.flock(f, fcntl.LOCK_UN)
        except ImportError:
            pass


# ---------------------------------------------------------------------------
# OpenMC run helper
# ---------------------------------------------------------------------------

def _export_and_run(
    model: openmc.Model,
    out_dir: Path,
    threads: int | None,
    mpi_args: list[str] | None,
) -> Path:
    """Export XML and run OpenMC with optional thread/MPI parallelism."""
    out_dir.mkdir(parents=True, exist_ok=True)
    model.geometry.export_to_xml( path=str(out_dir / 'geometry.xml'))
    model.materials.export_to_xml(path=str(out_dir / 'materials.xml'))
    model.settings.export_to_xml( path=str(out_dir / 'settings.xml'))
    model.tallies.export_to_xml(  path=str(out_dir / 'tallies.xml'))
    openmc_exec = str(Path(sys.executable).parent / 'openmc')
    openmc.run(
        cwd=str(out_dir),
        openmc_exec=openmc_exec,
        threads=threads,
        mpi_args=mpi_args,
    )
    return out_dir / f'statepoint.{model.settings.batches}.h5'


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_case(
    base_params,
    overrides: dict,
    tag: str,
    run_dir: Path,
    *,
    csv_path: Path | None = None,
    threads: int | None = None,
    mpi_args: list[str] | None = None,
    quick: bool = False,
    force: bool = False,
) -> dict:
    """Build and run one eigenvalue case; append result to CSV.

    Parameters
    ----------
    base_params : frozen MappingProxyType from load_params()
    overrides   : nested dict mirroring params structure; keys starting with '_'
                  map to build_model kwargs: '_state', '_stage', '_z_flood', '_water_density'
    tag         : unique identifier used as CSV key and directory label
    run_dir     : directory for OpenMC XML files and statepoint output
    csv_path    : results CSV (default: results/sweep.csv)
    threads     : OpenMC shared-memory threads per case (None → OMP_NUM_THREADS)
    mpi_args    : MPI launch command, e.g. ['mpiexec', '-n', '4']
    quick       : use reduced particles/batches to smoke-test sweep definitions
    force       : re-run even if tag already exists in CSV

    Returns
    -------
    dict with all result columns including 'overrides_json'
    """
    from furnace.model import build_model

    if csv_path is None:
        csv_path = _RESULTS_DIR / 'sweep.csv'

    if not force and _tag_in_csv(csv_path, tag):
        return _read_row(csv_path, tag)

    # Build modified params: quick defaults first so per-case overrides win.
    # (Otherwise a mass sweep's charge_mass_g would be wiped out by _QUICK.)
    params_dict = _unfreeze(base_params)
    if quick:
        _merge(params_dict, _QUICK)
    _merge(params_dict, overrides)
    params = _freeze(params_dict)

    # Build-model kwargs encoded as underscore-prefixed override keys
    state         = overrides.get('_state',         'fluidized')
    stage         = overrides.get('_stage',         'bare_kernel')
    z_flood       = overrides.get('_z_flood',       None)
    water_density = overrides.get('_water_density', 1.0)

    t0 = time.perf_counter()

    model, stats = build_model(
        params, state=state, stage=stage,
        z_flood=z_flood, water_density_gcc=water_density,
    )
    sp_path = _export_and_run(model, Path(run_dir), threads=threads, mpi_args=mpi_args)

    wall_time = time.perf_counter() - t0

    with openmc.StatePoint(str(sp_path)) as sp:
        keff = sp.keff
        k = float(keff.nominal_value)
        s = float(keff.std_dev)
        thermal_ff = _thermal_flux_fraction(sp)

    all_mats = {m.name: m for m in model.geometry.get_all_materials().values()}
    if z_flood is not None:
        fill_mat = next(m for name, m in all_mats.items() if name.startswith('water_'))
    else:
        fill_mat = all_mats['process_gas']

    h_per_u235, c_per_u235 = _compute_hc_ratios(stats, params, stage, fill_mat)

    n_slabs = int(params['model']['n_slabs'])
    dk_disc = _load_dk_disc(n_slabs)

    overrides_public = {kk: vv for kk, vv in overrides.items() if not kk.startswith('_')}

    row = {
        'tag':                        tag,
        'overrides_json':             json.dumps(overrides_public),
        'k_eff':                      k,
        'sigma':                      s,
        'k_plus_2sigma':              k + 2.0 * s,
        'k_plus_3sigma':              k + 3.0 * s,
        'dk_disc':                    dk_disc,
        'k_plus_2sigma_plus_dk_disc': k + 2.0 * s + dk_disc,
        'u235_mass_g':                stats.u235_mass_g,
        'state':                      state,
        'stage':                      stage,
        'z_flood_cm':                 z_flood if z_flood is not None else '',
        'water_density_gcc':          water_density if z_flood is not None else '',
        'n_slabs':                    n_slabs,
        'bed_volume_cm3':             stats.V_bulk_cm3,
        'pf_achieved':                stats.pf_achieved,
        'bed_height_cm':              stats.bed_height_cm,
        'h_per_u235':                 h_per_u235,
        'c_per_u235':                 c_per_u235,
        'thermal_flux_fraction':      thermal_ff,
        'wall_time_s':                wall_time,
        'seed':                       int(params['model']['seed']),
        'openmc_version':             openmc.__version__,
    }

    _append_to_csv(csv_path, row)
    return row
