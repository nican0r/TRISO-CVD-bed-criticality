"""OpenMC Model assembly: geometry + materials + settings + tallies. Step 5."""

from __future__ import annotations

import math
import shutil
import sys
from pathlib import Path
from typing import NamedTuple

import openmc
import openmc.model

from furnace.materials import (
    uco_kernel as _uco_kernel,
    buffer_pyc as _buffer_pyc,
    ipyc as _ipyc,
    sic as _sic,
    opyc as _opyc,
    graphite_structural as _graphite_structural,
    process_gas as _process_gas,
    water as _water,
)
from furnace.geometry import bed_region, furnace_shell_cells

_OPENMC_EXEC = shutil.which('openmc') or str(Path(sys.executable).parent / 'openmc')

# U-235 thermal Watt fission spectrum: a in eV, b in 1/eV (ENDF/B-VIII.0 defaults)
_WATT_A = 0.988e6    # eV
_WATT_B = 2.249e-6   # 1/eV

# Broad energy group boundaries [eV]: thermal / epithermal / fast
# Thermal:    0 – 0.625 eV  (conventional 2200 m/s boundary)
# Epithermal: 0.625 eV – 1 MeV
# Fast:       1 MeV – 20 MeV  (upper limit of ENDF/B library)
_ENERGY_BINS_EV = [0.0, 0.625, 1.0e6, 20.0e6]


class BedStats(NamedTuple):
    n_particles: int
    pf_achieved: float
    V_bulk_cm3: float
    bed_height_cm: float
    u235_mass_g: float


def build_model(
    params,
    *,
    state: str = 'fluidized',
    stage: str = 'bare_kernel',
    z_flood: float | None = None,
    water_density_gcc: float = 1.0,
    n_inactive: int | None = None,
    n_active: int | None = None,
    n_particles: int | None = None,
    seed: int | None = None,
    charge_mass_g: float | None = None,
) -> tuple[openmc.Model, BedStats]:
    """Assemble and return an openmc.Model for the CVD furnace NCS eigenvalue case.

    Parameters
    ----------
    params        : frozen params dict from load_params()
    state         : 'fluidized' (nominal) or 'collapsed'
    stage         : TRISO deposition stage; 'bare_kernel' for nominal
    z_flood          : bottom-up flood level (cm, absolute z). When set, all bed/cone cells
                       whose z_bot < z_flood receive water; cells above receive gas.
                       The gas-above-bed cell is also split at z_flood if needed.
                       When None, process gas fills the entire bed and retort interior.
    water_density_gcc: density of flood water (g/cm³); default 1.0 (liquid)
    n_inactive    : inactive batches; defaults to params['model']['inactive']
    n_active      : active batches; defaults to params['model']['batches'] - inactive
    n_particles   : particles per batch; defaults to params['model']['particles']
    seed          : OpenMC seed AND packing seed (both change together for seed study)
    charge_mass_g : overrides params bed charge mass for mass sweeps

    Returns
    -------
    (model, stats) — openmc.Model and BedStats namedtuple
    """
    mdl = params['model']
    dim = params['dimensions']

    _seed        = int(seed        if seed        is not None else mdl['seed'])
    _n_inactive  = int(n_inactive  if n_inactive  is not None else mdl['inactive'])
    _n_active    = int(n_active    if n_active    is not None else (mdl['batches'] - mdl['inactive']))
    _n_particles = int(n_particles if n_particles is not None else mdl['particles'])

    # ── Materials ────────────────────────────────────────────────────────────
    kernel_mat = _uco_kernel(params)
    buf_mat    = _buffer_pyc(params)
    ipyc_mat   = _ipyc(params)
    sic_mat    = _sic(params)
    opyc_mat   = _opyc(params)
    graph_mat  = _graphite_structural(params)
    water_mat  = _water(water_density_gcc)
    gas_mat    = _process_gas(params)

    bed_fill = water_mat if z_flood is not None else gas_mat

    # ── Bed ──────────────────────────────────────────────────────────────────
    bed = bed_region(
        params,
        state=state,
        stage=stage,
        n_slabs=int(mdl['n_slabs']),
        background_material=bed_fill,
        charge_mass_g=charge_mass_g,
        seed=_seed,
        z_flood=z_flood,
        dry_material=gas_mat if z_flood is not None else None,
    )

    # ── Furnace shell ─────────────────────────────────────────────────────────
    shell = furnace_shell_cells(params, graphite=graph_mat, water=water_mat)

    # ── Gas above the bed (retort cylinder interior above bed top) ───────────
    r_ret_in      = dim['retort']['id_cm'] / 2.0
    z_ct          = dim['cone']['vertical_drop_cm']
    z_rt          = z_ct + dim['retort']['height_cm']
    z_bed_top_val = bed['bed_height_cm']    # absolute z of bed top
    z_bed_bot     = bed['z_bed_bot_cm']     # 0.0: bed starts at cone base

    # Gas (or water) fills the retort cylinder above the bed.  The retort cylinder
    # only exists above z_ct, so start the cell there even if bed_top < z_ct.
    z_gas_above_bot = max(z_bed_top_val, z_ct)

    cyl_ret_in       = openmc.ZCylinder(r=r_ret_in)
    zp_gas_above_bot = openmc.ZPlane(z0=z_gas_above_bot)
    zp_rt_top        = openmc.ZPlane(z0=z_rt)   # interior top; vacuum is the graphite cap above

    # For bottom-up flood: split the above-bed cell at z_flood if the flood
    # level extends into the gas region.
    if z_flood is not None and z_flood > z_gas_above_bot:
        z_flood_capped = min(z_flood, z_rt)
        if z_flood_capped < z_rt:
            zp_flood_lvl = openmc.ZPlane(z0=z_flood_capped)
            gas_cells = [
                openmc.Cell(name='flood_above_bed', fill=water_mat,
                            region=-cyl_ret_in & +zp_gas_above_bot & -zp_flood_lvl),
                openmc.Cell(name='gas_above_bed', fill=gas_mat,
                            region=-cyl_ret_in & +zp_flood_lvl & -zp_rt_top),
            ]
        else:
            gas_cells = [
                openmc.Cell(name='gas_above_bed', fill=water_mat,
                            region=-cyl_ret_in & +zp_gas_above_bot & -zp_rt_top),
            ]
    else:
        gas_cells = [
            openmc.Cell(name='gas_above_bed', fill=gas_mat,
                        region=-cyl_ret_in & +zp_gas_above_bot & -zp_rt_top),
        ]

    # ── Geometry ──────────────────────────────────────────────────────────────
    all_cells = bed['cells'] + bed['cone_void_cells'] + gas_cells + shell['all']
    root_univ = openmc.Universe(cells=all_cells)
    geometry  = openmc.Geometry(root_univ)

    # Collect every material actually referenced by the geometry — this includes
    # the TRISO layer materials created inside particle_at_stage() (called from
    # bed_region()), which are different Python objects from the ones above and
    # would be missing from a manually assembled list.
    all_mats = list(geometry.get_all_materials().values())

    materials = openmc.Materials(all_mats)

    # ── Settings ──────────────────────────────────────────────────────────────
    # Shannon entropy mesh: 10×10×20 spanning cone base (z=0) to retort top.
    # Lateral cells 0.5×0.5 cm.  Cells inside the cone but outside the particle
    # region will be empty and contribute zero entropy — harmless.
    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left  = (-r_ret_in, -r_ret_in, 0.0)
    entropy_mesh.upper_right = ( r_ret_in,  r_ret_in, z_rt)
    entropy_mesh.dimension   = [10, 10, 20]

    # Watt source box spans the full bed (cone base z=0 to bed top).
    # constraints={'fissionable': True} rejects non-fissile sample positions;
    # only affects the first generation — active batches use the fission bank.
    src = openmc.IndependentSource(
        space=openmc.stats.Box(
            lower_left=(-r_ret_in, -r_ret_in, z_bed_bot),   # z_bed_bot = 0.0
            upper_right=( r_ret_in,  r_ret_in, z_bed_top_val),
        ),
        energy=openmc.stats.Watt(a=_WATT_A, b=_WATT_B),
        angle=openmc.stats.Isotropic(),
        constraints={'fissionable': True},
    )

    settings = openmc.Settings()
    settings.run_mode    = 'eigenvalue'
    settings.batches     = _n_inactive + _n_active
    settings.inactive    = _n_inactive
    settings.particles   = _n_particles
    settings.seed        = _seed
    settings.entropy_mesh = entropy_mesh
    settings.source      = [src]
    # Materials are set to 293.6 K (room-temperature bounding case; see materials.py).
    # 'nearest' snaps to the closest kT in the library for each nuclide:
    #   294 K  — U, H, O, N, Si, Cl, Ar isotopes (added via scripts/add_room_temp.py)
    #   900 K  — C12, C13 (nndc_hdf5 has only C0, no isotopic match; conservative)
    #   1200 K — c_Graphite S(α,β) (only temperature in endfb80_hdf5/thermal/)
    # 'default' = 900 K catches any nuclide below the library minimum.
    # 'range' must include 293.6 K so OpenMC accepts the material temperature before snapping.
    settings.temperature = {
        'method':  'nearest',
        'default': 900.0,
        'range':   [250.0, 3000.0],
    }
    # UCO kernels occupy ~3.2% of the source Box volume (pf=0.333 × kernel fraction
    # 0.123 × π/4 box-to-cylinder ratio).  The default source_rejection_fraction=0.05
    # (5%) would reject our valid geometry.  Setting to 0.005 allows acceptance rates
    # down to 0.5% — comfortably below the ~3.2% actual rate.
    settings.source_rejection_fraction = 0.005

    # ── Tallies ───────────────────────────────────────────────────────────────
    # MaterialFilter must reference the same objects that are in the geometry.
    # particle_at_stage() creates its own material instances (different IDs from
    # kernel_mat/buf_mat/… above), so look them up from the geometry by name.
    _mat_by_name = {m.name: m for m in all_mats}
    _fill_name   = bed_fill.name   # 'process_gas' or 'water_…'
    _bed_names   = ['uco_kernel', 'buffer_pyc', 'ipyc', 'sic', 'opyc',
                    _fill_name, 'graphite_structural']
    bed_mat_list = [_mat_by_name[n] for n in _bed_names if n in _mat_by_name]
    tallies = _build_tallies(r_ret_in, z_bed_bot, z_bed_top_val, bed_mat_list)

    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=tallies,
    )

    stats = BedStats(
        n_particles=bed['n_particles'],
        pf_achieved=bed['pf_achieved'],
        V_bulk_cm3=bed['V_bulk_cm3'],
        bed_height_cm=bed['bed_height_cm'],
        u235_mass_g=_u235_mass_g(params, bed['n_particles']),
    )
    return model, stats


def _build_tallies(
    r_ret_in: float,
    z_bed_bot: float,
    z_bed_top: float,
    bed_mat_list: list,
) -> openmc.Tallies:
    """Three tallies for step 5: flux spectrum, reaction rates, U-235 fission mesh."""

    # 1. Flux spectrum integrated over the bed — single spatial bin, 3 energy groups
    flux_mesh = openmc.RegularMesh()
    flux_mesh.lower_left  = [-r_ret_in, -r_ret_in, z_bed_bot]
    flux_mesh.upper_right = [ r_ret_in,  r_ret_in, z_bed_top]
    flux_mesh.dimension   = [1, 1, 1]

    t_flux = openmc.Tally(name='flux_spectrum_bed')
    t_flux.filters = [
        openmc.MeshFilter(flux_mesh),
        openmc.EnergyFilter(_ENERGY_BINS_EV),
    ]
    t_flux.scores = ['flux']

    # 2. Fission and absorption rates per material (all bed materials)
    t_rxn = openmc.Tally(name='reaction_rates_by_material')
    t_rxn.filters = [openmc.MaterialFilter(bed_mat_list)]
    t_rxn.scores  = ['fission', 'absorption']

    # 3. U-235 fission spatial distribution: 20×20×30 mesh over the bed volume
    fiss_mesh = openmc.RegularMesh()
    fiss_mesh.lower_left  = [-r_ret_in, -r_ret_in, z_bed_bot]
    fiss_mesh.upper_right = [ r_ret_in,  r_ret_in, z_bed_top]
    fiss_mesh.dimension   = [20, 20, 30]

    t_u235 = openmc.Tally(name='u235_fission_spatial')
    t_u235.filters  = [openmc.MeshFilter(fiss_mesh)]
    t_u235.nuclides = ['U235']
    t_u235.scores   = ['fission']

    return openmc.Tallies([t_flux, t_rxn, t_u235])


def _u235_mass_g(params, n_particles: int) -> float:
    """U-235 mass in the packed bed charge (g)."""
    kernel = _uco_kernel(params)
    nuc_dens = kernel.get_nuclide_atom_densities()   # atoms/b·cm

    N_A    = 6.02214076e23
    M_U235 = 235.044  # g/mol

    # atoms/b·cm × 1e24 cm²/b = atoms/cm³
    rho_U235_gcc = nuc_dens.get('U235', 0.0) * 1e24 * M_U235 / N_A

    r_k = params['triso']['kernel_diameter_cm'] / 2.0
    V_kernel = (4.0 / 3.0) * math.pi * r_k ** 3

    return rho_U235_gcc * V_kernel * n_particles


def export_and_run(model: openmc.Model, out_dir: Path) -> Path:
    """Export model XML files to out_dir and run OpenMC. Returns statepoint path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    model.geometry.export_to_xml( path=str(out_dir / 'geometry.xml'))
    model.materials.export_to_xml(path=str(out_dir / 'materials.xml'))
    model.settings.export_to_xml( path=str(out_dir / 'settings.xml'))
    model.tallies.export_to_xml(  path=str(out_dir / 'tallies.xml'))
    openmc.run(cwd=str(out_dir), openmc_exec=_OPENMC_EXEC)
    return out_dir / f'statepoint.{model.settings.batches}.h5'
