"""TRISO particle geometry and packed-bed universe."""

from __future__ import annotations

import hashlib
import json
import math
import types
from pathlib import Path

import numpy as np
import openmc
import openmc.model

from furnace.materials import (
    uco_kernel as _uco_kernel,
    buffer_pyc as _buffer_pyc,
    ipyc as _ipyc,
    sic as _sic,
    opyc as _opyc,
)

_N_A = 6.02214076e23
_M_C  = 12.011
_M_Si = 28.086
_M_O  = 15.999
_M_U235 = 235.044
_M_U238 = 238.051

STAGES = ('bare_kernel', 'buffered', 'ipyc', 'sic', 'full_triso')

# Sphere-center arrays are cached here; covered by the repo-root .gitignore 'cases/' entry.
_CACHE_DIR = Path(__file__).parent.parent / "cases" / ".triso_cache"


# ---------------------------------------------------------------------------
# Layer geometry helpers
# ---------------------------------------------------------------------------

def _layer_radii(params):
    """Return (r_kernel, r_buffer, r_ipyc, r_sic, r_opyc) outer radii in cm."""
    t = params['triso']
    r_k   = t['kernel_diameter_cm'] / 2
    r_buf = r_k   + t['buffer_thickness_cm']
    r_ipy = r_buf + t['ipyc_thickness_cm']
    r_sic = r_ipy + t['sic_thickness_cm']
    r_opy = r_sic + t['opyc_thickness_cm']
    return r_k, r_buf, r_ipy, r_sic, r_opy


def _shell_vol(r_outer, r_inner=0.0):
    """Volume (cm³) of a spherical shell."""
    return (4.0 / 3.0) * math.pi * (r_outer**3 - r_inner**3)


def _stage_outer_radius(stage, params):
    r_k, r_buf, r_ipy, r_sic, r_opy = _layer_radii(params)
    return {
        'bare_kernel': r_k,
        'buffered':    r_buf,
        'ipyc':        r_ipy,
        'sic':         r_sic,
        'full_triso':  r_opy,
    }[stage]


# ---------------------------------------------------------------------------
# Per-particle analytical quantities
# ---------------------------------------------------------------------------

def _u235_atom_fraction(params):
    w = params['materials']['enrichment_wt_pct'] / 100.0
    return (w / _M_U235) / (w / _M_U235 + (1.0 - w) / _M_U238)


def _kernel_formula(params):
    """Return (n_formula [units/cm³], M_formula [g/mol], M_U_avg [g/mol]) for UCO.

    UCO formula unit: U(C₀.₅O₀.₄)  →  M = M_U + 0.5·M_C + 0.4·M_O
    """
    x235 = _u235_atom_fraction(params)
    M_U = x235 * _M_U235 + (1.0 - x235) * _M_U238
    M_formula = M_U + 0.5 * _M_C + 0.4 * _M_O
    n = params['materials']['kernel_density_gcc'] * _N_A / M_formula
    return n, M_formula, M_U


def particle_mass_g(stage, params):
    """Total mass (g) of one particle at the given deposition stage."""
    m = params['materials']
    r_k, r_buf, r_ipy, r_sic, r_opy = _layer_radii(params)
    total = m['kernel_density_gcc'] * _shell_vol(r_k)
    if stage == 'bare_kernel':
        return total
    total += m['buffer_density_gcc'] * _shell_vol(r_buf, r_k)
    if stage == 'buffered':
        return total
    total += m['ipyc_density_gcc'] * _shell_vol(r_ipy, r_buf)
    if stage == 'ipyc':
        return total
    total += m['sic_density_gcc'] * _shell_vol(r_sic, r_ipy)
    if stage == 'sic':
        return total
    total += m['opyc_density_gcc'] * _shell_vol(r_opy, r_sic)
    return total  # full_triso


def u235_mass_per_particle_g(params):
    """U-235 mass (g) in one kernel; constant across all deposition stages."""
    n_formula, _, _ = _kernel_formula(params)
    x235 = _u235_atom_fraction(params)
    V_k = _shell_vol(_layer_radii(params)[0])
    return n_formula * x235 * V_k * _M_U235 / _N_A


def u_mass_per_particle_g(params):
    """Total uranium (heavy-metal) mass (g) in one kernel."""
    n_formula, _, M_U = _kernel_formula(params)
    V_k = _shell_vol(_layer_radii(params)[0])
    return n_formula * V_k * M_U / _N_A


def particle_atom_counts(stage, params):
    """Return (N_C, N_Si, N_U235) total atom counts for one particle at the given stage.

    Carbon sources: UCO kernel (0.5 C per formula unit), buffer / IPyC / OPyC (pure C),
    SiC layer (1:1 Si:C by atom count). Silicon is present only in the SiC layer.
    """
    m = params['materials']
    r_k, r_buf, r_ipy, r_sic, r_opy = _layer_radii(params)

    n_formula, _, _ = _kernel_formula(params)
    x235 = _u235_atom_fraction(params)
    V_k = _shell_vol(r_k)
    N_U235 = n_formula * x235 * V_k
    N_C    = n_formula * 0.5  * V_k   # 0.5 C per formula unit in UCO
    N_Si   = 0.0

    if stage == 'bare_kernel':
        return N_C, N_Si, N_U235

    N_C += m['buffer_density_gcc'] * _N_A / _M_C * _shell_vol(r_buf, r_k)
    if stage == 'buffered':
        return N_C, N_Si, N_U235

    N_C += m['ipyc_density_gcc'] * _N_A / _M_C * _shell_vol(r_ipy, r_buf)
    if stage == 'ipyc':
        return N_C, N_Si, N_U235

    # SiC: one C and one Si atom per formula unit; M_SiC = M_Si + M_C
    n_sic = m['sic_density_gcc'] * _N_A / (_M_Si + _M_C)
    V_sic = _shell_vol(r_sic, r_ipy)
    N_C  += n_sic * V_sic
    N_Si += n_sic * V_sic
    if stage == 'sic':
        return N_C, N_Si, N_U235

    N_C += m['opyc_density_gcc'] * _N_A / _M_C * _shell_vol(r_opy, r_sic)
    return N_C, N_Si, N_U235  # full_triso


# ---------------------------------------------------------------------------
# Particle universe builder
# ---------------------------------------------------------------------------

def particle_at_stage(stage, params):
    """Return (universe, outer_radius_cm) for the given deposition stage.

    The universe contains one concentric-sphere cell per layer present at this
    stage. All surfaces are centred at the origin; openmc.model.TRISO translates
    them to each particle's position when the bed is assembled.
    """
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {STAGES!r}, got {stage!r}")

    r_k, r_buf, r_ipy, r_sic, r_opy = _layer_radii(params)

    _cfg = {
        'bare_kernel': (
            [r_k],
            [_uco_kernel(params)],
        ),
        'buffered': (
            [r_k, r_buf],
            [_uco_kernel(params), _buffer_pyc(params)],
        ),
        'ipyc': (
            [r_k, r_buf, r_ipy],
            [_uco_kernel(params), _buffer_pyc(params), _ipyc(params)],
        ),
        'sic': (
            [r_k, r_buf, r_ipy, r_sic],
            [_uco_kernel(params), _buffer_pyc(params), _ipyc(params), _sic(params)],
        ),
        'full_triso': (
            [r_k, r_buf, r_ipy, r_sic, r_opy],
            [_uco_kernel(params), _buffer_pyc(params), _ipyc(params), _sic(params), _opyc(params)],
        ),
    }

    radii, mats = _cfg[stage]
    surfaces = [openmc.Sphere(r=r) for r in radii]

    cells = []
    for i, (surf, mat) in enumerate(zip(surfaces, mats)):
        region = -surf if i == 0 else +surfaces[i - 1] & -surf
        cells.append(openmc.Cell(fill=mat, region=region))

    return openmc.Universe(cells=cells), radii[-1]


# ---------------------------------------------------------------------------
# Packing cache
# ---------------------------------------------------------------------------

def _bbox_tuple(region):
    """Return JSON-serialisable (lower_left, upper_right) from a region's bounding box."""
    bb = region.bounding_box
    return (
        [round(float(v), 8) for v in bb.lower_left],
        [round(float(v), 8) for v in bb.upper_right],
    )


def _make_cache_key(region, packing_fraction, outer_radius, seed):
    ll, ur = _bbox_tuple(region)
    payload = json.dumps(
        {'ll': ll, 'ur': ur,
         'pf': round(packing_fraction, 8),
         'r':  round(outer_radius, 8),
         'seed': seed},
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# Packed-bed builder
# ---------------------------------------------------------------------------

def pack_bed(region, packing_fraction, outer_radius, fill_universe, seed, params):
    """Generate a random packed bed of TRISO particles within *region*.

    Sphere centres are written to cases/.triso_cache/<hash>.npz on first call
    and reloaded on subsequent calls with the same (region bbox, pf, r, seed).
    The cache key is derived from the region's axis-aligned bounding box; for
    non-convex regions that share a bounding box with a different shape, pass
    distinct seeds or verify manually that the cached geometry is appropriate.

    An operating fluidised bed runs at packing_fraction_fluidized (~0.333);
    a gravity-settled collapsed bed uses packing_fraction_static (~0.50), well
    below the random close-packing limit of 0.64.
    """
    max_pf = params['model']['max_packing_fraction']
    if packing_fraction > max_pf:
        raise ValueError(
            f"packing_fraction {packing_fraction:.4f} exceeds "
            f"max_packing_fraction {max_pf:.4f} (random close-packing limit). "
            "A collapsed or vibration-settled bed may approach this limit; "
            "an operating fluidised bed is far below it."
        )

    key = _make_cache_key(region, packing_fraction, outer_radius, seed)
    cache_path = _CACHE_DIR / f"{key}.npz"

    if cache_path.exists():
        centers = np.load(cache_path)['centers']
    else:
        # Guard: skip regions too small to fit even one sphere (partial-slab
        # edge case where num_spheres=0 causes a ZeroDivisionError inside
        # pack_spheres' RSP mesh initialisation).
        V_sphere = (4.0 / 3.0) * math.pi * outer_radius**3
        bb = region.bounding_box
        if int(packing_fraction * float(bb.volume) / V_sphere) == 0:
            return []

        # Thin-slab / narrow-cylinder guard: FBP degenerates in two regimes:
        #   (a) near-monolayer geometry: h < 2.5d (partial last slab at h ≈ 1.5d).
        #   (b) narrow confined cylinder: r < 15 particle radii (slab 0 at
        #       r=14.12 radii — near-wall ordering pushes achievable pf near 0.50).
        # Fix: extend the region to 5d and reduce pf proportionally so that
        # pf_ext = pf × z_extent / (5d) ≤ 0.25 across the full trigger range
        # (worst case: z=2.5d, pf=0.50 → pf_ext=0.25), forcing RSP (pf<0.30).
        # Guard only fires when min_h > z_extent; for tall slabs (z_extent ≥ min_h)
        # FBP converges normally and the extension path is skipped.
        # Yield ≈ z_extent / min_h (~30 % for the 95 g / n_slabs=32 partial slab).
        _TRIGGER_DIAMETERS = 2.5     # height trigger threshold
        _NARROW_R_RADII = 15.0       # radius trigger: r < 15 particle radii
        _MIN_H_DIAMETERS = 5.0       # RSP extension height (both triggers)

        z_bot = float(bb.lower_left[2])
        z_top = float(bb.upper_right[2])
        z_extent = z_top - z_bot
        r_slab = float(bb.upper_right[0])  # cylinder radius from bounding box

        narrow_cylinder = r_slab / outer_radius < _NARROW_R_RADII
        trigger_h = _TRIGGER_DIAMETERS * 2.0 * outer_radius
        min_h = _MIN_H_DIAMETERS * 2.0 * outer_radius

        if (z_extent < trigger_h or narrow_cylinder) and min_h > z_extent:
            z_top_ext = z_bot + min_h
            ext_region = (
                -openmc.ZCylinder(r=r_slab)
                & +openmc.ZPlane(z0=z_bot)
                & -openmc.ZPlane(z0=z_top_ext)
            )
            pf_ext = packing_fraction * z_extent / min_h
            all_centers = openmc.model.pack_spheres(
                radius=outer_radius,
                region=ext_region,
                pf=pf_ext,
                seed=seed,
            )
            # Keep only centers whose sphere lies fully within the original slab.
            mask = (all_centers[:, 2] - outer_radius >= z_bot) & \
                   (all_centers[:, 2] + outer_radius <= z_top)
            centers = all_centers[mask]
        else:
            centers = openmc.model.pack_spheres(
                radius=outer_radius,
                region=region,
                pf=packing_fraction,
                seed=seed,
            )

        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(cache_path, centers=centers)

    return [openmc.model.TRISO(outer_radius, fill_universe, tuple(c)) for c in centers]


# ---------------------------------------------------------------------------
# Lattice wrapper
# ---------------------------------------------------------------------------

def lattice_bed(trisos, lower_left, upper_right, background_material, pitch_target):
    """Wrap TRISO objects in an optimised lookup lattice for efficient tracking.

    *pitch_target* is the desired cubic cell side length (cm). The actual pitch
    is adjusted to fit an integer number of cells between lower_left and
    upper_right in each dimension.

    Recommended pitch_target = 4 × outer_radius (= 2 × particle diameter):
    each cubic cell holds at most one particle even at the densest packing,
    while the corner dead-space is small. At full TRISO outer radius ≈ 0.046 cm
    this gives pitch ≈ 0.185 cm (~2 cells/mm).
    """
    ll = np.asarray(lower_left,  dtype=float)
    ur = np.asarray(upper_right, dtype=float)
    span  = ur - ll
    shape = tuple(max(1, int(math.ceil(span[i] / pitch_target))) for i in range(3))
    pitch = tuple(float(span[i] / shape[i]) for i in range(3))

    return openmc.model.create_triso_lattice(
        trisos,
        lower_left=tuple(ll),
        pitch=pitch,
        shape=shape,
        background=background_material,
    )


# ---------------------------------------------------------------------------
# Bed statistics
# ---------------------------------------------------------------------------

def bed_stats(trisos, region_volume, stage, params):
    """Return a dict of bed-level and per-particle statistics.

    Parameters
    ----------
    trisos        : list of openmc.model.TRISO returned by pack_bed
    region_volume : bulk volume of the bed region (cm³)
    stage         : deposition stage string; one of STAGES
    params        : frozen params dict from load_params()

    Returns
    -------
    dict with keys:
        stage, n_particles, packing_fraction, outer_radius_cm,
        particle_mass_g, hm_mass_g, u235_mass_g, u235_mass_fraction,
        smeared_density_gcc, c_per_u235, si_per_u235, c_mass_g
    """
    outer_r = _stage_outer_radius(stage, params)
    n       = len(trisos)
    pf_ach  = n * _shell_vol(outer_r) / region_volume

    p_mass  = particle_mass_g(stage, params)
    u235_m  = u235_mass_per_particle_g(params)
    hm_m    = u_mass_per_particle_g(params)
    N_C, N_Si, N_U235 = particle_atom_counts(stage, params)

    return {
        'stage':               stage,
        'n_particles':         n,
        'packing_fraction':    pf_ach,
        'outer_radius_cm':     outer_r,
        'particle_mass_g':     p_mass,
        'hm_mass_g':           n * hm_m,
        'u235_mass_g':         n * u235_m,
        'u235_mass_fraction':  u235_m / p_mass,
        'smeared_density_gcc': n * p_mass / region_volume,
        'c_per_u235':          N_C  / N_U235,
        'si_per_u235':         N_Si / N_U235,
        'c_mass_g':            n * N_C * _M_C / _N_A,
    }


def print_bed_stats(stats):
    """Pretty-print the output of bed_stats()."""
    w = 52
    print(f"\n{'='*w}")
    print(f"  Bed stats — stage: {stats['stage']}")
    print(f"{'='*w}")
    print(f"  Particles               : {stats['n_particles']}")
    print(f"  Achieved PF             : {stats['packing_fraction']:.4f}")
    print(f"  Outer radius (cm)       : {stats['outer_radius_cm']:.6f}")
    print(f"  Particle mass (g)       : {stats['particle_mass_g']:.4e}")
    print(f"  U-235 mass fraction     : {stats['u235_mass_fraction']:.4f}")
    print(f"  C / U-235  (atoms)      : {stats['c_per_u235']:.2f}")
    print(f"  Si / U-235 (atoms)      : {stats['si_per_u235']:.2f}")
    print(f"  Total HM mass (g)       : {stats['hm_mass_g']:.4e}")
    print(f"  Total U-235 mass (g)    : {stats['u235_mass_g']:.4e}")
    print(f"  Total carbon mass (g)   : {stats['c_mass_g']:.4e}")
    print(f"  Smeared density (g/cm³) : {stats['smeared_density_gcc']:.4e}")


# ---------------------------------------------------------------------------
# Analytical stage-progression table (no transport)
# ---------------------------------------------------------------------------

def stage_progression_table(params, n_particles_fixed, region_volume_fixed, pf_fixed):
    """Print stage-by-stage analytical stats under both fixed-mass and fixed-volume framings.

    Fixed-mass framing: the number of kernels is fixed (the physical charge).
    As layers are deposited each particle grows, so the bed volume expands if
    the packing fraction is held constant; in a rigid vessel the packing fraction
    would fall instead.  The uranium inventory is the same at every stage.

    Fixed-volume framing: bed volume and packing fraction are fixed.  Fewer
    (larger) particles fit as layers accumulate, so the uranium inventory
    decreases at each stage.

    The fixed-mass framing is more meaningful for NCS: a real coating run
    starts with a fixed charge of kernels and adds material to each particle.
    The question is which stage yields the highest k-eff, not how many particles
    are present.  The fixed-volume framing confounds the effect of particle
    geometry with a change in total fissile inventory, making trends harder to
    interpret physically.
    """
    hdr_fm = (f"{'Stage':<14} {'r_out(cm)':>10} {'V_bed(cm³)':>11} "
              f"{'C/U235':>8} {'Si/U235':>8} {'HM(g)':>8} {'U235(g)':>10}")
    print("\n--- FIXED MASS (N kernels fixed, V_bed expands at fixed PF) ---")
    print(hdr_fm)
    print('-' * len(hdr_fm))
    for stage in STAGES:
        outer_r  = _stage_outer_radius(stage, params)
        v_bed    = n_particles_fixed * _shell_vol(outer_r) / pf_fixed
        hm_m     = u_mass_per_particle_g(params)
        u235_m   = u235_mass_per_particle_g(params)
        N_C, N_Si, N_U235 = particle_atom_counts(stage, params)
        print(f"{stage:<14} {outer_r:>10.5f} {v_bed:>11.3f} "
              f"{N_C/N_U235:>8.2f} {N_Si/N_U235:>8.2f} "
              f"{n_particles_fixed*hm_m:>8.4f} {n_particles_fixed*u235_m:>10.4e}")

    hdr_fv = (f"{'Stage':<14} {'r_out(cm)':>10} {'N_particles':>12} "
              f"{'C/U235':>8} {'Si/U235':>8} {'HM(g)':>8} {'U235(g)':>10}")
    print("\n--- FIXED VOLUME (V_bed and PF fixed, N decreases as particle grows) ---")
    print(hdr_fv)
    print('-' * len(hdr_fv))
    for stage in STAGES:
        outer_r  = _stage_outer_radius(stage, params)
        n        = int(pf_fixed * region_volume_fixed / _shell_vol(outer_r))
        hm_m     = u_mass_per_particle_g(params)
        u235_m   = u235_mass_per_particle_g(params)
        N_C, N_Si, N_U235 = particle_atom_counts(stage, params)
        print(f"{stage:<14} {outer_r:>10.5f} {n:>12} "
              f"{N_C/N_U235:>8.2f} {N_Si/N_U235:>8.2f} "
              f"{n*hm_m:>8.4f} {n*u235_m:>10.4e}")


# ---------------------------------------------------------------------------
# Diagnostic entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    from furnace.params import check_env, load_params

    check_env()
    params = load_params()

    # Build a small bare-kernel test bed inside a 0.25 cm sphere.
    # At r_kernel ≈ 0.02125 cm and PF = 0.39 this gives ~630 particles —
    # enough for a meaningful mass check while pack_spheres finishes quickly.
    stage = 'bare_kernel'
    fill_univ, outer_r = particle_at_stage(stage, params)

    test_r_cm = 0.25
    container  = openmc.Sphere(r=test_r_cm)
    region     = -container
    V_region   = _shell_vol(test_r_cm)

    pf_target  = float(params['model']['packing_fraction_fluidized'])
    seed       = int(params['model']['seed'])

    print(f"\nPacking {stage} (r_particle = {outer_r:.5f} cm) into sphere "
          f"r = {test_r_cm} cm, PF target = {pf_target} ...")
    trisos = pack_bed(region, pf_target, outer_r, fill_univ, seed, params)

    stats = bed_stats(trisos, V_region, stage, params)
    print_bed_stats(stats)

    # --- Mass verification ---
    # Analytic estimate uses target PF; deviation reflects how closely
    # pack_spheres achieved the requested packing fraction.
    m_analytic = pf_target * V_region * params['materials']['kernel_density_gcc']
    m_computed = stats['n_particles'] * particle_mass_g(stage, params)
    pct_diff   = 100.0 * (m_computed - m_analytic) / m_analytic
    print(f"\nMass check (bare_kernel):")
    print(f"  Analytic  PF_target × V × ρ_kernel : {m_analytic:.4e} g")
    print(f"  Computed  N × m_particle             : {m_computed:.4e} g")
    print(f"  % difference                         : {pct_diff:+.2f}%")
    print(f"  (Achieved PF {stats['packing_fraction']:.4f} vs target {pf_target:.4f})")

    # --- Stage progression (analytical, no transport) ---
    # N_fixed: number of bare kernels that fill the nominal static bed at PF_static.
    V_nominal  = params['dimensions']['bed']['static_bulk_cc']
    pf_static  = float(params['model']['packing_fraction_static'])
    r_k        = _layer_radii(params)[0]
    N_fixed    = int(pf_static * V_nominal / _shell_vol(r_k))

    print(f"\nNominal bed: V = {V_nominal} cm³, PF_static = {pf_static}, "
          f"N_fixed (bare kernels) = {N_fixed}")
    stage_progression_table(params, N_fixed, V_nominal, pf_static)
