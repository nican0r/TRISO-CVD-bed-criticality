"""Furnace geometry: bed, cone, retort cylinder, nozzle throat. Steps 3-4.

Coordinate system: z=0 at the nozzle-to-cone junction plane, z increases upward.
  Throat:           z in [-h_throat, 0]             r = r_throat (graphite, no particles)
  Cone (staircase): z in [0, z_cone_top]            r from r_throat to r_retort
  Retort cylinder:  z in [z_cone_top, z_retort_top]  r = r_retort
"""

from __future__ import annotations

import hashlib
import json
import math
import warnings

import numpy as np
import openmc
import openmc.model

from furnace.triso import (
    particle_at_stage,
    pack_bed,
    lattice_bed,
    tiled_bed,
    _stage_outer_radius,
    _shell_vol,
    particle_mass_g,
    _CACHE_DIR,
)


# ---------------------------------------------------------------------------
# Primary surface constructors
# ---------------------------------------------------------------------------

def retort_inner_cylinder(params):
    """ZCylinder at retort inner radius (r = r_retort)."""
    r_retort = params['dimensions']['retort']['id_cm'] / 2.0
    return openmc.ZCylinder(r=r_retort)


def cone_wall(params):
    """Return (ZCone, z_plane_bottom, z_plane_top) for the structural graphite cone.

    60° included angle = 30° half-angle → r2 = tan²(30°) = 1/3.
    Algebra: ZCone equation x²+y² = r2*(z-z0)²
    Apex z0 = -r_throat/tan(30°) (below z=0 origin).
    Verification: at z=0, r = tan(30°)*|z_apex| = r_throat ✓
                  at z=z_cone_top, r = tan(30°)*(z_cone_top - z_apex) = r_retort ✓
    Truncated by ZPlanes at z=0 (throat junction) and z=z_cone_top (retort junction).
    """
    dim = params['dimensions']
    r_throat = dim['nozzle']['throat_diameter_cm'] / 2.0
    z_cone_top = dim['cone']['vertical_drop_cm']

    half_angle_rad = math.radians(dim['cone']['included_angle_deg'] / 2.0)
    tan_theta = math.tan(half_angle_rad)
    r2 = tan_theta ** 2          # = 1/3 for 30° half-angle

    z_apex = -r_throat / tan_theta   # apex is below z=0

    cone = openmc.ZCone(x0=0.0, y0=0.0, z0=z_apex, r2=r2)
    z_bot = openmc.ZPlane(z0=0.0)
    z_top = openmc.ZPlane(z0=z_cone_top)
    return cone, z_bot, z_top


def nozzle_throat(params):
    """Return (ZCylinder, z_plane_bottom, z_plane_top) for the throat cylinder.

    Throat: r = r_throat, axial extent z in [-h_throat, 0].
    Filled with graphite_structural (no particles ever enter this region).
    """
    dim = params['dimensions']
    r_throat = dim['nozzle']['throat_diameter_cm'] / 2.0
    h_throat = dim['nozzle']['nozzle_height_cm']

    cyl = openmc.ZCylinder(r=r_throat)
    z_bot = openmc.ZPlane(z0=-h_throat)
    z_top = openmc.ZPlane(z0=0.0)
    return cyl, z_bot, z_top


# ---------------------------------------------------------------------------
# Volume helpers
# ---------------------------------------------------------------------------

def frustum_volume(r_top, r_bottom, half_angle_deg=30.0):
    """Analytic truncated-cone volume: V = (π/3)*(r_top³ - r_bottom³)/tan(θ).

    Derivation: V = ∫π r(z)² dz; let u=r(z)=tan(θ)(z-z_apex), dz=du/tan(θ)
    → V = (π/tan θ) ∫_{r_bot}^{r_top} u² du = (π/3)(r_top³-r_bot³)/tan θ
    """
    tan_theta = math.tan(math.radians(half_angle_deg))
    return (math.pi / 3.0) * (r_top**3 - r_bottom**3) / tan_theta


# ---------------------------------------------------------------------------
# Staircase approximation
# ---------------------------------------------------------------------------

def staircase_bed(params, n_slabs):
    """Inscribed cylindrical staircase approximation of the cone region.

    Divides the full cone vertical extent (0 to z_cone_top) into n_slabs
    equal-height slabs. Each slab radius = cone radius at the slab's bottom
    (inscribed fit: guaranteed no particle protrudes through cone wall; a
    midpoint or circumscribed fit would produce overlapping-cell geometry).

    Returns:
        slabs: list of dicts {z_bot, z_top, radius, volume}
        total_volume_cm3: sum of slab volumes (strictly < frustum_volume)
        volume_error_frac: (frustum_volume - total_volume) / frustum_volume
    """
    dim = params['dimensions']
    r_throat = dim['nozzle']['throat_diameter_cm'] / 2.0
    z_cone_top = dim['cone']['vertical_drop_cm']

    half_angle_rad = math.radians(dim['cone']['included_angle_deg'] / 2.0)
    tan_theta = math.tan(half_angle_rad)
    z_apex = -r_throat / tan_theta

    dz = z_cone_top / n_slabs

    slabs = []
    for i in range(n_slabs):
        z_bot = i * dz
        z_top = (i + 1) * dz
        # Inscribed: radius at the bottom of the slab
        r = tan_theta * (z_bot - z_apex)
        vol = math.pi * r**2 * dz
        slabs.append({'z_bot': z_bot, 'z_top': z_top, 'radius': r, 'volume': vol})

    total_vol = sum(s['volume'] for s in slabs)

    r_retort = dim['retort']['id_cm'] / 2.0
    v_frustum = frustum_volume(r_retort, r_throat, half_angle_deg=dim['cone']['included_angle_deg'] / 2.0)
    vol_error_frac = (v_frustum - total_vol) / v_frustum

    return slabs, total_vol, vol_error_frac


# ---------------------------------------------------------------------------
# Per-particle effective density
# ---------------------------------------------------------------------------

def _particle_effective_density(stage, params):
    """ρ_eff = m_particle / V_particle_outer_sphere (g/cm³)."""
    r_out = _stage_outer_radius(stage, params)
    m = particle_mass_g(stage, params)
    return m / _shell_vol(r_out)


def _flood_fill(z_bot, wet, dry, z_flood):
    """Select wet or dry material based on whether z_bot is below the flood level."""
    if z_flood is None or dry is None:
        return wet
    return wet if z_bot < z_flood else dry


# ---------------------------------------------------------------------------
# Bed region assembly
# ---------------------------------------------------------------------------

def bed_region(params, state, stage, n_slabs, background_material, charge_mass_g=None, seed=None, z_flood=None, dry_material=None, *, use_tiled_bed=False, tile_size_cm=0.5):
    """Assemble the particle bed for the given state and deposition stage.

    Parameters
    ----------
    params              : frozen params dict from load_params()
    state               : 'fluidized' or 'collapsed'
    stage               : deposition stage (one of furnace.triso.STAGES)
    n_slabs             : number of equal-height staircase slabs in the cone
    background_material : openmc.Material for interstitial space (e.g. process_gas or water)
    charge_mass_g       : total particle charge mass (g). If None, reads from
                          params['dimensions']['bed']['charge_mass_g'].
                          Pass explicitly to sweep charge mass without mutating params.
    z_flood             : absolute z (cm) of the bottom-up flood level. Cells whose
                          bottom face is below z_flood receive background_material (wet);
                          cells above receive dry_material. None → uniform background.
    dry_material        : openmc.Material used above z_flood (typically process_gas).
                          Ignored when z_flood is None.

    Packing fractions:
        collapsed  : params['model']['packing_fraction_static']   (~0.50)
        fluidized  : pf_static / bed_expansion_ratio               (~0.333)

    Both states use the same inscribed staircase geometry: cone slabs filled from
    bottom to top at the state's packing fraction, then overflow into the retort
    cylinder. For the fluidized state the overflow height is checked against
    fluidized_height_cm (available cylinder zone). The cone interior is never
    left empty of particles unless the charge mass is too small to reach the top.

    Volume computation (mass-conserving):
        rho_eff = m_particle / V_particle   (stage-dependent)
        V_solid = charge_mass_g / rho_eff   (conserved solid volume)
        V_bulk  = V_solid / pf              (bulk volume for this state)

    Returns
    -------
    dict with keys:
        'cells'          : list of openmc.Cell for the bed
        'trisos'         : list of openmc.model.TRISO (all regions combined)
        'pf_achieved'    : achieved packing fraction (n_particles * V_particle / V_bulk)
        'V_bulk_cm3'     : bulk volume used
        'V_solid_cm3'    : solid particle volume
        'bed_height_cm'  : absolute z of bed top
        'z_bed_bot_cm'   : absolute z of bed bottom (always 0.0 — cone base)
        'n_particles'    : total number of particles packed
        'background_mat' : background_material (passed through for model assembly)
    """
    if state not in ('fluidized', 'collapsed'):
        raise ValueError(f"state must be 'fluidized' or 'collapsed', got {state!r}")

    dim = params['dimensions']
    mdl = params['model']

    if charge_mass_g is None:
        charge_mass_g = float(dim['bed']['charge_mass_g'])

    pf_static = float(mdl['packing_fraction_static'])
    bed_expansion_ratio = float(mdl['bed_expansion_ratio'])
    pf_fluidized = pf_static / bed_expansion_ratio
    max_pf = float(mdl['max_packing_fraction'])

    pf = pf_static if state == 'collapsed' else pf_fluidized

    r_retort = dim['retort']['id_cm'] / 2.0
    z_cone_top = dim['cone']['vertical_drop_cm']
    fluidized_height_cm = dim['bed']['fluidized_height_cm']

    rho_eff = _particle_effective_density(stage, params)
    V_solid = charge_mass_g / rho_eff
    V_bulk = V_solid / pf

    outer_r = _stage_outer_radius(stage, params)
    fill_univ, _ = particle_at_stage(stage, params)
    # pitch_target: aim for ~1 particle per lattice cell, but cap at 0.20 cm.
    # 4×r for bare_kernel (r≈0.021 cm) → 0.085 cm → 275 k-cell overflow lattice → 211 MB XML.
    # Cap at 0.20 cm keeps the overflow lattice ≈ 25×25×34 ≈ 21 k cells → ~16 MB XML,
    # with ~10-40 particles/cell — acceptable for tracking in a deeply subcritical geometry.
    pitch_target = max(4.0 * outer_r, 0.20)

    base_seed = int(seed) if seed is not None else int(mdl['seed'])

    all_cells = []
    all_trisos = []
    total_particles = 0

    # Cone surface used for void cells in both states.
    # Same parameters as cone_wall(); separate surface object is fine in OpenMC.
    dim_c = params['dimensions']
    _r_throat = dim_c['nozzle']['throat_diameter_cm'] / 2.0
    _half_rad = math.radians(dim_c['cone']['included_angle_deg'] / 2.0)
    _tan_th = math.tan(_half_rad)
    _z_apex = -_r_throat / _tan_th
    _cone_surf = openmc.ZCone(z0=_z_apex, r2=_tan_th**2)
    _zp_cone_bot = openmc.ZPlane(z0=0.0)
    _zp_cone_top = openmc.ZPlane(z0=z_cone_top)

    # --- staircase fill: cone bottom → top, then overflow into retort cylinder ---
    # Both states use the same inscribed staircase geometry; the only difference is pf.
    slabs, staircase_vol, vol_error_frac = staircase_bed(params, n_slabs)

    remaining_bulk = V_bulk
    slab_index = 0
    z_last_fill_top = 0.0   # tracks the bed top within the cone axial extent
    cone_void_cells = []    # annular voids + unfilled cone interior

    for slab in slabs:
        if remaining_bulk <= 1e-12:
            break

        r_slab = slab['radius']
        z_bot = slab['z_bot']
        z_top = slab['z_top']
        slab_vol = slab['volume']

        if remaining_bulk >= slab_vol:
            used_vol = slab_vol
            z_fill_top = z_top
        else:
            # Partial slab: solve π*r²*h = remaining_bulk
            h_partial = remaining_bulk / (math.pi * r_slab**2)
            z_fill_top = z_bot + h_partial
            used_vol = remaining_bulk

        cyl = openmc.ZCylinder(r=r_slab)
        zp_bot_s = openmc.ZPlane(z0=z_bot)
        zp_top_s = openmc.ZPlane(z0=z_fill_top)

        # Bed cell: inscribed cylinder packed with particles
        slab_fill = _flood_fill(z_bot, background_material, dry_material, z_flood)
        region = -cyl & +zp_bot_s & -zp_top_s
        seed = base_seed + slab_index
        if use_tiled_bed:
            lattice, tile_pf_ach, _n_per_tile = tiled_bed(
                lower_left=(-r_slab, -r_slab, z_bot),
                upper_right=(r_slab, r_slab, z_fill_top),
                packing_fraction=pf,
                outer_radius=outer_r,
                fill_universe=fill_univ,
                background_material=slab_fill,
                seed=seed,
                params=params,
                tile_size_cm=tile_size_cm,
            )
            all_cells.append(openmc.Cell(fill=lattice, region=region))
            v_slab_cyl = math.pi * r_slab ** 2 * (z_fill_top - z_bot)
            total_particles += int(round(tile_pf_ach * v_slab_cyl / _shell_vol(outer_r)))
        else:
            trisos = pack_bed(region, pf, outer_r, fill_univ, seed, params)
            lattice = lattice_bed(
                trisos,
                lower_left=(-r_slab, -r_slab, z_bot),
                upper_right=(r_slab, r_slab, z_fill_top),
                background_material=slab_fill,
                pitch_target=pitch_target,
            )
            all_cells.append(openmc.Cell(fill=lattice, region=region))
            all_trisos.extend(trisos)
            total_particles += len(trisos)

        # Annular void: inside cone, outside inscribed cylinder, at this slab height.
        # At z_bot the inscribed radius equals the cone radius (inscribed fit), so the
        # annular void has zero width there and grows toward z_fill_top. This cell is
        # required so every point inside the cone is covered by exactly one cell.
        annular_region = +cyl & -_cone_surf & +zp_bot_s & -zp_top_s
        cone_void_cells.append(
            openmc.Cell(fill=slab_fill, region=annular_region)
        )

        z_last_fill_top = z_fill_top
        remaining_bulk -= used_vol
        slab_index += 1

    # Unfilled cone interior above bed (if the charge doesn't fill the full cone).
    # When z_flood falls inside this region, split exactly at z_flood so the flood
    # level is accurate (avoids the slab-height over-approximation for this void).
    if z_last_fill_top < z_cone_top - 1e-10:
        zp_bed_last = openmc.ZPlane(z0=z_last_fill_top)
        if (z_flood is not None and dry_material is not None
                and z_last_fill_top < z_flood < z_cone_top):
            zp_flood_in_cone = openmc.ZPlane(z0=z_flood)
            cone_void_cells.append(
                openmc.Cell(fill=background_material,
                            region=-_cone_surf & +zp_bed_last & -zp_flood_in_cone)
            )
            cone_void_cells.append(
                openmc.Cell(fill=dry_material,
                            region=-_cone_surf & +zp_flood_in_cone & -_zp_cone_top)
            )
        else:
            void_fill = _flood_fill(z_last_fill_top, background_material, dry_material, z_flood)
            cone_void_cells.append(
                openmc.Cell(
                    fill=void_fill,
                    region=-_cone_surf & +zp_bed_last & -_zp_cone_top,
                )
            )

    # Overflow into retort cylinder above cone if charge exceeds staircase volume
    if remaining_bulk > 1e-10:
        h_overflow = remaining_bulk / (math.pi * r_retort**2)

        if state == 'fluidized' and h_overflow > fluidized_height_cm:
            raise ValueError(
                f"Fluidized bed overflow height {h_overflow:.3f} cm exceeds available zone "
                f"{fluidized_height_cm:.3f} cm. Reduce charge_mass_g or check bed_expansion_ratio."
            )

        if h_overflow < 2.0 * outer_r:
            # Overflow thinner than one particle diameter — no TRISO can fit.
            # This arises when V_bulk barely exceeds V_staircase due to the
            # inscribed-staircase volume shortfall (~4.7% at n_slabs=32).
            # Absorb the residual bulk into the last cone slab: the gas-above-bed
            # cell in model.py starts at max(bed_height, z_cone_top), so setting
            # bed_height_cm = z_last_fill_top leaves no geometry gap.
            bed_height_cm = z_last_fill_top
        else:
            z_ov_bot = z_cone_top
            z_ov_top = z_cone_top + h_overflow

            cyl_ret = openmc.ZCylinder(r=r_retort)
            zp_ov_bot = openmc.ZPlane(z0=z_ov_bot)
            zp_ov_top = openmc.ZPlane(z0=z_ov_top)
            region_ov = -cyl_ret & +zp_ov_bot & -zp_ov_top

            ov_fill = _flood_fill(z_ov_bot, background_material, dry_material, z_flood)
            seed = base_seed + slab_index
            if use_tiled_bed:
                lattice_ov, tile_pf_ach, _n_per_tile = tiled_bed(
                    lower_left=(-r_retort, -r_retort, z_ov_bot),
                    upper_right=(r_retort, r_retort, z_ov_top),
                    packing_fraction=pf,
                    outer_radius=outer_r,
                    fill_universe=fill_univ,
                    background_material=ov_fill,
                    seed=seed,
                    params=params,
                    tile_size_cm=tile_size_cm,
                )
                all_cells.append(openmc.Cell(fill=lattice_ov, region=region_ov))
                v_ov_cyl = math.pi * r_retort ** 2 * (z_ov_top - z_ov_bot)
                total_particles += int(round(tile_pf_ach * v_ov_cyl / _shell_vol(outer_r)))
            else:
                trisos_ov = pack_bed(region_ov, pf, outer_r, fill_univ, seed, params)
                lattice_ov = lattice_bed(
                    trisos_ov,
                    lower_left=(-r_retort, -r_retort, z_ov_bot),
                    upper_right=(r_retort, r_retort, z_ov_top),
                    background_material=ov_fill,
                    pitch_target=pitch_target,
                )
                all_cells.append(openmc.Cell(fill=lattice_ov, region=region_ov))
                all_trisos.extend(trisos_ov)
                total_particles += len(trisos_ov)

            bed_height_cm = z_ov_top
    else:
        bed_height_cm = z_last_fill_top

    # In tiled mode `all_trisos` is intentionally empty (unique cells live in the
    # single tile universe); fall back to `total_particles`, which the tiled path
    # populates analytically from the tile pf × region volume.
    n_actual = len(all_trisos) if all_trisos else total_particles
    V_actual_solid = n_actual * _shell_vol(outer_r)
    pf_achieved = V_actual_solid / V_bulk

    assert pf_achieved <= max_pf, (
        f"Achieved PF {pf_achieved:.4f} exceeds max_packing_fraction {max_pf:.4f}"
    )

    return {
        'cells': all_cells,
        'cone_void_cells': cone_void_cells,
        'trisos': all_trisos,
        'pf_achieved': pf_achieved,
        'V_bulk_cm3': V_bulk,
        'V_solid_cm3': V_solid,
        'bed_height_cm': bed_height_cm,   # absolute z of bed top
        'z_bed_bot_cm': 0.0,              # bed always starts at cone base (z=0)
        'n_particles': total_particles,
        'background_mat': background_material,
    }


# ---------------------------------------------------------------------------
# Exact-cone reference geometry (validation only)
# ---------------------------------------------------------------------------

def _exact_cone_cache_key(params, target_n: int, outer_radius: float, base_seed: int) -> str:
    dim = params['dimensions']
    payload = json.dumps({
        'type': 'exact_cone_v1',
        'r_throat_cm': round(dim['nozzle']['throat_diameter_cm'] / 2.0, 8),
        'r_retort_cm': round(dim['retort']['id_cm'] / 2.0, 8),
        'z_cone_top_cm': round(dim['cone']['vertical_drop_cm'], 8),
        'half_angle_deg': round(dim['cone']['included_angle_deg'] / 2.0, 8),
        'target_n': int(target_n),
        'r_outer': round(outer_radius, 8),
        'seed': int(base_seed),
    }, sort_keys=True)
    return 'ec_' + hashlib.sha256(payload.encode()).hexdigest()[:16]


def _cone_acceptance_mask(centers: np.ndarray, outer_radius: float,
                          z_apex: float, tan_theta: float, z_cone_top: float) -> np.ndarray:
    """Boolean mask: True if the sphere body clears the cone wall and both planes.

    Perpendicular clearance to the cone surface:
        d_perp = cos(θ) × (r_cone(z) − ρ)   where r_cone(z) = tan(θ)·(z − z_apex)
    A sphere of radius r clears the wall when d_perp ≥ r.
    It clears the floor (z=0) when z ≥ r, and the ceiling (z=z_cone_top) when z ≤ z_cone_top − r.
    """
    z = centers[:, 2]
    rho = np.hypot(centers[:, 0], centers[:, 1])
    cos_theta = 1.0 / np.sqrt(1.0 + tan_theta ** 2)
    d_perp = cos_theta * (tan_theta * (z - z_apex) - rho)
    return (
        (d_perp >= outer_radius)
        & (z >= outer_radius)
        & (z <= z_cone_top - outer_radius)
    )


def exact_cone_bed(params, state: str, stage: str, background_material,
                   charge_mass_g: float | None = None, seed: int | None = None,
                   z_flood: float | None = None, dry_material=None,
                   *, use_tiled_bed: bool = False, tile_size_cm: float = 0.5) -> dict:
    """Pack the true frustum by rejection sampling — validation-only reference geometry.

    Does not use the staircase approximation.  Sphere centres are placed in the
    bounding cylinder (r=r_retort, z in [0, z_cone_top]) by pack_spheres, then
    any sphere whose body intersects the cone surface or the floor/ceiling planes
    is discarded.  pf_trial in the bounding cylinder is iterated until the
    surviving count matches the charge mass to within 0.5 %.

    The accepted centres are cached under cases/.triso_cache/ with an 'ec_' key
    prefix distinct from staircase-packing cache keys.

    Returns the same dict structure as bed_region(), with cone_void_cells=[].
    """
    if state not in ('fluidized', 'collapsed'):
        raise ValueError(f"state must be 'fluidized' or 'collapsed', got {state!r}")

    dim = params['dimensions']
    mdl = params['model']

    if charge_mass_g is None:
        charge_mass_g = float(dim['bed']['charge_mass_g'])

    pf_static = float(mdl['packing_fraction_static'])
    pf = pf_static if state == 'collapsed' else pf_static / float(mdl['bed_expansion_ratio'])
    max_pf = float(mdl['max_packing_fraction'])

    r_throat = dim['nozzle']['throat_diameter_cm'] / 2.0
    r_retort = dim['retort']['id_cm'] / 2.0
    z_cone_top = dim['cone']['vertical_drop_cm']
    half_angle_rad = math.radians(dim['cone']['included_angle_deg'] / 2.0)
    tan_theta = math.tan(half_angle_rad)
    z_apex = -r_throat / tan_theta

    outer_r = _stage_outer_radius(stage, params)
    rho_eff = _particle_effective_density(stage, params)
    V_solid = charge_mass_g / rho_eff
    V_bulk = V_solid / pf

    base_seed = int(seed) if seed is not None else int(mdl['seed'])
    fill_univ, _ = particle_at_stage(stage, params)
    pitch_target = max(4.0 * outer_r, 0.20)

    target_n = int(round(charge_mass_g / particle_mass_g(stage, params)))
    v_sphere = _shell_vol(outer_r)
    v_bounding_cyl = math.pi * r_retort ** 2 * z_cone_top
    v_frustum = frustum_volume(r_retort, r_throat,
                               half_angle_deg=dim['cone']['included_angle_deg'] / 2.0)

    cache_key = _exact_cone_cache_key(params, target_n, outer_r, base_seed)
    cache_path = _CACHE_DIR / f'{cache_key}.npz'

    if cache_path.exists():
        accepted = np.load(cache_path)['centers']
    else:
        bound_cyl  = openmc.ZCylinder(r=r_retort)
        bound_zbot = openmc.ZPlane(z0=0.0)
        bound_ztop = openmc.ZPlane(z0=z_cone_top)
        bounding_region = -bound_cyl & +bound_zbot & -bound_ztop

        # Initial pf_trial: scale up from frustum density to bounding-cylinder density,
        # correcting for geometric acceptance (~V_frustum/V_cyl) and a 30% wall-depletion margin.
        acceptance_est = v_frustum / v_bounding_cyl
        pf_trial = min(
            target_n * v_sphere / v_bounding_cyl / acceptance_est * 1.3,
            max_pf * 0.95,
        )

        accepted = np.empty((0, 3))
        current_seed = base_seed
        n_acc = 0
        for attempt in range(30):
            all_centers_list = openmc.model.pack_spheres(
                radius=outer_r,
                region=bounding_region,
                pf=pf_trial,
                seed=current_seed,
            )
            if len(all_centers_list) == 0:
                pf_trial = min(pf_trial * 2.0, max_pf * 0.99)
                current_seed += 1
                continue

            all_centers = np.array([[c[0], c[1], c[2]] for c in all_centers_list])
            mask = _cone_acceptance_mask(all_centers, outer_r, z_apex, tan_theta, z_cone_top)
            accepted = all_centers[mask]
            n_acc = len(accepted)

            tol = max(1, int(round(0.005 * target_n)))
            if n_acc >= target_n - tol:
                if n_acc > target_n + tol:
                    rng = np.random.default_rng(base_seed + 9999)
                    idx = rng.choice(n_acc, target_n, replace=False)
                    accepted = accepted[idx]
                break

            # Scale pf_trial proportionally to the shortfall
            pf_trial = min(pf_trial * (target_n / max(n_acc, 1)) * 1.15, max_pf * 0.99)
            current_seed += 1
        else:
            warnings.warn(
                f"exact_cone_bed: did not converge after 30 attempts; "
                f"accepted {n_acc}/{target_n} particles. Using available particles.",
                stacklevel=2,
            )

        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        np.savez(cache_path, centers=accepted)

    # Build the frustum cell — one bed cell spanning the entire cone interior.
    cone_surf = openmc.ZCone(z0=z_apex, r2=tan_theta ** 2)
    zp_cone_bot = openmc.ZPlane(z0=0.0)
    zp_cone_top_plane = openmc.ZPlane(z0=z_cone_top)

    bed_fill = _flood_fill(0.0, background_material, dry_material, z_flood)
    frustum_region = -cone_surf & +zp_cone_bot & -zp_cone_top_plane

    if use_tiled_bed:
        lattice, tile_pf_ach, _ = tiled_bed(
            lower_left=(-r_retort, -r_retort, 0.0),
            upper_right=(r_retort, r_retort, z_cone_top),
            packing_fraction=pf,
            outer_radius=outer_r,
            fill_universe=fill_univ,
            background_material=bed_fill,
            seed=base_seed,
            params=params,
            tile_size_cm=tile_size_cm,
        )
        frustum_cell = openmc.Cell(fill=lattice, region=frustum_region)
        trisos = []  # unique cells live in the tile universe, not per-particle
        n_actual = int(round(tile_pf_ach * v_frustum / v_sphere))
        V_actual_solid = n_actual * v_sphere
        pf_achieved = V_actual_solid / V_bulk
    else:
        trisos = [openmc.model.TRISO(outer_r, fill_univ, tuple(c.tolist())) for c in accepted]
        lattice = lattice_bed(
            trisos,
            lower_left=(-r_retort, -r_retort, 0.0),
            upper_right=(r_retort, r_retort, z_cone_top),
            background_material=bed_fill,
            pitch_target=pitch_target,
        )
        frustum_cell = openmc.Cell(fill=lattice, region=frustum_region)

        n_actual = len(trisos)
        V_actual_solid = n_actual * v_sphere
        pf_achieved = V_actual_solid / V_bulk

    # Overflow into retort cylinder if charge exceeds cone volume
    overflow_cells: list[openmc.Cell] = []
    extra_trisos: list = []
    bed_height_cm = z_cone_top
    remaining_solid = V_solid - V_actual_solid

    if remaining_solid > 1e-10 and V_bulk > v_frustum:
        remaining_bulk = remaining_solid / pf
        h_overflow = remaining_bulk / (math.pi * r_retort ** 2)
        z_ov_bot = z_cone_top
        z_ov_top = z_cone_top + h_overflow
        cyl_ret = openmc.ZCylinder(r=r_retort)
        zp_ov_bot = openmc.ZPlane(z0=z_ov_bot)
        zp_ov_top = openmc.ZPlane(z0=z_ov_top)
        region_ov = -cyl_ret & +zp_ov_bot & -zp_ov_top
        ov_fill = _flood_fill(z_ov_bot, background_material, dry_material, z_flood)
        if use_tiled_bed:
            lattice_ov, tile_pf_ach, _ = tiled_bed(
                lower_left=(-r_retort, -r_retort, z_ov_bot),
                upper_right=(r_retort, r_retort, z_ov_top),
                packing_fraction=pf,
                outer_radius=outer_r,
                fill_universe=fill_univ,
                background_material=ov_fill,
                seed=base_seed + 9000,
                params=params,
                tile_size_cm=tile_size_cm,
            )
            overflow_cells.append(openmc.Cell(fill=lattice_ov, region=region_ov))
            v_ov_cyl = math.pi * r_retort ** 2 * h_overflow
            n_extra_est = int(round(tile_pf_ach * v_ov_cyl / v_sphere))
            # Represent the tile-mode count in n_total via a placeholder list length.
            extra_trisos.extend([None] * n_extra_est)
        else:
            trisos_ov = pack_bed(region_ov, pf, outer_r, fill_univ, base_seed + 9000, params)
            lattice_ov = lattice_bed(
                trisos_ov,
                lower_left=(-r_retort, -r_retort, z_ov_bot),
                upper_right=(r_retort, r_retort, z_ov_top),
                background_material=ov_fill,
                pitch_target=pitch_target,
            )
            overflow_cells.append(openmc.Cell(fill=lattice_ov, region=region_ov))
            extra_trisos.extend(trisos_ov)
        bed_height_cm = z_ov_top

    # In tiled mode `trisos` is empty and `extra_trisos` holds None placeholders;
    # combine into `all_trisos` for the 'trisos' export (which no caller consumes).
    all_trisos = trisos + [t for t in extra_trisos if t is not None]
    n_total = n_actual + len(extra_trisos)
    pf_achieved = (n_total * v_sphere) / V_bulk

    return {
        'cells': [frustum_cell] + overflow_cells,
        'cone_void_cells': [],    # no staircase gaps: frustum surface is the cell boundary
        'trisos': all_trisos,
        'pf_achieved': pf_achieved,
        'V_bulk_cm3': V_bulk,
        'V_solid_cm3': V_solid,
        'bed_height_cm': bed_height_cm,
        'z_bed_bot_cm': 0.0,
        'n_particles': n_total,
        'background_mat': background_material,
    }


# ---------------------------------------------------------------------------
# Step 4 — Retort wall, cone wall, heater element, injector, outer boundary
# ---------------------------------------------------------------------------

def retort_outer_cylinder(params):
    """ZCylinder at retort outer radius (r = r_retort_outer)."""
    r = params['dimensions']['retort']['od_cm'] / 2.0
    return openmc.ZCylinder(r=r)


def outer_cone_surface(params):
    """ZCone for the outer wall of the graphite retort cone.

    Same half-angle as the inner cone; apex placed so that the surface
    radius equals r_retort_outer exactly at z = z_cone_top, giving a
    smooth transition into the retort outer cylinder.
    """
    dim = params['dimensions']
    r_out = dim['retort']['od_cm'] / 2.0
    z_ct = dim['cone']['vertical_drop_cm']
    half_rad = math.radians(dim['cone']['included_angle_deg'] / 2.0)
    tan_th = math.tan(half_rad)
    z0 = z_ct - r_out / tan_th
    return openmc.ZCone(x0=0.0, y0=0.0, z0=z0, r2=tan_th ** 2)


def heater_surfaces(params):
    """Return (inner_cyl, outer_cyl, z_plane_bot, z_plane_top) for the heater.

    Bottom of heater is flush with the retort base (z = z_cone_top).
    """
    dim = params['dimensions']
    z_ct = dim['cone']['vertical_drop_cm']
    inner = openmc.ZCylinder(r=dim['heater']['id_cm'] / 2.0)
    outer = openmc.ZCylinder(r=dim['heater']['od_cm'] / 2.0)
    z_bot = openmc.ZPlane(z0=z_ct)
    z_top = openmc.ZPlane(z0=z_ct + dim['heater']['heated_length_cm'])
    return inner, outer, z_bot, z_top


def injector_surfaces(params):
    """Return (cyls, z_top_plane, z_bot_plane) for the water-cooled injector assembly.

    The assembly spans the full cooled_length below the cone base (z = 0):
    z_top = 0 (cone base), z_bot = −cooled_length.
    The water jacket wraps the gas bore from bottom to cone base.

    cyls keys: 'gas', 'body_inner', 'cool_inner', 'cool_outer', 'body_outer'
    """
    dim = params['dimensions']
    inj = dim['injector']
    h_cl = inj['cooled_length_cm']
    cyls = {
        'gas':        openmc.ZCylinder(r=inj['gas_throat_cm'] / 2.0),
        'body_inner': openmc.ZCylinder(r=inj['body_id_cm'] / 2.0),
        'cool_inner': openmc.ZCylinder(r=inj['coolant_id_cm'] / 2.0),
        'cool_outer': openmc.ZCylinder(r=inj['coolant_od_cm'] / 2.0),
        'body_outer': openmc.ZCylinder(r=inj['body_od_cm'] / 2.0),
    }
    return cyls, openmc.ZPlane(z0=0.0), openmc.ZPlane(z0=-h_cl)


def outer_boundary_surfaces(params):
    """Return (radial_cyl, z_top, z_bot) vacuum boundary surfaces.

    Radial: heater OD (graphite felt insulation removed per step 4 revision).
    Top:    top face of graphite retort cap (z_cone_top + retort height + cap thickness).
            Cap thickness equals the cylindrical wall thickness (od − id) / 2.
    Bottom: bottom of injector body (−cooled_length).
    """
    dim = params['dimensions']
    r_out     = dim['heater']['od_cm'] / 2.0
    r_ret_in  = dim['retort']['id_cm'] / 2.0
    r_ret_out = dim['retort']['od_cm'] / 2.0
    t_wall    = r_ret_out - r_ret_in
    z_top = dim['cone']['vertical_drop_cm'] + dim['retort']['height_cm'] + t_wall
    z_bot = -dim['injector']['cooled_length_cm']
    return (
        openmc.ZCylinder(r=r_out, boundary_type='vacuum'),
        openmc.ZPlane(z0=z_top, boundary_type='vacuum'),
        openmc.ZPlane(z0=z_bot, boundary_type='vacuum'),
    )


def furnace_shell_cells(params, graphite, water):
    """Assemble all structural and void cells outboard of the bed/gas interior.

    Does NOT include the bed interior (r < r_retort_inner above z_cone_top,
    and inside the inner cone for z in [0, z_cone_top]) — those are built by
    bed_region().  This function covers everything else: retort wall, cone
    wall, vacuum gap, heater element, exterior void regions, throat, injector.

    Gas is assumed to be entirely inside the retort; all bores below the cone
    are void.  A graphite nozzle plate (r_bore to r_retort_inner, nozzle_height
    tall) closes the base of the cone interior at z = 0.  The water-cooled
    injector body sits below this plate for the remaining cooled_length −
    nozzle_height.

    Parameters
    ----------
    params       : frozen params mapping from load_params()
    graphite     : openmc.Material  graphite_structural
    water        : openmc.Material  water at coolant density (1.0 g/cm³)

    Returns
    -------
    dict with keys:
        'retort_wall'         list of openmc.Cell
        'cone_wall'           list of openmc.Cell
        'heater'              list of openmc.Cell
        'vacuum_gap'          list of openmc.Cell
        'exterior_void'       list of openmc.Cell
        'nozzle'              list of openmc.Cell  (graphite bore cap, z = −h_nozzle to 0; h_nozzle = params nozzle_height_cm)
        'injector'            list of openmc.Cell  (water-cooled body, z = −h_cl to −h_nozzle)
        'boundary_surfaces'   dict {'radial', 'top', 'bottom'}
        'all'                 flat list of all cells
    """
    dim = params['dimensions']

    # Derived scalars
    z_ct      = dim['cone']['vertical_drop_cm']
    z_rt      = z_ct + dim['retort']['height_cm']
    h_th      = dim['nozzle']['nozzle_height_cm']   # nozzle plate height
    h_cl      = dim['injector']['cooled_length_cm']
    r_ret_in  = dim['retort']['id_cm'] / 2.0
    r_ret_out = dim['retort']['od_cm'] / 2.0
    t_wall    = r_ret_out - r_ret_in    # wall thickness = top cap thickness
    z_cap_top = z_rt + t_wall           # outer (vacuum) face of the top cap
    r_hi   = dim['heater']['id_cm'] / 2.0
    r_ho   = dim['heater']['od_cm'] / 2.0
    r_thr  = dim['nozzle']['throat_diameter_cm'] / 2.0
    h_heat = dim['heater']['heated_length_cm']

    half_rad   = math.radians(dim['cone']['included_angle_deg'] / 2.0)
    tan_th     = math.tan(half_rad)
    z_apex_in  = -r_thr / tan_th
    z_apex_out = z_ct - r_ret_out / tan_th

    inj        = dim['injector']
    r_inj_gas  = inj['gas_throat_cm'] / 2.0
    r_inj_ci   = inj['coolant_id_cm'] / 2.0
    r_inj_co   = inj['coolant_od_cm'] / 2.0
    r_inj_bo   = inj['body_od_cm'] / 2.0

    # ----- Surfaces -----------------------------------------------------------
    s_cone_in  = openmc.ZCone(z0=z_apex_in,  r2=tan_th ** 2)
    s_cone_out = openmc.ZCone(z0=z_apex_out, r2=tan_th ** 2)

    s_cyl_ret_in  = openmc.ZCylinder(r=r_ret_in)
    s_cyl_ret_out = openmc.ZCylinder(r=r_ret_out)
    s_cyl_hi      = openmc.ZCylinder(r=r_hi)
    s_cyl_ho      = openmc.ZCylinder(r=r_ho, boundary_type='vacuum')

    s_inj_gas = openmc.ZCylinder(r=r_inj_gas)
    s_inj_ci  = openmc.ZCylinder(r=r_inj_ci)
    s_inj_co  = openmc.ZCylinder(r=r_inj_co)
    s_inj_bo  = openmc.ZCylinder(r=r_inj_bo)

    zp_ct      = openmc.ZPlane(z0=z_ct)
    zp_rt      = openmc.ZPlane(z0=z_rt)                             # retort interior top
    zp_cap_top = openmc.ZPlane(z0=z_cap_top, boundary_type='vacuum')  # outer cap face
    zp_ht      = openmc.ZPlane(z0=z_ct + h_heat)
    zp_cb = openmc.ZPlane(z0=0.0)           # cone base = nozzle plate top
    zp_tb = openmc.ZPlane(z0=-h_th)         # nozzle plate bottom = injector body top
    zp_ib = openmc.ZPlane(z0=-h_cl, boundary_type='vacuum')

    # ----- Cells --------------------------------------------------------------

    retort_wall = [
        openmc.Cell(
            name='retort_cyl_wall',
            fill=graphite,
            region=+s_cyl_ret_in & -s_cyl_ret_out & +zp_ct & -zp_rt,
        )
    ]

    cone_wall = [
        openmc.Cell(
            name='cone_graphite_wall',
            fill=graphite,
            region=+s_cone_in & -s_cone_out & +zp_cb & -zp_ct,
        )
    ]

    vacuum_gap = [
        openmc.Cell(
            name='heater_vacuum_gap',
            fill=None,
            region=+s_cyl_ret_out & -s_cyl_hi & +zp_ct & -zp_ht,
        )
    ]

    heater = [
        openmc.Cell(
            name='heater_element',
            fill=graphite,
            region=+s_cyl_hi & -s_cyl_ho & +zp_ct & -zp_ht,
        )
    ]

    retort_top_cap = [
        openmc.Cell(
            name='retort_top_cap',
            fill=graphite,
            region=-s_cyl_ret_out & +zp_rt & -zp_cap_top,
        )
    ]

    exterior_void = [
        openmc.Cell(
            name='void_above_heater',
            fill=None,
            region=+s_cyl_ret_out & -s_cyl_ho & +zp_ht & -zp_cap_top,
        ),
        openmc.Cell(
            name='void_outside_cone',
            fill=None,
            region=+s_cone_out & -s_cyl_ho & +zp_cb & -zp_ct,
        ),
    ]

    # Nozzle bore cap: graphite plug that closes the bore at the cone base.
    # Tracked via params['dimensions']['nozzle']['nozzle_height_cm'] (h_th = 1.5 cm).
    # The bore (r < r_thr) is void below z = -h_th; graphite above it to z = 0.
    nozzle = [
        openmc.Cell(
            name='nozzle_bore_cap',
            fill=graphite,
            region=-s_inj_gas & +zp_tb & -zp_cb,
        ),
    ]

    # Water-cooled injector assembly.
    # The bore (r < r_thr) is void only in z ∈ [-h_cl, -h_th]; it is closed by
    # nozzle_bore_cap above.  All other annular cells run the full cooled length
    # z ∈ [-h_cl, 0] so the water jacket extends into the nozzle zone.
    injector = [
        openmc.Cell(
            name='injector_gas_bore',
            fill=None,
            region=-s_inj_gas & +zp_ib & -zp_tb,
        ),
        openmc.Cell(
            name='injector_inner_graphite',
            fill=graphite,
            region=+s_inj_gas & -s_inj_ci & +zp_ib & -zp_cb,
        ),
        openmc.Cell(
            name='injector_coolant',
            fill=water,
            region=+s_inj_ci & -s_inj_co & +zp_ib & -zp_cb,
        ),
        openmc.Cell(
            name='injector_outer_graphite',
            fill=graphite,
            region=+s_inj_co & -s_inj_bo & +zp_ib & -zp_cb,
        ),
        openmc.Cell(
            name='void_injector_exterior',
            fill=None,
            region=+s_inj_bo & -s_cyl_ho & +zp_ib & -zp_cb,
        ),
    ]

    all_cells = (
        retort_wall + retort_top_cap + cone_wall + vacuum_gap + heater
        + exterior_void + nozzle + injector
    )

    return {
        'retort_wall':    retort_wall,
        'retort_top_cap': retort_top_cap,
        'cone_wall':      cone_wall,
        'vacuum_gap':     vacuum_gap,
        'heater':         heater,
        'exterior_void':  exterior_void,
        'nozzle':         nozzle,
        'injector':       injector,
        'boundary_surfaces': {
            'radial': s_cyl_ho,
            'top':    zp_cap_top,
            'bottom': zp_ib,
        },
        'all': all_cells,
    }
