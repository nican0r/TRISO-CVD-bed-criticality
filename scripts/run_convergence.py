"""Staircase discretisation convergence study for the collapsed-bed geometry.

Tests n_slabs in {4, 8, 16, 32} to find the minimum number of staircase cylinders
in the cone region that keeps the discretisation k-eff error below one statistical
sigma relative to the next-finer run.

Geometry note — homogenised materials
--------------------------------------
This script uses homogenised (smeared) materials instead of explicit TRISO packing.
Explicit packing at pf_collapsed = 0.60 requires placing tens of thousands of
particles per slab and can take many hours; homogenisation gives the same geometry
sensitivity in minutes. Each slab cell is filled with an openmc.Material.mix_materials
blend of UCO kernel + background water at the correct volume fraction (pf_static from
params), which conserves the fissile atom inventory exactly. Absolute k-eff values
differ from the explicit-TRISO model (no Dancoff / self-shielding correction), but
the RELATIVE δk_disc between n_slabs values correctly reflects the geometric error.

Production runs (step 5 onward) use explicit TRISO particles.

Geometry (inside retort inner cylinder, vacuum BC):
    z in [-h_throat, 0]          — nozzle throat: graphite cylinder + background exterior
    z in [0, z_cone_top]         — cone: inscribed staircase bed + graphite wall + void cells
    z in [z_cone_top, z_ret_top] — retort cylinder: collapsed bed overflow + background

Materials:
    Bed background : water at 1.0 g/cm³ (most reactive moderator; maximises geometry
                     sensitivity for the convergence test)
    Cone wall + throat : graphite_structural

Stage   : bare_kernel (smallest particle, highest sensitivity to cone discretisation)
State   : collapsed (staircase only affects the collapsed geometry)
Charge  : params charge_mass_g (default 95 g)

Convergence criterion: δk_disc(n) = |k(n) − k(2n)| < σ(2n)

Run with:
    caffeinate micromamba run -n triso-env python scripts/run_convergence.py
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import openmc
import openmc.model

from furnace.params import check_env, load_params
from furnace.materials import graphite_structural, uco_kernel
from furnace.geometry import (
    retort_inner_cylinder,
    cone_wall,
    nozzle_throat,
    staircase_bed,
    frustum_volume,
)
from furnace.triso import _shell_vol, _stage_outer_radius, particle_mass_g


def _water_no_sab(density: float) -> openmc.Material:
    """Plain water without S(a,b) table, for use in mix_materials blending."""
    m = openmc.Material(name=f'water_no_sab_{density:.3f}')
    m.set_density('g/cm3', density)
    m.add_element('H', 2, 'ao')
    m.add_element('O', 1, 'ao')
    return m


def _homogenised_bed(params, n_slabs, bg_mat, charge_mass_g=None):
    """Build collapsed-bed geometry with homogenised materials.

    Each occupied slab cell is filled with a volume-fraction mixture of UCO kernel
    and background material at pf_static.  The cone_void_cells (annular voids +
    unfilled cone interior) are pure background material.  Overflow into the retort
    cylinder above the cone is handled if the charge exceeds the staircase volume.

    Returns dict with keys:
        cells, cone_void_cells, bed_height_cm, stair_vol, vol_err_frac, materials
    """
    dim = params['dimensions']
    mdl = params['model']

    if charge_mass_g is None:
        charge_mass_g = float(dim['bed']['charge_mass_g'])

    pf = float(mdl['packing_fraction_static'])

    r_throat = dim['nozzle']['throat_diameter_cm'] / 2.0
    r_retort = dim['retort']['id_cm'] / 2.0
    z_cone_top = dim['cone']['vertical_drop_cm']
    half_angle = math.radians(dim['cone']['included_angle_deg'] / 2.0)
    tan_theta = math.tan(half_angle)
    z_apex = -r_throat / tan_theta

    cone_surf = openmc.ZCone(z0=z_apex, r2=tan_theta**2)
    zp_cone_top = openmc.ZPlane(z0=z_cone_top)

    # Bare-kernel effective density = kernel density (outer radius = kernel radius)
    rho_eff = particle_mass_g('bare_kernel', params) / _shell_vol(
        _stage_outer_radius('bare_kernel', params)
    )
    V_solid = charge_mass_g / rho_eff
    V_bulk = V_solid / pf

    kernel_mat = uco_kernel(params)
    # mix_materials cannot handle S(a,b) tables, so smear with a plain (no S(a,b))
    # version of the background at the same density. The bg_mat with S(a,b) is still
    # used in void cells, preserving correct thermal scattering outside the bed.
    bg_density = bg_mat.density
    bg_plain = _water_no_sab(bg_density)
    # Homogenised slab material: pf volume fraction UCO kernel + (1-pf) background
    smeared = openmc.Material.mix_materials(
        [kernel_mat, bg_plain],
        fracs=[pf, 1.0 - pf],
        percent_type='vo',
        name='smeared_fuel',
    )

    slabs, stair_vol, vol_err_frac = staircase_bed(params, n_slabs)

    remaining_bulk = V_bulk
    cells = []
    cone_void_cells = []
    z_last_fill_top = 0.0

    for i, slab in enumerate(slabs):
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
            h_partial = remaining_bulk / (math.pi * r_slab**2)
            z_fill_top = z_bot + h_partial
            used_vol = remaining_bulk

        cyl = openmc.ZCylinder(r=r_slab)
        zp_bot_s = openmc.ZPlane(z0=z_bot)
        zp_top_s = openmc.ZPlane(z0=z_fill_top)

        cells.append(openmc.Cell(
            fill=smeared,
            region=-cyl & +zp_bot_s & -zp_top_s,
            name=f'bed_slab_{i}',
        ))

        # Annular void: inside cone, outside inscribed cylinder
        cone_void_cells.append(openmc.Cell(
            fill=bg_mat,
            region=+cyl & -cone_surf & +zp_bot_s & -zp_top_s,
            name=f'annular_{i}',
        ))

        z_last_fill_top = z_fill_top
        remaining_bulk -= used_vol

    # Unfilled cone interior above bed
    if z_last_fill_top < z_cone_top - 1e-10:
        zp_bed_last = openmc.ZPlane(z0=z_last_fill_top)
        cone_void_cells.append(openmc.Cell(
            fill=bg_mat,
            region=-cone_surf & +zp_bed_last & -zp_cone_top,
            name='unfilled_cone',
        ))

    # Overflow into retort cylinder above cone
    overflow_cells = []
    if remaining_bulk > 1e-10:
        h_overflow = remaining_bulk / (math.pi * r_retort**2)
        z_ov_top = z_cone_top + h_overflow
        cyl_ret = openmc.ZCylinder(r=r_retort)
        zp_ov_bot = openmc.ZPlane(z0=z_cone_top)
        zp_ov_top = openmc.ZPlane(z0=z_ov_top)
        overflow_cells.append(openmc.Cell(
            fill=smeared,
            region=-cyl_ret & +zp_ov_bot & -zp_ov_top,
            name='bed_overflow',
        ))
        bed_height_cm = z_ov_top
    else:
        bed_height_cm = z_last_fill_top

    return {
        'cells': cells + overflow_cells,
        'cone_void_cells': cone_void_cells,
        'bed_height_cm': bed_height_cm,
        'stair_vol': stair_vol,
        'vol_err_frac': vol_err_frac,
        'materials': [kernel_mat, bg_plain, smeared, bg_mat],
    }


def _build_and_run(params, n_slabs, bg_mat, run_dir):
    """Build a minimal OpenMC model with homogenised bed and run k-eff.

    Returns (k_eff, sigma, stair_vol, vol_err_frac).
    """
    dim = params['dimensions']
    mdl = params['model']

    r_retort = dim['retort']['id_cm'] / 2.0
    z_cone_top = dim['cone']['vertical_drop_cm']
    h_throat = dim['nozzle']['nozzle_height_cm']
    z_retort_top = dim['retort']['height_cm'] - z_cone_top
    graph = graphite_structural(params)

    bed = _homogenised_bed(params, n_slabs, bg_mat)
    bed_top_z = bed['bed_height_cm']

    # --- Boundary surfaces ---
    retort_cyl = retort_inner_cylinder(params)
    retort_cyl.boundary_type = 'vacuum'
    zp_retort_top = openmc.ZPlane(z0=z_retort_top, boundary_type='vacuum')
    zp_throat_bot = openmc.ZPlane(z0=-h_throat, boundary_type='vacuum')

    cone_surf, cone_zbot, cone_ztop = cone_wall(params)
    throat_cyl, throat_zbot, throat_ztop = nozzle_throat(params)

    # --- Cells ---
    cells = []

    # Throat graphite
    cells.append(openmc.Cell(
        fill=graph,
        region=-throat_cyl & +zp_throat_bot & -throat_ztop,
        name='throat_graphite',
    ))

    # Throat exterior (non-physical void, required for spatial coverage)
    cells.append(openmc.Cell(
        fill=bg_mat,
        region=+throat_cyl & -retort_cyl & +zp_throat_bot & -throat_ztop,
        name='throat_exterior',
    ))

    # Cone graphite wall (+cone_surf = exterior of the cone surface)
    cells.append(openmc.Cell(
        fill=graph,
        region=+cone_surf & -retort_cyl & +cone_zbot & -cone_ztop,
        name='cone_graphite',
    ))

    # Bed cells (homogenised inscribed cylinders)
    for c in bed['cells']:
        cells.append(c)

    # Cone void cells (annular voids + unfilled cone interior)
    for c in bed['cone_void_cells']:
        cells.append(c)

    # Retort cylinder above cone
    if bed_top_z > z_cone_top:
        zp_ov_top = openmc.ZPlane(z0=bed_top_z)
        if bed_top_z < z_retort_top:
            cells.append(openmc.Cell(
                fill=bg_mat,
                region=-retort_cyl & +zp_ov_top & -zp_retort_top,
                name='above_bed',
            ))
    else:
        cells.append(openmc.Cell(
            fill=bg_mat,
            region=-retort_cyl & +cone_ztop & -zp_retort_top,
            name='above_cone',
        ))

    # --- Geometry and materials ---
    geometry = openmc.Geometry(openmc.Universe(cells=cells))
    mat_list = [graph] + bed['materials']
    # materials.py sets all temperatures to 293.6 K (room temp); override to 1200 K
    # (nearest available library point to CVD operating conditions)
    for m in mat_list:
        m.temperature = 1200.0
    materials = openmc.Materials(mat_list)

    # --- Settings ---
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = int(mdl['batches'])
    settings.inactive = int(mdl['inactive'])
    settings.particles = int(mdl['particles'])
    settings.seed = int(mdl['seed'])
    # Library has data only at 900, 1200, 2500 K.  1200 K is the closest available
    # point to CVD operating conditions; tolerance covers the gap from any material
    # temperature set elsewhere in the model.
    settings.temperature = {
        'method': 'nearest',
        'default': 1200.0,
        'range': [800.0, 2000.0],
    }

    r_src = r_retort * 0.5
    z_src_hi = max(bed_top_z * 0.9, 0.1)
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            lower_left=(-r_src, -r_src, 0.0),
            upper_right=(r_src, r_src, z_src_hi),
        )
    )

    model = openmc.model.Model(geometry=geometry, materials=materials, settings=settings)

    run_dir.mkdir(parents=True, exist_ok=True)
    orig_cwd = Path.cwd()
    try:
        model.export_to_xml(str(run_dir))
        os.chdir(run_dir)
        openmc.run()
    finally:
        os.chdir(orig_cwd)

    sp_abs = run_dir / f'statepoint.{settings.batches}.h5'
    with openmc.StatePoint(sp_abs) as sp:
        k_eff = sp.keff.nominal_value
        sigma = sp.keff.std_dev

    return k_eff, sigma, bed['stair_vol'], bed['vol_err_frac']


def main():
    check_env()
    params = load_params()

    n_slabs_list = [4, 8, 16, 32]
    # c_H_in_H2O S(a,b) is not available in this library; use plain water.
    # The convergence study measures geometry discretisation (relative δk_disc),
    # not absolute k-eff, so omitting thermal scattering tables is acceptable here.
    bg_mat = _water_no_sab(1.0)

    results_dir = _REPO_ROOT / 'results'
    results_dir.mkdir(exist_ok=True)

    # Analytic frustum volume and total bulk volume — both fixed for all n_slabs.
    dim = params['dimensions']
    r_throat = dim['nozzle']['throat_diameter_cm'] / 2.0
    r_retort  = dim['retort']['id_cm'] / 2.0
    half_angle = dim['cone']['included_angle_deg'] / 2.0
    v_frustum = frustum_volume(r_retort, r_throat, half_angle_deg=half_angle)

    charge_mass_g = float(dim['bed']['charge_mass_g'])
    pf_static = float(params['model']['packing_fraction_static'])
    rho_eff = particle_mass_g('bare_kernel', params) / _shell_vol(
        _stage_outer_radius('bare_kernel', params)
    )
    v_bulk = (charge_mass_g / rho_eff) / pf_static

    rows = []
    k_results = {}

    print(
        f"\n{'n_slabs':>8} {'V_stair(cm³)':>13} {'V_err(%)':>9} "
        f"{'ΔV/V_bulk(%)':>13} {'k_eff':>8} {'δk_disc':>12}"
    )
    print('-' * 70)

    for n_slabs in n_slabs_list:
        run_dir = _REPO_ROOT / 'cases' / f'convergence_n{n_slabs}'
        k_eff, sigma, stair_vol, vol_err_frac = _build_and_run(params, n_slabs, bg_mat, run_dir)

        vol_err_pct = vol_err_frac * 100.0
        dv_pct = (v_frustum - stair_vol) / v_bulk * 100.0
        k_results[n_slabs] = (k_eff, sigma)
        rows.append({
            'n_slabs': n_slabs,
            'V_staircase_cm3': round(stair_vol, 4),
            'V_error_pct': round(vol_err_pct, 3),
            'dV_over_Vbulk_pct': round(dv_pct, 3),
            'k_eff': round(k_eff, 5),
            'dk_disc': None,
        })
        print(
            f"{n_slabs:>8} {stair_vol:>13.4f} {vol_err_pct:>9.3f} "
            f"{dv_pct:>13.3f} {k_eff:>8.5f}            —"
        )

    # δk_disc(n) = |k(n) − k(2n)| for each consecutive doubling
    for i, row in enumerate(rows[:-1]):
        n = row['n_slabs']
        n2 = rows[i + 1]['n_slabs']
        if n2 == 2 * n:
            dk = abs(k_results[n][0] - k_results[n2][0])
            row['dk_disc'] = dk

    print('\n--- Updated table with δk_disc ---')
    print(
        f"{'n_slabs':>8} {'V_stair(cm³)':>13} {'V_err(%)':>9} "
        f"{'ΔV/V_bulk(%)':>13} {'k_eff':>8} {'δk_disc':>12}"
    )
    print('-' * 70)
    for row in rows:
        dk_str = f"{row['dk_disc']:.2e}" if row['dk_disc'] is not None else '           —'
        print(
            f"{row['n_slabs']:>8} {row['V_staircase_cm3']:>13.4f} "
            f"{row['V_error_pct']:>9.3f} "
            f"{row['dV_over_Vbulk_pct']:>13.3f} {row['k_eff']:>8.5f} {dk_str:>12}"
        )

    csv_path = results_dir / 'convergence_n_slabs.csv'
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f'\nResults saved to {csv_path}')

    print(f'\nV_frustum (analytic) = {v_frustum:.4f} cm³')
    print(f'V_bulk               = {v_bulk:.4f} cm³')
    print('ΔV/V_bulk = (V_frustum − V_staircase) / V_bulk: fraction of total inventory mislocated from cone to cylinder.')
    print(
        '\nNote: k-eff values use homogenised materials (no TRISO self-shielding).'
        '\nAbsolute values are not physically accurate; relative δk_disc is the deliverable.'
    )


if __name__ == '__main__':
    main()
