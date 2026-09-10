"""Step 4 geometry check and XZ-slice plots for the full furnace shell.

Bed interior is filled with process_gas (no TRISO packing) for fast validation.
OpenMC runs in plot mode: overlapping cells produce an error; undefined regions
appear as black pixels in the output PNGs.

Run:
    caffeinate python scripts/check_geometry.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import openmc

# openmc binary lives next to the Python interpreter in the same conda env
_OPENMC_EXEC = str(Path(sys.executable).parent / 'openmc')

from furnace.params import check_env, load_params
from furnace.materials import graphite_structural, process_gas, water
import furnace.geometry as geom


def _simplified_bed_cells(params, pg_mat, shell_bounds):
    """Two cells that fill the bed interior without TRISO packing.

    Covers:
      - Cone interior: inside inner ZCone, z in [0, z_cone_top]
      - Retort cylinder interior: r < r_retort_inner, z in [z_cone_top, z_retort_top]
    """
    dim = params['dimensions']
    z_ct = dim['cone']['vertical_drop_cm']
    r_ret_in = dim['retort']['id_cm'] / 2.0

    half_rad = math.radians(dim['cone']['included_angle_deg'] / 2.0)
    tan_th = math.tan(half_rad)
    r_thr = dim['nozzle']['throat_diameter_cm'] / 2.0
    z_apex_in = -r_thr / tan_th

    s_cone_in = openmc.ZCone(z0=z_apex_in, r2=tan_th ** 2)
    zp_cb = openmc.ZPlane(z0=0.0)
    zp_ct = openmc.ZPlane(z0=z_ct)
    s_cyl_ret = openmc.ZCylinder(r=r_ret_in)

    cone_interior = openmc.Cell(
        name='cone_interior_simplified',
        fill=pg_mat,
        region=-s_cone_in & +zp_cb & -zp_ct,
    )
    retort_interior = openmc.Cell(
        name='retort_interior_simplified',
        fill=pg_mat,
        region=-s_cyl_ret & +zp_ct & -shell_bounds['top'],
    )
    return [cone_interior, retort_interior]


def main():
    check_env()
    params = load_params()

    g_mat = graphite_structural(params)
    w_mat = water(1.0)
    pg_mat = process_gas(params)

    shell = geom.furnace_shell_cells(params, graphite=g_mat, water=w_mat)
    bed_cells = _simplified_bed_cells(params, pg_mat, shell['boundary_surfaces'])

    universe = openmc.Universe(cells=bed_cells + shell['all'])
    geometry = openmc.Geometry(universe)

    materials = openmc.Materials([g_mat, w_mat, pg_mat])

    settings = openmc.Settings()
    settings.run_mode = 'plot'

    # Plot geometry: centre of model axially
    dim = params['dimensions']
    z_ct = dim['cone']['vertical_drop_cm']
    z_rt = z_ct + dim['retort']['height_cm']
    h_cl = dim['injector']['cooled_length_cm']
    z_mid = (z_rt - h_cl) / 2.0

    # Full-model XZ: covers heater OD (9.4 cm diameter) and full axial span
    p_full = openmc.Plot(name='full_model')
    p_full.filename = 'full_model_xz'
    p_full.origin = (0.0, 0.0, z_mid)
    p_full.width = (11.0, z_rt + h_cl + 2.0)   # small margin top and bottom
    p_full.pixels = (550, round(550 * (z_rt + h_cl + 2.0) / 11.0))
    p_full.basis = 'xz'
    p_full.color_by = 'material'

    # Zoom: cone + throat + injector (z ≈ -7 to +7 cm)
    p_zoom = openmc.Plot(name='cone_injector_zoom')
    p_zoom.filename = 'cone_injector_xz'
    p_zoom.origin = (0.0, 0.0, 0.0)
    p_zoom.width = (11.0, 14.0)
    p_zoom.pixels = (550, 700)
    p_zoom.basis = 'xz'
    p_zoom.color_by = 'material'

    plots = openmc.Plots([p_full, p_zoom])

    out_dir = Path('results/step4_geometry_check')
    out_dir.mkdir(parents=True, exist_ok=True)

    geometry.export_to_xml(path=str(out_dir / 'geometry.xml'))
    materials.export_to_xml(path=str(out_dir / 'materials.xml'))
    settings.export_to_xml(path=str(out_dir / 'settings.xml'))
    plots.export_to_xml(path=str(out_dir / 'plots.xml'))

    print('Running OpenMC in plot mode (geometry check)...')
    openmc.run(cwd=str(out_dir), openmc_exec=_OPENMC_EXEC)
    print(f'\nGeometry check complete.')
    print(f'Plots written to: {out_dir}/')
    print('  full_model_xz.png     — full furnace, XZ slice')
    print('  cone_injector_xz.png  — cone + throat + injector zoom')
    print('Any black pixels indicate undefined geometry regions.')
    print('An OpenMC error above would indicate overlapping cells.')


if __name__ == '__main__':
    main()
