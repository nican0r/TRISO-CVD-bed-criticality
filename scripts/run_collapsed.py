"""Step 5 collapsed-bed case: collapsed bed, bare UCO kernel, process gas, 293.6 K.

Companion to run_nominal.py. Uses the same charge, materials, and settings, but
with state='collapsed' so the bed occupies only the cone at pf_static (denser,
lower bed height). Together the two runs bound the k-eff between the fluidized
(operating) and collapsed (settled) bed geometries.

Run:
    caffeinate python scripts/run_collapsed.py
"""
from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import openmc

from furnace.params import check_env, load_params
from furnace.model import build_model, export_and_run

_OUT_DIR = Path('results/step5_collapsed')


def main() -> None:
    check_env()
    params = load_params()

    model, stats = build_model(params, state='collapsed')

    mdl = params['model']
    n_inactive = int(mdl['inactive'])
    n_active   = int(mdl['batches']) - n_inactive

    print('\nCollapsed case — settled bed, bare UCO kernel, 293.6 K, process gas')
    print(f'  U-235 mass   : {stats.u235_mass_g:.4f} g')
    print(f'  Charge mass  : {params["dimensions"]["bed"]["charge_mass_g"]:.1f} g')
    print(f'  Particles    : {stats.n_particles}')
    print(f'  PF achieved  : {stats.pf_achieved:.4f}')
    print(f'  Bed height   : {stats.bed_height_cm:.3f} cm')
    print(f'  Batches      : {n_inactive} inactive + {n_active} active')
    print(f'  Particles/gen: {mdl["particles"]}')
    print(f'  Seed         : {mdl["seed"]}')
    print(f'\nRunning OpenMC eigenvalue calculation → {_OUT_DIR}/')

    sp_path = export_and_run(model, _OUT_DIR)

    with openmc.StatePoint(str(sp_path)) as sp:
        keff = sp.keff

    print(f'\n{"=" * 52}')
    print(f'U-235 mass : {stats.u235_mass_g:.4f} g')
    print(f'k-eff      : {keff.nominal_value:.5f} ± {keff.std_dev:.5f}')
    print(f'k-eff + 2σ : {keff.nominal_value + 2 * keff.std_dev:.5f}')
    print(f'k-eff + 3σ : {keff.nominal_value + 3 * keff.std_dev:.5f}')
    print(f'{"=" * 52}')
    print(f'\nStagepoint written to: {sp_path}')


if __name__ == '__main__':
    main()
