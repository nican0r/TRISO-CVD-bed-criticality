"""Smoke test: exercise the full run_nominal code path with a tiny geometry.

Uses charge_mass_g=0.5 (~500 TRISO particles, ~1 MB geometry XML) and
the minimum possible batch/particle count so the full pipeline completes
in under 2 minutes.  Catches cross-section, geometry, and tally errors
before committing to a full production run.

Run from repo root:
    python3 scripts/smoke_test.py
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

_OUT_DIR = Path('results/smoke_test')


def main() -> None:
    check_env()
    params = load_params()

    print('Smoke test — 5.0 g charge, 1 inactive + 1 active batch, 10 000 particles')
    print('Catches: cross-section errors, geometry errors, tally errors')
    print('Does NOT validate k-eff or convergence.\n')

    model, stats = build_model(
        params,
        charge_mass_g=5.0,   # minimum that gives pack_spheres enough headroom (~14 particle diameters)
        n_inactive=1,
        n_active=1,
        n_particles=10_000,  # Watt source at 1–2 MeV has low fission p per particle in tiny kernel
                             # (~0.2%); 10 000 particles gives ~20 expected fissions/gen to bank
    )

    print(f'  Particles in bed : {stats.n_particles}')
    print(f'  U-235 mass       : {stats.u235_mass_g:.5f} g')
    print(f'  Bed height       : {stats.bed_height_cm:.3f} cm')
    print(f'\nExporting XML and running OpenMC → {_OUT_DIR}/')

    sp_path = export_and_run(model, _OUT_DIR)

    with openmc.StatePoint(str(sp_path)) as sp:
        keff = sp.keff

    print(f'\nSmoke test PASSED')
    print(f'  k-eff (meaningless at 10k particles, 5g): {keff.nominal_value:.4f} ± {keff.std_dev:.4f}')
    print(f'\nFull run is safe to start:')
    print('  caffeinate python3 scripts/run_nominal.py')


if __name__ == '__main__':
    main()
