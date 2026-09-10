"""Step 5 convergence study: batch-count comparison (250 vs 500 active) + 3-seed test.

Run AFTER run_nominal.py (nominal provides the seed=42, 250-active reference point).

Cases produced:
    results/step5_convergence/batch_500/  — seed=42, 500 active batches
    results/step5_convergence/seed_43/   — seed=43, 250 active batches, different packing
    results/step5_convergence/seed_44/   — seed=44, 250 active batches, different packing

Batch comparison: k(250 active) vs k(500 active) — must agree within combined σ.
Seed comparison : k across seeds 42, 43, 44 — must agree within combined σ.

Run:
    caffeinate python scripts/run_step5_convergence.py
"""
from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import openmc

from furnace.params import check_env, load_params
from furnace.model import build_model, export_and_run

_BASE_DIR  = Path('results/step5_convergence')
_NOM_DIR   = Path('results/step5_nominal')

_CASES = [
    {'tag': 'batch_500', 'seed': 42, 'n_active': 500},
    {'tag': 'seed_43',   'seed': 43, 'n_active': 250},
    {'tag': 'seed_44',   'seed': 44, 'n_active': 250},
]


def _run_case(params, case: dict) -> dict:
    seed     = case['seed']
    n_active = case['n_active']
    out_dir  = _BASE_DIR / case['tag']

    n_inactive = int(params['model']['inactive'])
    model, stats = build_model(params, seed=seed, n_active=n_active)

    print(f"\n[{case['tag']}]  seed={seed}  inactive={n_inactive}  active={n_active}")
    print(f"  U-235 mass : {stats.u235_mass_g:.4f} g  |  particles: {stats.n_particles}")

    sp_path = export_and_run(model, out_dir)

    with openmc.StatePoint(str(sp_path)) as sp:
        keff = sp.keff

    print(f"  k-eff      : {keff.nominal_value:.5f} ± {keff.std_dev:.5f}  "
          f"(U-235 mass {stats.u235_mass_g:.4f} g)")

    return {
        'tag':      case['tag'],
        'seed':     seed,
        'n_active': n_active,
        'k_nom':    keff.nominal_value,
        'k_sig':    keff.std_dev,
        'sp_path':  str(sp_path),
        'u235_g':   stats.u235_mass_g,
    }


def _report_batch_comparison(params, results: list[dict]) -> None:
    """Compare k-eff at 250 vs 500 active batches (seed=42 only)."""
    # 250-active result from the nominal run statepoint
    n_inactive = int(params['model']['inactive'])
    nom_sp = _NOM_DIR / f'statepoint.{n_inactive + 250}.h5'
    if not nom_sp.exists():
        print('\n[Batch comparison] Nominal statepoint not found — run run_nominal.py first.')
        return

    with openmc.StatePoint(str(nom_sp)) as sp:
        k250 = sp.keff

    # 500-active result
    b500 = next((r for r in results if r['tag'] == 'batch_500'), None)
    if b500 is None:
        return

    k500_nom = b500['k_nom']
    k500_sig = b500['k_sig']

    diff   = abs(k250.nominal_value - k500_nom)
    combined_sig = (k250.std_dev**2 + k500_sig**2) ** 0.5

    print('\n── Batch comparison (seed=42) ───────────────────────────')
    print(f'  k(250 active) : {k250.nominal_value:.5f} ± {k250.std_dev:.5f}')
    print(f'  k(500 active) : {k500_nom:.5f} ± {k500_sig:.5f}')
    print(f'  |Δk|          : {diff:.5f}   combined σ : {combined_sig:.5f}')
    status = 'PASS — no drift detected' if diff < combined_sig else 'FAIL — k-eff drifting, investigate'
    print(f'  Criterion (|Δk| < combined σ): {status}')


def _report_seed_comparison(params, results: list[dict]) -> None:
    """Compare k-eff across seeds 42, 43, 44."""
    n_inactive = int(params['model']['inactive'])
    nom_sp = _NOM_DIR / f'statepoint.{n_inactive + 250}.h5'
    if not nom_sp.exists():
        print('\n[Seed comparison] Nominal statepoint not found.')
        return

    with openmc.StatePoint(str(nom_sp)) as sp:
        k42 = sp.keff

    seed_results = [
        ('seed=42 (nominal)', k42.nominal_value, k42.std_dev),
    ]
    for r in results:
        if r['tag'] in ('seed_43', 'seed_44'):
            seed_results.append((f"seed={r['seed']}", r['k_nom'], r['k_sig']))

    print('\n── Seed comparison (three independent geometries) ────────')
    for label, k, sig in seed_results:
        print(f'  {label:<22}: {k:.5f} ± {sig:.5f}')

    # Check all pairs agree within combined σ
    pass_all = True
    for i in range(len(seed_results)):
        for j in range(i + 1, len(seed_results)):
            la, ka, sa = seed_results[i]
            lb, kb, sb = seed_results[j]
            diff = abs(ka - kb)
            comb = (sa**2 + sb**2)**0.5
            ok   = diff < comb
            pass_all = pass_all and ok
            print(f'  |k({i+1}) - k({j+1})| = {diff:.5f},  combined σ = {comb:.5f}  '
                  f'→ {"OK" if ok else "FAIL"}')
    print(f'  Overall: {"PASS" if pass_all else "FAIL — seeds disagree, investigate geometry or source convergence"}')


def main() -> None:
    check_env()
    params = load_params()
    _BASE_DIR.mkdir(parents=True, exist_ok=True)

    results = [_run_case(params, case) for case in _CASES]

    print('\n' + '=' * 60)
    print('CONVERGENCE SUMMARY')
    print('=' * 60)

    _report_batch_comparison(params, results)
    _report_seed_comparison(params, results)

    print('\nNext: generate convergence plots:')
    print('    python scripts/plot_step5_convergence.py')


if __name__ == '__main__':
    main()
