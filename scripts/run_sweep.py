#!/usr/bin/env python
"""Parametric sweep driver. Step 6: seed-check sweep. Steps 7+ add mass/packing sweeps.

Usage:
    # Smoke-test the sweep infrastructure (3 quick cases varying only seed):
    caffeinate python scripts/run_sweep.py --seed-check --quick

    # Full seed-check at nominal settings:
    caffeinate python scripts/run_sweep.py --seed-check

    # Parallel cases (N concurrent processes):
    caffeinate python scripts/run_sweep.py --seed-check --quick --jobs 3

    # Re-run cases even if already in CSV:
    caffeinate python scripts/run_sweep.py --seed-check --quick --force
"""
from __future__ import annotations

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from furnace.params import check_env, load_params
from furnace.sweeps import run_case, vessel_capacity_g

_RESULTS_DIR = _REPO_ROOT / 'results'
_SEED_CHECK_CSV = _RESULTS_DIR / 'seed_check.csv'
_SEED_CHECK_DIR = _RESULTS_DIR / 'seed_check'

_MASS_SWEEP_CSV = _RESULTS_DIR / 'mass_sweep.csv'
_MASS_SWEEP_DIR = _RESULTS_DIR / 'mass_sweep'

# Nominal charge is 95 g (params.yaml).  Multipliers span 1× to 40× nominal.
# Grid spans 2× – 38× the nominal 95 g charge (190 g – 3610 g).
# 1× and 1.5× are excluded: the nominal collapsed-bed case is already the
# step-6 reference run, and the engineering question is about accumulation
# above the operating charge, not below it.
# 38× (3610 g) is the practical ceiling: the retort's usable volume at
# pf_static = 0.50 holds at most 3653 g of bare kernel, so 38× leaves only a
# ~43 g margin before the bed overflows the vacuum boundary.  Any multiplier
# above 38× (e.g. 40× = 3800 g) would be auto-skipped by vessel_capacity_g.
# See docs/steps/step-7-mass-packing-sweeps.md §Packing-fraction change for
# why pf is 0.50 rather than 0.60.
_MASS_MULTIPLIERS = (2.0, 3.0, 5.0, 10.0, 20.0, 38.0)

_SEED_CHECK_CASES = [
    {'tag': 'seed_42', 'overrides': {'model': {'seed': 42}}},
    {'tag': 'seed_43', 'overrides': {'model': {'seed': 43}}},
    {'tag': 'seed_44', 'overrides': {'model': {'seed': 44}}},
]


def _mass_sweep_cases(params) -> list[dict]:
    """Build the mass sweep case list from _MASS_MULTIPLIERS × nominal charge."""
    nominal_g = float(params['dimensions']['bed']['charge_mass_g'])
    cases = []
    for mult in _MASS_MULTIPLIERS:
        mass_g = round(nominal_g * mult, 3)
        tag = f'mass_{mult:g}x_{mass_g:g}g'.replace('.', 'p')
        cases.append({
            'tag': tag,
            'multiplier': mult,
            'charge_mass_g': mass_g,
            'overrides': {
                'dimensions': {'bed': {'charge_mass_g': mass_g}},
                '_state': 'collapsed',
                '_stage': 'bare_kernel',
                '_background': 'gas',
            },
        })
    return cases


# ---------------------------------------------------------------------------
# Worker (module-level for ProcessPoolExecutor pickling)
# ---------------------------------------------------------------------------

def _worker(args: tuple) -> dict:
    """Run one case in a subprocess; reloads params from disk to avoid pickling issues."""
    params_path, overrides, tag, run_dir_str, csv_path_str, threads, mpi_args, quick, force = args
    base_params = load_params(params_path)
    return run_case(
        base_params, overrides, tag, Path(run_dir_str),
        csv_path=Path(csv_path_str),
        threads=threads,
        mpi_args=mpi_args,
        quick=quick,
        force=force,
    )


# ---------------------------------------------------------------------------
# Sweep definitions
# ---------------------------------------------------------------------------

def seed_check_sweep(params, args) -> list[dict]:
    """3-case sweep varying only seed; validates the pipeline and derived ratios."""
    params_path = str(_REPO_ROOT / 'params.yaml')
    cases = [
        (
            params_path,
            c['overrides'],
            c['tag'],
            str(_SEED_CHECK_DIR / c['tag']),
            str(_SEED_CHECK_CSV),
            args.threads,
            None,           # mpi_args: not implemented for seed check
            args.quick,
            args.force,
        )
        for c in _SEED_CHECK_CASES
    ]

    results = []
    if args.jobs == 1:
        for case_args in cases:
            results.append(_worker(case_args))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(_worker, a): a[2] for a in cases}
            for fut in as_completed(futures):
                tag = futures[fut]
                try:
                    results.append(fut.result())
                except Exception as exc:
                    print(f'  [ERROR] {tag}: {exc}')
        # Restore original order
        order = {c['tag']: i for i, c in enumerate(_SEED_CHECK_CASES)}
        results.sort(key=lambda r: order.get(r.get('tag', ''), 999))

    return results


def mass_sweep(params, args) -> tuple[list[dict], list[dict]]:
    """Mass sweep at collapsed packing, gas background, bare_kernel stage.

    Returns (results, skipped) where `skipped` is the list of cases dropped
    because the charge would exceed vessel capacity at collapsed packing.
    """
    params_path = str(_REPO_ROOT / 'params.yaml')
    cases = _mass_sweep_cases(params)

    max_mass_g = vessel_capacity_g(params, state='collapsed', stage='bare_kernel')
    print(f'  Vessel capacity at collapsed packing (bare_kernel): {max_mass_g:.1f} g')

    to_run, skipped = [], []
    for c in cases:
        if c['charge_mass_g'] > max_mass_g:
            skipped.append(c)
        else:
            to_run.append(c)

    if skipped:
        print(f'  Skipping {len(skipped)} case(s) that exceed vessel capacity:')
        for c in skipped:
            print(f'    {c["tag"]}: {c["charge_mass_g"]:.1f} g > {max_mass_g:.1f} g')

    print(f'  Running {len(to_run)} case(s): '
          f'{[c["tag"] for c in to_run]}')

    worker_cases = [
        (
            params_path,
            c['overrides'],
            c['tag'],
            str(_MASS_SWEEP_DIR / c['tag']),
            str(_MASS_SWEEP_CSV),
            args.threads,
            None,
            args.quick,
            args.force,
        )
        for c in to_run
    ]

    results = []
    if args.jobs == 1:
        for wc in worker_cases:
            results.append(_worker(wc))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(_worker, a): a[2] for a in worker_cases}
            for fut in as_completed(futures):
                tag = futures[fut]
                try:
                    results.append(fut.result())
                except Exception as exc:
                    print(f'  [ERROR] {tag}: {exc}')
        order = {c['tag']: i for i, c in enumerate(to_run)}
        results.sort(key=lambda r: order.get(r.get('tag', ''), 999))

    return results, skipped


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _print_csv(csv_path: Path) -> None:
    if not csv_path.exists():
        print(f'  (no CSV at {csv_path})')
        return
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print('  (empty CSV)')
        return
    # Print a readable subset of columns
    show = ['tag', 'seed', 'u235_mass_g', 'k_eff', 'sigma', 'k_plus_2sigma',
            'k_plus_2sigma_plus_dk_disc', 'bed_height_cm', 'pf_achieved',
            'h_per_u235', 'c_per_u235', 'thermal_flux_fraction', 'wall_time_s']
    available = [c for c in show if c in rows[0]]
    col_w = 22
    header = '  ' + '  '.join(f'{c:<{col_w}}' for c in available)
    print(header)
    print('  ' + '-' * (len(header) - 2))
    for row in rows:
        line = '  ' + '  '.join(f'{row.get(c, ""):<{col_w}}' for c in available)
        print(line)


def main() -> None:
    parser = argparse.ArgumentParser(description='CVD furnace NCS parametric sweep driver')
    parser.add_argument('--seed-check', action='store_true',
                        help='run 3-case seed-variation sweep to validate the pipeline')
    parser.add_argument('--mass-sweep', action='store_true',
                        help='run step-7 mass sweep (collapsed packing, gas, bare kernel)')
    parser.add_argument('--quick', action='store_true',
                        help='use reduced particles/batches (1 000 particles, 12 batches)')
    parser.add_argument('--force', action='store_true',
                        help='re-run cases even if tag already exists in CSV')
    parser.add_argument('-s', '--threads', type=int, default=None, metavar='N',
                        help='OpenMC shared-memory threads per case (default: OMP_NUM_THREADS)')
    parser.add_argument('-j', '--jobs', type=int, default=1, metavar='N',
                        help='number of cases to run in parallel (default: 1)')
    parser.add_argument('--csv', type=str, default=None,
                        help='path to results CSV (overrides per-sweep default)')
    args = parser.parse_args()

    check_env()
    params = load_params()

    if not (args.seed_check or args.mass_sweep):
        parser.print_help()
        print('\nNo sweep selected. Pass --seed-check or --mass-sweep.')
        return

    if args.seed_check:
        mode = 'quick' if args.quick else 'full'
        print(f'\nSeed-check sweep ({mode} mode, {args.jobs} job(s))')
        print(f'  CSV → {_SEED_CHECK_CSV}')
        print(f'  Cases: {[c["tag"] for c in _SEED_CHECK_CASES]}')

        results = seed_check_sweep(params, args)

        print(f'\nResults ({len(results)} cases):')
        _print_csv(_SEED_CHECK_CSV)
        print(f'\nCSV written to: {_SEED_CHECK_CSV}')

    if args.mass_sweep:
        mode = 'quick' if args.quick else 'full'
        print(f'\nMass sweep ({mode} mode, {args.jobs} job(s))')
        print(f'  CSV → {_MASS_SWEEP_CSV}')

        results, skipped = mass_sweep(params, args)

        print(f'\nResults ({len(results)} cases run, {len(skipped)} skipped):')
        _print_csv(_MASS_SWEEP_CSV)
        print(f'\nCSV written to: {_MASS_SWEEP_CSV}')
        print('Generate plot with:  python scripts/plot_mass_sweep.py')


if __name__ == '__main__':
    main()
