#!/usr/bin/env python
"""Step 8 flooding accident sweep driver — bottom-up flood from injector throat.

The injector coolant leak is the bounding water source for this design: water
enters at the cone throat and fills upward.  A partial fill from the bottom is
more reactive than full flood when the bed is collapsed, because the cone is
where mass concentrates — moderator arrives exactly where the fuel is densest.

Sweep parameters:
    - Bed state   : collapsed (pf_static), fluidized (pf_fluidized)
    - Water density: liquid (1.0 g/cm³), vapor (0.001 g/cm³ — near-atmospheric steam)
    - z_flood     : 20 levels from throat to retort top

Case count: 2 states × 20 levels × 2 densities = 80 cases total
Dry baselines (no water) were already run in step 5 (run_nominal.py / run_collapsed.py).
CSV → results/flood_bottomup.csv

Usage:
    # Smoke-test (1 000 particles, 12 batches, ~30 s/case):
    caffeinate python scripts/run_flood.py --quick

    # Full run:
    caffeinate python scripts/run_flood.py

    # Parallel (N cases concurrently):
    caffeinate python scripts/run_flood.py -j 4

    # Re-run cases already in CSV:
    caffeinate python scripts/run_flood.py --force

    # Plot after CSV is written:
    caffeinate python scripts/plot_flood_sweep.py
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
from furnace.sweeps import run_case

_RESULTS_DIR = _REPO_ROOT / 'results'
_BU_CSV      = _RESULTS_DIR / 'flood_bottomup.csv'
_BU_DIR      = _RESULTS_DIR / 'flood_bottomup'

_STATES = ['collapsed', 'fluidized']

_WATER_DENSITIES = {
    'liquid': 1.0,
    'vapor':  0.001,  # near-atmospheric steam (g/cm³)
}


# ---------------------------------------------------------------------------
# Case definitions
# ---------------------------------------------------------------------------

def _bottomup_cases(params) -> list[dict]:
    """2 states × 20 z_flood levels × 2 densities = 80 cases.

    Dry baselines are omitted — those were already run in step 5.
    z_flood levels are computed from params so they stay in sync with any
    geometry edits to params.yaml.
    """
    dim = params['dimensions']
    z_rt = dim['cone']['vertical_drop_cm'] + dim['retort']['height_cm']
    n_levels = 20
    levels = [round(z_rt * (i + 1) / n_levels, 3) for i in range(n_levels)]

    cases = []
    for state in _STATES:
        base = {'_state': state, '_stage': 'bare_kernel'}
        for density_label, density_gcc in _WATER_DENSITIES.items():
            for i, z in enumerate(levels):
                cases.append({
                    'tag': f'{state}_bottomup_{density_label}_level{i+1:02d}_z{z:.3f}cm',
                    'overrides': {
                        **base,
                        '_background': 'gas',
                        '_flood_extent': 'none',
                        '_z_flood': z,
                        '_water_density': density_gcc,
                    },
                })

    return cases


# ---------------------------------------------------------------------------
# Worker (module-level for ProcessPoolExecutor pickling)
# ---------------------------------------------------------------------------

def _worker(args: tuple) -> dict:
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


def _run_cases(cases, run_base_dir, csv_path, params_path, args) -> list[dict]:
    worker_args = [
        (
            params_path,
            c['overrides'],
            c['tag'],
            str(run_base_dir / c['tag']),
            str(csv_path),
            args.threads,
            None,
            args.quick,
            args.force,
        )
        for c in cases
    ]

    results = []
    if args.jobs == 1:
        for wa in worker_args:
            results.append(_worker(wa))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(_worker, a): a[2] for a in worker_args}
            for fut in as_completed(futures):
                tag = futures[fut]
                try:
                    results.append(fut.result())
                except Exception as exc:
                    print(f'  [ERROR] {tag}: {exc}')
        order = {c['tag']: i for i, c in enumerate(cases)}
        results.sort(key=lambda r: order.get(r.get('tag', ''), 999))

    return results


# ---------------------------------------------------------------------------
# CLI helpers
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
    show = ['tag', 'state', 'z_flood_cm', 'water_density_gcc',
            'k_eff', 'sigma', 'k_plus_2sigma', 'k_plus_2sigma_plus_dk_disc',
            'h_per_u235', 'wall_time_s']
    available = [c for c in show if c in rows[0]]
    col_w = 22
    header = '  ' + '  '.join(f'{c:<{col_w}}' for c in available)
    print(header)
    print('  ' + '-' * (len(header) - 2))
    for row in rows:
        line = '  ' + '  '.join(f'{row.get(c, ""):<{col_w}}' for c in available)
        print(line)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='CVD furnace step-8 bottom-up flood sweep driver'
    )
    parser.add_argument('--quick', action='store_true',
                        help='use reduced particles/batches (1 000 particles, 12 batches)')
    parser.add_argument('--force', action='store_true',
                        help='re-run cases even if tag already exists in CSV')
    parser.add_argument('-s', '--threads', type=int, default=None, metavar='N',
                        help='OpenMC shared-memory threads per case')
    parser.add_argument('-j', '--jobs', type=int, default=1, metavar='N',
                        help='number of cases to run in parallel')
    args = parser.parse_args()

    check_env()
    params = load_params()
    params_path = str(_REPO_ROOT / 'params.yaml')
    mode = 'quick' if args.quick else 'full'

    cases = _bottomup_cases(params)
    print(f'\nBottom-up flood sweep ({mode} mode, {args.jobs} job(s))')
    print(f'  Bed states     : {_STATES}')
    print(f'  Water densities: {_WATER_DENSITIES}')
    print(f'  Cases          : {len(cases)} (dry baselines run in step 5)')
    print(f'  CSV           → {_BU_CSV}')

    results = _run_cases(cases, _BU_DIR, _BU_CSV, params_path, args)

    print(f'\nResults ({len(results)} cases):')
    _print_csv(_BU_CSV)
    print(f'\nCSV written to: {_BU_CSV}')
    print('\nGenerate plots with:  caffeinate python scripts/plot_flood_sweep.py')


if __name__ == '__main__':
    main()
