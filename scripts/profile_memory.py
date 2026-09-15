"""Memory profiler for the bed geometry pipeline.

Runs ``build_model`` at a series of charge masses in isolated subprocesses,
measuring peak RSS at four checkpoints:

    baseline           — interpreter start, after openmc import
    after_params       — after check_env() + load_params()
    after_pack         — after every pack_bed() call has returned
    after_lattice      — after every lattice_bed() call has returned
    after_build_model  — after build_model() returns (geometry finalised)

The point of this diagnostic is to answer: at high mass (many TRISO particles),
is memory dominated by (a) the per-particle openmc.Cell objects that pack_bed +
lattice_bed create, or by (b) the RectLattice / XML output layer? The answer
determines whether a tiled-unit-cell refactor is worth the effort.

Each mass point is run in its own subprocess so the parent doesn't accumulate
memory between points. Requires psutil.

Usage:
    python3 scripts/profile_memory.py                          # default sweep
    python3 scripts/profile_memory.py --masses 95 190 475 950  # custom sweep
    python3 scripts/profile_memory.py --child --mass 190       # (internal)
"""
from __future__ import annotations

import argparse
import gc
import json
import subprocess
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))


_DEFAULT_MASSES = [95.0, 190.0, 475.0, 950.0]  # 1×, 2×, 5×, 10× nominal


def _rss_mb() -> float:
    import psutil
    return psutil.Process().memory_info().rss / (1024 ** 2)


def _run_child(mass_g: float) -> dict:
    """In-process: build model at *mass_g*, return checkpoint dict."""
    import openmc  # noqa: F401 — ensure openmc allocations counted in baseline
    baseline = _rss_mb()

    from furnace.params import check_env, load_params
    from furnace import triso as _triso_mod
    from furnace import geometry as _geom_mod

    check_env()
    params = load_params()
    after_params = _rss_mb()

    # Wrap pack_bed / lattice_bed to record max RSS observed after each return.
    counters = {'after_pack': 0.0, 'after_lattice': 0.0, 'n_particles': 0}

    _orig_pack = _triso_mod.pack_bed
    _orig_lat  = _triso_mod.lattice_bed

    def pack_bed_traced(*args, **kwargs):
        trisos = _orig_pack(*args, **kwargs)
        counters['n_particles'] += len(trisos)
        counters['after_pack'] = max(counters['after_pack'], _rss_mb())
        return trisos

    def lattice_bed_traced(*args, **kwargs):
        lat = _orig_lat(*args, **kwargs)
        counters['after_lattice'] = max(counters['after_lattice'], _rss_mb())
        return lat

    _triso_mod.pack_bed = pack_bed_traced
    _triso_mod.lattice_bed = lattice_bed_traced
    # geometry.py imported these by name at module load, so patch there too.
    _geom_mod.pack_bed = pack_bed_traced
    _geom_mod.lattice_bed = lattice_bed_traced

    from furnace.model import build_model

    model, stats = build_model(
        params,
        state='collapsed',
        stage='bare_kernel',
        charge_mass_g=mass_g,
        n_inactive=1,
        n_active=1,
        n_particles=1_000,
    )
    gc.collect()
    after_build = _rss_mb()

    # Count unique openmc.Cell instances actually registered in the geometry.
    n_cells = len(model.geometry.get_all_cells())

    return {
        'mass_g': mass_g,
        'n_particles': counters['n_particles'],
        'n_cells': n_cells,
        'baseline_mb': baseline,
        'after_params_mb': after_params,
        'after_pack_mb': counters['after_pack'],
        'after_lattice_mb': counters['after_lattice'],
        'after_build_mb': after_build,
    }


def _child_main(mass_g: float) -> None:
    result = _run_child(mass_g)
    print('__RESULT__' + json.dumps(result))


def _parent_main(masses: list[float]) -> None:
    rows = []
    for m in masses:
        print(f'\n── running child at {m:g} g ──', flush=True)
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), '--child', '--mass', str(m)],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            print(f'  FAILED (rc={proc.returncode})')
            print('  stderr:', proc.stderr[-2000:])
            continue
        for line in proc.stdout.splitlines():
            if line.startswith('__RESULT__'):
                rows.append(json.loads(line[len('__RESULT__'):]))
                break

    if not rows:
        print('\nNo successful runs.')
        return

    print('\n' + '=' * 96)
    print(f"{'mass_g':>8} {'N_particles':>12} {'N_cells':>10} "
          f"{'base':>7} {'params':>7} {'pack':>7} {'lattice':>8} {'build':>7} "
          f"{'Δpack':>7} {'Δlat':>7}")
    print('-' * 96)
    for r in rows:
        d_pack = r['after_pack_mb'] - r['after_params_mb']
        d_lat  = r['after_lattice_mb'] - r['after_pack_mb']
        print(f"{r['mass_g']:>8.1f} {r['n_particles']:>12d} {r['n_cells']:>10d} "
              f"{r['baseline_mb']:>7.1f} {r['after_params_mb']:>7.1f} "
              f"{r['after_pack_mb']:>7.1f} {r['after_lattice_mb']:>8.1f} "
              f"{r['after_build_mb']:>7.1f} "
              f"{d_pack:>7.1f} {d_lat:>7.1f}")
    print('=' * 96)
    print('\nAll values in MB (RSS). Δpack = pack − params. Δlat = lattice − pack.')
    print('If Δpack scales linearly with N_particles → per-particle Cell objects dominate.')
    print('If Δlat scales linearly instead → lattice/XML layer dominates.')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--child', action='store_true', help=argparse.SUPPRESS)
    ap.add_argument('--mass', type=float, help='(child-only) single mass in g')
    ap.add_argument('--masses', type=float, nargs='+', default=_DEFAULT_MASSES,
                    help='charge masses in g to profile (parent mode)')
    args = ap.parse_args()

    if args.child:
        if args.mass is None:
            ap.error('--child requires --mass')
        _child_main(args.mass)
    else:
        _parent_main(args.masses)


if __name__ == '__main__':
    main()
