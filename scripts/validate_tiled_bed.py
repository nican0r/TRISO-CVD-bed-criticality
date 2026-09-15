"""Validate the tiled bed against the fully-random pack_bed reference.

Runs a k-eff calculation at the same charge mass in both modes:

    random  : fully random pack_spheres per slab + create_triso_lattice
    tiled   : one small cube packed once, tiled as a shared-universe RectLattice

Reports k-eff, σ(k), Δk = k_tiled − k_random, and peak Python RSS during build.
If |Δk| is within a couple σ (few × 10⁻³ at the default statistics), the tile
periodicity is not physically distorting the answer and the tiled path is safe
to use for the mass sweep.

Runs each config in a subprocess so RSS peaks are cleanly isolated.

Usage:
    python3 scripts/validate_tiled_bed.py                        # 95 g, quick
    python3 scripts/validate_tiled_bed.py --mass 190             # 2× nominal
    python3 scripts/validate_tiled_bed.py --particles 20000 \\    # tighter stats
        --inactive 50 --active 200
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


def _child(mass_g: float, mode: str, tile_size_cm: float,
           n_particles: int, n_inactive: int, n_active: int, seed: int,
           use_exact_cone: bool = False) -> dict:
    import psutil
    p = psutil.Process()

    from furnace.params import check_env, load_params
    from furnace.model import build_model, export_and_run

    check_env()
    params = load_params()
    rss_before = p.memory_info().rss / (1024 ** 2)

    use_tiled = (mode == 'tiled')
    model, stats = build_model(
        params,
        state='collapsed',
        stage='bare_kernel',
        charge_mass_g=mass_g,
        n_inactive=n_inactive,
        n_active=n_active,
        n_particles=n_particles,
        seed=seed,
        use_tiled_bed=use_tiled,
        tile_size_cm=tile_size_cm,
        use_exact_cone=use_exact_cone,
    )
    gc.collect()
    rss_build = p.memory_info().rss / (1024 ** 2)
    n_cells = len(model.geometry.get_all_cells())

    geom_tag = 'exactcone' if use_exact_cone else 'staircase'
    out_dir = _REPO_ROOT / 'results' / 'validate_tiled_bed' / f'{geom_tag}_{mode}_{int(mass_g)}g'
    out_dir.mkdir(parents=True, exist_ok=True)

    # Reuse a completed statepoint if present (avoids re-running the 20-min transport
    # after a post-processing crash). Matches on batch count in the filename.
    sp_pattern = f'statepoint.{n_inactive + n_active}.h5'
    existing = list(out_dir.glob(sp_pattern))
    if existing:
        sp_path = existing[0]
        print(f'  reusing {sp_path.name} (skipping transport)', flush=True)
    else:
        sp_path = export_and_run(model, out_dir)

    # autolink=False skips summary.h5 — that file scales with unique-Cell count and
    # takes minutes to parse in random mode. keff is stored directly in the statepoint.
    import openmc
    with openmc.StatePoint(sp_path, autolink=False) as sp:
        k_combined = sp.keff
        k = float(k_combined.nominal_value)
        sig = float(k_combined.std_dev)

    return {
        'mode': mode,
        'mass_g': mass_g,
        'tile_size_cm': tile_size_cm if use_tiled else None,
        'n_particles': stats.n_particles,
        'pf_achieved': stats.pf_achieved,
        'bed_height_cm': stats.bed_height_cm,
        'u235_mass_g': stats.u235_mass_g,
        'n_unique_cells': n_cells,
        'rss_before_mb': rss_before,
        'rss_after_build_mb': rss_build,
        'k': k,
        'sigma': sig,
    }


def _child_main(args) -> None:
    result = _child(args.mass, args.mode, args.tile_size,
                    args.particles, args.inactive, args.active, args.seed,
                    use_exact_cone=args.exact_cone)
    print('__RESULT__' + json.dumps(result))


def _parent_main(args) -> None:
    rows = []
    for mode in ('random', 'tiled'):
        print(f'\n── running {mode} at {args.mass:g} g '
              f'({args.particles} p/gen, {args.inactive} inactive + {args.active} active) ──',
              flush=True)
        child_argv = [sys.executable, str(Path(__file__).resolve()),
                      '--child', '--mode', mode,
                      '--mass', str(args.mass),
                      '--tile-size', str(args.tile_size),
                      '--particles', str(args.particles),
                      '--inactive', str(args.inactive),
                      '--active', str(args.active),
                      '--seed', str(args.seed)]
        if args.exact_cone:
            child_argv.append('--exact-cone')
        proc = subprocess.run(child_argv, capture_output=True, text=True)
        if proc.returncode != 0:
            print(f'  FAILED (rc={proc.returncode})')
            print('  stderr tail:', proc.stderr[-2000:])
            return
        for line in proc.stdout.splitlines():
            if line.startswith('__RESULT__'):
                rows.append(json.loads(line[len('__RESULT__'):]))
                break

    if len(rows) != 2:
        print('Missing results')
        return
    r_rand = next(r for r in rows if r['mode'] == 'random')
    r_tile = next(r for r in rows if r['mode'] == 'tiled')

    dk = r_tile['k'] - r_rand['k']
    sig_diff = (r_rand['sigma'] ** 2 + r_tile['sigma'] ** 2) ** 0.5
    z = dk / sig_diff if sig_diff > 0 else float('inf')

    geom_label = 'exact-cone' if args.exact_cone else 'staircase'
    print('\n' + '=' * 78)
    print(f"  Validation: tiled vs random {geom_label} bed at charge_mass = {args.mass:g} g")
    print('=' * 78)
    print(f"{'':>18} {'random':>18} {'tiled':>18}")
    for key, fmt in [
        ('n_particles',       '{:>18d}'),
        ('pf_achieved',       '{:>18.4f}'),
        ('bed_height_cm',     '{:>18.4f}'),
        ('n_unique_cells',    '{:>18d}'),
        ('rss_after_build_mb','{:>18.1f}'),
        ('k',                 '{:>18.5f}'),
        ('sigma',             '{:>18.5f}'),
    ]:
        vr, vt = r_rand[key], r_tile[key]
        print(f"  {key:<16} " + fmt.format(vr) + ' ' + fmt.format(vt))
    print('-' * 78)
    print(f"  Δk  (tiled − random) : {dk:+.5f}")
    print(f"  σ(Δk)                : {sig_diff:.5f}")
    print(f"  |Δk| / σ(Δk)         : {abs(z):.2f}")
    print('=' * 78)
    if abs(z) < 3:
        print('  → within 3σ: tiled path OK to use for the mass sweep.')
    else:
        print('  → >3σ deviation: tile periodicity is affecting k-eff. Investigate.')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--mass', type=float, default=95.0,
                    help='charge mass in g (default 95 = 1× nominal)')
    ap.add_argument('--tile-size', type=float, default=0.5,
                    help='tiled-bed unit-cell edge in cm (default 0.5)')
    ap.add_argument('--particles', type=int, default=10_000,
                    help='particles per generation (default 10 000)')
    ap.add_argument('--inactive', type=int, default=30,
                    help='inactive batches (default 30)')
    ap.add_argument('--active', type=int, default=100,
                    help='active batches (default 100)')
    ap.add_argument('--seed', type=int, default=42, help='OpenMC seed (default 42)')
    ap.add_argument('--exact-cone', action='store_true',
                    help='validate the exact-cone reference geometry instead of the staircase')
    ap.add_argument('--child', action='store_true', help=argparse.SUPPRESS)
    ap.add_argument('--mode', choices=('random', 'tiled'), help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.child:
        if args.mode is None:
            ap.error('--child requires --mode')
        _child_main(args)
    else:
        _parent_main(args)


if __name__ == '__main__':
    main()
