"""Step 5 convergence plots: Shannon entropy, batch comparison, seed comparison.

Reads statepoints from step5_nominal/ and step5_convergence/ and writes
PNG plots to results/step5_convergence/plots/.

Run AFTER both run_nominal.py and run_step5_convergence.py:
    python scripts/plot_step5_convergence.py
"""
from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use('Agg')  # non-interactive; safe for headless environments
import matplotlib.pyplot as plt
import numpy as np
import openmc

from furnace.params import load_params

_NOM_DIR   = Path('results/step5_nominal')
_CONV_DIR  = Path('results/step5_convergence')
_PLOT_DIR  = _CONV_DIR / 'plots'


def _load_sp(path: Path) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Return (k_generation, entropy, k_nom, k_sig) from a statepoint."""
    with openmc.StatePoint(str(path)) as sp:
        k_gen   = np.array(sp.k_generation)
        entropy = np.array(sp.entropy)
        k_nom   = sp.keff.nominal_value
        k_sig   = sp.keff.std_dev
    return k_gen, entropy, k_nom, k_sig


def plot_entropy(params: object) -> None:
    """Shannon entropy vs generation from the 500-active batch run."""
    n_inactive = int(params['model']['inactive'])
    sp_path    = _CONV_DIR / 'batch_500' / f'statepoint.{n_inactive + 500}.h5'
    if not sp_path.exists():
        print(f'[entropy] statepoint not found: {sp_path}')
        return

    _, entropy, _, _ = _load_sp(sp_path)
    n_total = len(entropy)
    gens    = np.arange(1, n_total + 1)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(gens, entropy, lw=0.8, color='steelblue', label='Shannon entropy H')
    ax.axvline(n_inactive, color='red', lw=1.2, ls='--', label=f'Inactive→Active (gen {n_inactive})')

    # Mark apparent convergence: first generation where |H - H_tail| < 0.02 H_tail
    H_tail = np.mean(entropy[n_inactive:])
    converged_gen = None
    for i in range(n_inactive):
        if abs(entropy[i] - H_tail) < 0.02 * H_tail:
            converged_gen = i + 1
            break
    if converged_gen:
        ax.axvline(converged_gen, color='green', lw=1.0, ls=':', label=f'≈Converged gen {converged_gen}')

    ax.set_xlabel('Generation')
    ax.set_ylabel('Shannon entropy H')
    ax.set_title('Source convergence: Shannon entropy vs generation\n(seed=42, 500 active batches)')
    ax.legend(fontsize=9)
    ax.grid(True, lw=0.3)
    fig.tight_layout()

    out = _PLOT_DIR / 'entropy_vs_generation.png'
    fig.savefig(str(out), dpi=150)
    plt.close(fig)
    print(f'  Entropy plot → {out}')
    if converged_gen:
        print(f'    Apparent convergence at generation {converged_gen} '
              f'(|H - H_tail| < 2% H_tail)')
    print(f'    H_tail (active batch mean) = {H_tail:.4f}')


def plot_batch_comparison(params: object) -> None:
    """k-eff vs active-batch count: 250 and 500 active at seed=42."""
    n_inactive = int(params['model']['inactive'])
    nom_sp  = _NOM_DIR  / f'statepoint.{n_inactive + 250}.h5'
    b500_sp = _CONV_DIR / 'batch_500' / f'statepoint.{n_inactive + 500}.h5'

    missing = [p for p in (nom_sp, b500_sp) if not p.exists()]
    if missing:
        print(f'[batch comparison] missing statepoints: {missing}')
        return

    # k-eff from final combined estimate in each statepoint
    _, _, k250, s250 = _load_sp(nom_sp)
    _, _, k500, s500 = _load_sp(b500_sp)

    # For the 500-active statepoint, also extract the 250-active sub-estimate
    k_gen_500, _, _, _ = _load_sp(b500_sp)
    k_sub250 = np.mean(k_gen_500[n_inactive : n_inactive + 250])
    s_sub250 = np.std( k_gen_500[n_inactive : n_inactive + 250], ddof=1) / np.sqrt(250)

    n_active_vals = [250, 500]
    k_vals = [k250, k500]
    s_vals = [s250, s500]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(n_active_vals, k_vals, yerr=[2*s for s in s_vals],
                fmt='o-', capsize=5, color='steelblue', label='k-eff ± 2σ')
    ax.set_xlabel('Active batches')
    ax.set_ylabel('k-eff')
    ax.set_title('Batch-count convergence check (seed=42)\n20 000 particles/gen, 50 inactive')
    ax.legend(fontsize=9)
    ax.grid(True, lw=0.3)

    # Annotate values
    for n, k, s in zip(n_active_vals, k_vals, s_vals):
        ax.annotate(f'{k:.4f}±{s:.4f}', (n, k), textcoords='offset points',
                    xytext=(0, 10), ha='center', fontsize=8)

    fig.tight_layout()
    out = _PLOT_DIR / 'keff_batch_comparison.png'
    fig.savefig(str(out), dpi=150)
    plt.close(fig)

    diff = abs(k250 - k500)
    comb = (s250**2 + s500**2)**0.5
    print(f'  Batch comparison → {out}')
    print(f'    k(250) = {k250:.5f} ± {s250:.5f}')
    print(f'    k(500) = {k500:.5f} ± {s500:.5f}')
    print(f'    |Δk| = {diff:.5f}   combined σ = {comb:.5f}   '
          f'{"PASS" if diff < comb else "FAIL"}')


def plot_seed_comparison(params: object) -> None:
    """k-eff across seeds 42, 43, 44 at 250 active batches."""
    n_inactive = int(params['model']['inactive'])
    statepoints = {
        42: _NOM_DIR  / f'statepoint.{n_inactive + 250}.h5',
        43: _CONV_DIR / 'seed_43' / f'statepoint.{n_inactive + 250}.h5',
        44: _CONV_DIR / 'seed_44' / f'statepoint.{n_inactive + 250}.h5',
    }

    missing = [p for p in statepoints.values() if not p.exists()]
    if missing:
        print(f'[seed comparison] missing statepoints: {missing}')
        return

    seeds, k_vals, s_vals = [], [], []
    for seed, sp_path in statepoints.items():
        _, _, k, s = _load_sp(sp_path)
        seeds.append(seed)
        k_vals.append(k)
        s_vals.append(s)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(seeds, k_vals, yerr=[2*s for s in s_vals],
                fmt='o', capsize=6, markersize=7, color='darkorange', label='k-eff ± 2σ')
    ax.axhline(np.mean(k_vals), ls='--', color='gray', lw=1, label='Mean k-eff')
    ax.set_xticks(seeds)
    ax.set_xlabel('RNG seed (geometry + transport)')
    ax.set_ylabel('k-eff')
    ax.set_title('Seed-to-seed reproducibility (3 independent geometries)\n250 active batches, 20 000 particles/gen')
    ax.legend(fontsize=9)
    ax.grid(True, lw=0.3, axis='y')
    for seed, k, s in zip(seeds, k_vals, s_vals):
        ax.annotate(f'{k:.4f}±{s:.4f}', (seed, k), textcoords='offset points',
                    xytext=(0, 12), ha='center', fontsize=8)
    fig.tight_layout()

    out = _PLOT_DIR / 'keff_seed_comparison.png'
    fig.savefig(str(out), dpi=150)
    plt.close(fig)
    print(f'  Seed comparison → {out}')


def main() -> None:
    params = load_params()
    _PLOT_DIR.mkdir(parents=True, exist_ok=True)

    print('Generating step 5 convergence plots...')
    plot_entropy(params)
    plot_batch_comparison(params)
    plot_seed_comparison(params)
    print('Done.')


if __name__ == '__main__':
    main()
