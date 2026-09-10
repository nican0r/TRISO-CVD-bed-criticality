"""Step 7 mass sweep plot: k+2σ vs charge mass with 0.95 subcritical line.

Reads results/mass_sweep.csv (written by run_sweep.py --mass-sweep) and writes
PNG to results/mass_sweep.png.

Run AFTER the sweep completes:
    caffeinate python scripts/plot_mass_sweep.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_CSV = _REPO_ROOT / 'results' / 'mass_sweep.csv'
_PNG = _REPO_ROOT / 'results' / 'mass_sweep.png'


def _load_rows(csv_path: Path) -> list[dict]:
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r['u235_mass_g']              = float(r['u235_mass_g'])
        r['k_eff']                    = float(r['k_eff'])
        r['sigma']                    = float(r['sigma'])
        r['k_plus_2sigma']            = float(r['k_plus_2sigma'])
        r['k_plus_2sigma_plus_dk_disc'] = float(r['k_plus_2sigma_plus_dk_disc'])
        # charge mass is stored inside overrides_json
        import json as _json
        ov = _json.loads(r.get('overrides_json', '{}'))
        r['charge_mass_g'] = float(
            ov.get('dimensions', {}).get('bed', {}).get('charge_mass_g', 0.0)
        )
    rows.sort(key=lambda r: r['charge_mass_g'])
    return rows


def plot() -> None:
    if not _CSV.exists():
        print(f'CSV not found: {_CSV}')
        print('Run:  caffeinate python scripts/run_sweep.py --mass-sweep')
        return

    rows = _load_rows(_CSV)
    if not rows:
        print(f'No rows in {_CSV}')
        return

    mass  = [r['charge_mass_g'] for r in rows]
    k     = [r['k_eff'] for r in rows]
    sig   = [r['sigma'] for r in rows]
    k2s   = [r['k_plus_2sigma'] for r in rows]
    k2sd  = [r['k_plus_2sigma_plus_dk_disc'] for r in rows]
    u235  = [r['u235_mass_g'] for r in rows]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(mass, k, yerr=[2*s for s in sig], fmt='o-', capsize=4,
                color='steelblue', label='k-eff ± 2σ')
    ax.plot(mass, k2sd, 's--', color='firebrick',
            label='k + 2σ + δk$_{disc}$ (conservative)')
    ax.axhline(0.95, color='black', lw=1.2, ls=':', label='0.95 subcritical limit')
    ax.set_xscale('log')
    ax.set_xlabel('Charge mass (g, bare UCO kernels)')
    ax.set_ylabel('k-eff')
    ax.set_title('Mass sweep — collapsed bed, gas atmosphere, bare kernels\n'
                 '(UNVALIDATED — SCREENING ONLY)')
    ax.grid(True, which='both', lw=0.3)
    ax.legend(fontsize=9, loc='best')

    for m, kk, ss, um in zip(mass, k, sig, u235):
        ax.annotate(f'{kk:.3f}±{ss:.3f}\n({um:.1f} g $^{{235}}$U)',
                    (m, kk), textcoords='offset points',
                    xytext=(6, 8), fontsize=7)

    fig.text(0.5, 0.01,
             'UNVALIDATED — SCREENING ONLY. Not benchmarked against ICSBEP.',
             ha='center', fontsize=8, style='italic', color='dimgray')
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    fig.savefig(str(_PNG), dpi=150)
    plt.close(fig)
    print(f'Plot → {_PNG}')

    print('\nSummary table:')
    print(f'  {"charge_g":>10}  {"U235_g":>8}  {"k":>7}  {"±2σ":>7}  {"+δk_disc":>9}')
    for m, kk, ss, k2, k2d, um in zip(mass, k, sig, k2s, k2sd, u235):
        print(f'  {m:>10.2f}  {um:>8.2f}  {kk:>7.4f}  {ss:>7.4f}  {k2d:>9.4f}')


if __name__ == '__main__':
    plot()
