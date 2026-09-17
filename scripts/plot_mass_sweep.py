"""Step 7 mass sweep plot: k+2σ vs charge mass with 0.95 subcritical line.

Reads either
  * results/mass_sweep.csv (written by run_sweep.py --mass-sweep, local runs), or
  * results/<sweep>/<submission>/summary.csv (written by scripts/batch/pull_results.py)

Auto-detects the schema by column names. Pass --csv to point at a specific file
and --png to control the output path.

Run AFTER the sweep completes:
    caffeinate python scripts/plot_mass_sweep.py
    python scripts/plot_mass_sweep.py --csv results/step7a_mass_tiled/<sub>/summary.csv
"""
from __future__ import annotations

import argparse
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
    if not rows:
        return []
    # Detect schema by column presence.
    is_batch = 'keff' in rows[0] and 'row.charge_mass_g' in rows[0]
    for r in rows:
        if is_batch:
            r['charge_mass_g']              = float(r['row.charge_mass_g'])
            r['u235_mass_g']                = float(r['u235_mass_g'])
            r['k_eff']                      = float(r['keff'])
            r['sigma']                      = float(r['keff_std'])
            r['k_plus_2sigma']              = float(r['keff_plus_2sigma'])
            # Batch results don't carry a δk_disc column; conservative bound falls
            # back to k+2σ. Documented as a Stage-0 approximation in step 7 doc.
            r['k_plus_2sigma_plus_dk_disc'] = float(r['keff_plus_2sigma'])
        else:
            r['u235_mass_g']              = float(r['u235_mass_g'])
            r['k_eff']                    = float(r['k_eff'])
            r['sigma']                    = float(r['sigma'])
            r['k_plus_2sigma']            = float(r['k_plus_2sigma'])
            r['k_plus_2sigma_plus_dk_disc'] = float(r['k_plus_2sigma_plus_dk_disc'])
            import json as _json
            ov = _json.loads(r.get('overrides_json', '{}'))
            r['charge_mass_g'] = float(
                ov.get('dimensions', {}).get('bed', {}).get('charge_mass_g', 0.0)
            )
    rows.sort(key=lambda r: r['charge_mass_g'])
    return rows


def plot(csv_path: Path, png_path: Path, ymax: float | None = None,
         show_limit_line: bool = False, unvalidated_title: bool = True,
         dk_disc: float | None = None, show_dk_disc_line: bool = True) -> None:
    if not csv_path.exists():
        print(f'CSV not found: {csv_path}')
        print('Run:  caffeinate python scripts/run_sweep.py --mass-sweep')
        print('  or: python scripts/batch/pull_results.py --sweep <name>')
        return

    rows = _load_rows(csv_path)
    if dk_disc is not None:
        # Batch schema doesn't carry a δk_disc column, so any external estimate
        # (e.g. from the step-3 convergence study) is applied here on top of k+2σ.
        for r in rows:
            r['k_plus_2sigma_plus_dk_disc'] = r['k_plus_2sigma'] + dk_disc
    if not rows:
        print(f'No rows in {csv_path}')
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
    if show_dk_disc_line:
        ax.plot(mass, k2sd, 's--', color='firebrick',
                label='k + 2σ + δk$_{disc}$ (conservative)')
    if show_limit_line:
        ax.axhline(0.95, color='black', lw=1.2, ls=':', label='0.95 subcritical limit')
    ax.set_xscale('log')
    # Widen the x-range beyond the data span so end-point annotations don't get
    # clipped by the plot edge or crossed by the enclosing vertical gridlines.
    # ~0.2 decades of log-space padding on each side.
    if mass:
        pad = 10 ** 0.2
        ax.set_xlim(min(mass) / pad, max(mass) * pad)
    if ymax is not None:
        ax.set_ylim(0, ymax)
    ax.set_xlabel('Charge mass (g, bare UCO kernels)')
    ax.set_ylabel('k-eff')
    title = 'Mass sweep — collapsed bed, gas atmosphere, bare kernels'
    if unvalidated_title:
        title += '\n(UNVALIDATED — SCREENING ONLY)'
    ax.set_title(title)
    ax.grid(True, which='both', lw=0.3)
    ax.legend(fontsize=9, loc='best')

    # Alternate labels above/below the line to avoid overlap on the rising curve.
    # Label shows k to 4 d.p. and 2σ in scientific notation (σ ~ 1e-5, invisible at
    # plot scale, but meaningful to record).
    for i, (m, kk, ss, um) in enumerate(zip(mass, k, sig, u235)):
        above = (i % 2 == 0)
        yoffset = 22 if above else -30
        va = 'bottom' if above else 'top'
        ax.annotate(f'{kk:.4f}±{2*ss:.1e}\n({m:.0f} g / {um:.1f} g $^{{235}}$U)',
                    xy=(m, kk),
                    xytext=(0, yoffset), textcoords='offset points',
                    ha='center', va=va, fontsize=7,
                    arrowprops=dict(arrowstyle='-', lw=0.4, color='gray',
                                    shrinkA=0, shrinkB=2))

    fig.text(0.5, 0.01,
             'UNVALIDATED — SCREENING ONLY. Not benchmarked against ICSBEP.',
             ha='center', fontsize=8, style='italic', color='dimgray')
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    fig.savefig(str(png_path), dpi=150)
    plt.close(fig)
    print(f'Plot → {png_path}')

    print('\nSummary table:')
    print(f'  {"charge_g":>10}  {"U235_g":>8}  {"k":>7}  {"±2σ":>7}  {"+δk_disc":>9}')
    for m, kk, ss, k2, k2d, um in zip(mass, k, sig, k2s, k2sd, u235):
        print(f'  {m:>10.2f}  {um:>8.2f}  {kk:>7.4f}  {ss:>7.4f}  {k2d:>9.4f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', type=Path, default=_CSV,
                    help='Input CSV (auto-detects local vs batch schema).')
    ap.add_argument('--png', type=Path, default=None,
                    help='Output PNG path. Defaults to <csv-dir>/mass_sweep.png.')
    ap.add_argument('--ymax', type=float, default=None,
                    help='Upper y-axis limit for k-eff (bottom pinned at 0). '
                         'Auto if omitted; pass e.g. 0.2 to zoom in on a deeply '
                         'subcritical range.')
    ap.add_argument('--limit-line', action='store_true',
                    help='Show the 0.95 subcritical horizontal reference line (hidden by default).')
    ap.add_argument('--no-unvalidated-title', action='store_true',
                    help='Omit the "(UNVALIDATED — SCREENING ONLY)" subtitle line '
                         '(the footer caveat is untouched).')
    ap.add_argument('--no-dk-disc-line', action='store_true',
                    help='Suppress the k+2σ+δk_disc conservative trace. '
                         'Appropriate when using exact-cone tiled geometry where '
                         'the staircase discretisation bias (δk_disc) is zero.')
    ap.add_argument('--dk-disc', type=float, default=None,
                    help='Discretisation bias (δk_disc) to add on top of k+2σ '
                         'for the conservative trace. Batch summary.csv has no '
                         'δk_disc column; pass e.g. 0.00230 from the step-3 '
                         'exact-cone-vs-n=32 packed-particle convergence study.')
    args = ap.parse_args()
    png = args.png if args.png is not None else args.csv.parent / 'mass_sweep.png'
    plot(args.csv, png, ymax=args.ymax,
         show_limit_line=args.limit_line,
         unvalidated_title=not args.no_unvalidated_title,
         dk_disc=args.dk_disc,
         show_dk_disc_line=not args.no_dk_disc_line)
