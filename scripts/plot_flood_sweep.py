"""Step 8 flood sweep plots.

Reads:
  results/flood_bottomup.csv — bottom-up z_flood sweep

Writes:
  results/flood_bottomup.png  — k+2σ+δk vs z_flood for bottom-up series
  results/flood_top10.txt     — top-10 most reactive cases

Run after the sweep completes:
    caffeinate python scripts/plot_flood_sweep.py
"""
from __future__ import annotations

import csv
import sys
import textwrap
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

_BU_CSV    = _REPO_ROOT / 'results' / 'flood_bottomup.csv'
_BU_PNG    = _REPO_ROOT / 'results' / 'flood_bottomup.png'
_TOP10_TXT = _REPO_ROOT / 'results' / 'flood_top10.txt'

_SUBCRIT = 0.95
_FOOTER  = 'UNVALIDATED — SCREENING ONLY. Not benchmarked against ICSBEP.'


def _load(csv_path: Path) -> list[dict]:
    if not csv_path.exists():
        return []
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for col in ('k_eff', 'sigma', 'k_plus_2sigma', 'k_plus_2sigma_plus_dk_disc',
                    'bed_height_cm', 'h_per_u235', 'u235_mass_g'):
            try:
                r[col] = float(r[col])
            except (ValueError, KeyError):
                r[col] = float('nan')
        for col in ('z_flood_cm', 'water_density_gcc'):
            try:
                r[col] = float(r[col]) if r.get(col, '') else None
            except ValueError:
                r[col] = None
    return rows


# ---------------------------------------------------------------------------
# Bottom-up flood line plot
# ---------------------------------------------------------------------------

def plot_bottomup(rows: list[dict]) -> None:
    if not rows:
        print(f'  No data in {_BU_CSV} — skipping bottom-up plot')
        return

    flood_rows = sorted(
        [r for r in rows if r.get('z_flood_cm') is not None],
        key=lambda r: r['z_flood_cm'],
    )
    dry_rows = [r for r in rows if r.get('z_flood_cm') is None]

    if not flood_rows:
        print('  No z_flood rows found — skipping bottom-up plot')
        return

    z    = [r['z_flood_cm']              for r in flood_rows]
    k    = [r['k_eff']                   for r in flood_rows]
    sig  = [r['sigma']                   for r in flood_rows]
    k2sd = [r['k_plus_2sigma_plus_dk_disc'] for r in flood_rows]

    dry_k2sd = dry_rows[0]['k_plus_2sigma_plus_dk_disc'] if dry_rows else None

    try:
        from furnace.params import load_params
        params = load_params()
        dim = params['dimensions']
        z_cone = dim['cone']['vertical_drop_cm']
        z_rt   = z_cone + dim['retort']['height_cm']
    except Exception:
        z_cone, z_rt = None, None

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar(z, k, yerr=[2 * s for s in sig], fmt='o-', capsize=4,
                color='steelblue', label='k-eff ± 2σ')
    ax.plot(z, k2sd, 's--', color='firebrick',
            label='k + 2σ + δk_disc (conservative)')
    ax.axhline(_SUBCRIT, color='black', lw=1.2, ls=':', label=f'{_SUBCRIT} subcritical limit')

    if dry_k2sd is not None:
        ax.axhline(dry_k2sd, color='gray', lw=1.0, ls='--',
                   label=f'dry baseline k+2σ+δk = {dry_k2sd:.4f}')

    if z_cone is not None:
        ax.axvline(z_cone, color='darkorange', lw=1.0, ls=':', alpha=0.7,
                   label=f'cone top z={z_cone:.1f} cm')

    if k2sd:
        i_peak = int(np.nanargmax(k2sd))
        ax.annotate(f'peak\nk+2σ+δk={k2sd[i_peak]:.4f}\nz={z[i_peak]:.2f} cm',
                    xy=(z[i_peak], k2sd[i_peak]),
                    xytext=(10, 12), textcoords='offset points',
                    fontsize=8, arrowprops=dict(arrowstyle='->', color='firebrick'))

    ax.set_xlabel('Bottom-up flood level z_flood (cm above cone base / throat junction)')
    ax.set_ylabel('k-eff')
    ax.set_title('Bottom-up flood sweep — collapsed bed, bare kernels, liquid water\n'
                 '(UNVALIDATED — SCREENING ONLY)', fontsize=10)
    ax.grid(True, lw=0.3)
    ax.legend(fontsize=8, loc='best')
    fig.text(0.5, 0.01, _FOOTER, ha='center', fontsize=8, style='italic', color='dimgray')
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(str(_BU_PNG), dpi=150)
    plt.close(fig)
    print(f'Bottom-up plot → {_BU_PNG}')

    if k2sd:
        i_peak = int(np.nanargmax(k2sd))
        print(f'\nBottom-up flood peak: z_flood={z[i_peak]:.3f} cm, '
              f'k-eff={k[i_peak]:.4f}±{sig[i_peak]:.4f}, '
              f'k+2σ+δk={k2sd[i_peak]:.4f}')
        if dry_k2sd is not None:
            delta = k2sd[i_peak] - dry_k2sd
            print(f'  Δ(k+2σ+δk) vs dry baseline: {delta:+.4f}')


# ---------------------------------------------------------------------------
# Top-10 table + double-contingency check
# ---------------------------------------------------------------------------

def top10_table(rows: list[dict]) -> None:
    if not rows:
        print('  No data for top-10 table')
        return

    ranked = sorted(
        [r for r in rows if not np.isnan(r['k_plus_2sigma_plus_dk_disc'])],
        key=lambda r: r['k_plus_2sigma_plus_dk_disc'],
        reverse=True,
    )

    top = ranked[:10]
    best = top[0] if top else None

    lines = []
    lines.append('Top-10 most reactive flood cases (k + 2σ + δk_disc)')
    lines.append('=' * 100)
    hdr = (f"{'Rank':>4}  {'Tag':<40}  {'z_flood_cm':>10}  "
           f"{'k_eff':>7}  {'±σ':>6}  {'k+2σ+δk':>9}")
    lines.append(hdr)
    lines.append('-' * 100)
    for i, r in enumerate(top, 1):
        zf   = f"{r['z_flood_cm']:.3f}" if r['z_flood_cm'] is not None else '— (dry)'
        flag = ' ⚠' if r['k_plus_2sigma_plus_dk_disc'] >= _SUBCRIT else ''
        lines.append(
            f"{i:>4}  {r['tag']:<40}  "
            f"{zf:>10}  {r['k_eff']:>7.4f}  {r['sigma']:>6.4f}  "
            f"{r['k_plus_2sigma_plus_dk_disc']:>9.4f}{flag}"
        )
    lines.append('')

    if best:
        lines.append('Most reactive case — full specification:')
        lines.append('-' * 60)
        lines.append(f"  Tag              : {best['tag']}")
        lines.append(f"  State            : {best.get('state', '—')}")
        lines.append(f"  Stage            : {best.get('stage', '—')}")
        lines.append(f"  z_flood_cm       : {best['z_flood_cm'] if best['z_flood_cm'] is not None else '— (dry)'}")
        lines.append(f"  Water density    : {best.get('water_density_gcc', '—')} g/cm³")
        lines.append(f"  k-eff            : {best['k_eff']:.4f} ± {best['sigma']:.4f}")
        lines.append(f"  k + 2σ + δk_disc : {best['k_plus_2sigma_plus_dk_disc']:.4f}")
        lines.append(f"  U-235 mass       : {best['u235_mass_g']:.2f} g")
        lines.append(f"  Bed height       : {best['bed_height_cm']:.3f} cm")
        lines.append(f"  H/²³⁵U ratio     : {best['h_per_u235']:.1f}")
        lines.append('')
        lines.append('Double-contingency check for most reactive case:')
        lines.append('-' * 60)
        lines.append(textwrap.dedent("""\
          For this case to reach its k-eff, the following independent parameters
          must be simultaneously off-normal:

            1. Flood water present
               Nominal: dry CVD atmosphere (no water in retort interior)
               Required: liquid water at 1.0 g/cm³ in the retort
               Initiator: injector coolant leak (bounding source per step-8 spec),
                          fire sprinkler discharge, hydrostatic test water retention,
                          or steam ingress on cooldown.
               Likelihood: single failure (loss of coolant boundary integrity).

            2. Bed in collapsed state
               Nominal: fluidized bed during CVD operation (pf ≈ 0.333)
               Required: collapsed/settled bed (pf ≈ 0.50, ~1.5× denser)
               Initiator: loss of process gas flow (pump failure, valve closure,
                          power outage) — simultaneous with or prior to flooding.
               Likelihood: independent failure of process gas system.

            3. Charge mass at accumulated inventory level
               Nominal: 95 g per loading cycle
               Required: whatever mass produced this run's k-eff
               Initiator: improper inventory management or unintended accumulation
                          between cycles.
               Likelihood: administrative failure (procedure violation).

          Reaching the reported k-eff requires at minimum failures (1) and (2)
          simultaneously. Adding (3) is only relevant if the charge is above nominal.
          Two independent hardware failures are required; this is a double-contingency
          scenario. If k-eff < 0.95 under double contingency, the double-contingency
          argument supports subcriticality.

          NOTE: Results are unvalidated Stage-0 screening calculations.
        """))

    output = '\n'.join(lines)
    print(output)
    _TOP10_TXT.parent.mkdir(parents=True, exist_ok=True)
    _TOP10_TXT.write_text(output)
    print(f'\nTop-10 table → {_TOP10_TXT}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print('Loading CSV...')
    bu_rows = _load(_BU_CSV)
    print(f'  flood_bottomup.csv: {len(bu_rows)} rows')

    print('\nGenerating bottom-up flood plot...')
    plot_bottomup(bu_rows)

    print('\nGenerating top-10 table...')
    top10_table(bu_rows)


if __name__ == '__main__':
    main()
