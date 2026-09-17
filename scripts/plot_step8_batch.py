"""Plot step8_flood results from AWS Batch summary.csv.

Reads:
  results/step8_flood/<submission>/summary.csv

Writes:
  results/step8_flood/flood_sweep.png
  results/step8_flood/flood_top10.txt

Usage:
    python scripts/plot_step8_batch.py [--submission <id>]
"""
from __future__ import annotations

import argparse
import csv
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
_RESULTS   = _REPO_ROOT / "results" / "step8_flood"
_SUBCRIT   = 0.95
_FOOTER    = (
    "UNVALIDATED — Stage-0 screening calculation. Cross-section model not benchmarked against ICSBEP critical experiments.\n"
    "c_H_in_H₂O S(α,β) omitted (non-conservative for flooding). Results bound geometry and reactivity trends only — not suitable for licensing."
)

# Geometry constants from furnace.params
_Z_CONE_TOP = 3.8    # cm — cone top above throat junction
_Z_RT_TOP   = 37.8   # cm — retort top (3.8 + 34.0)


def _latest_submission(results_dir: Path) -> Path:
    subs = sorted(d for d in results_dir.iterdir() if d.is_dir())
    if not subs:
        raise SystemExit(f"No submission directories found in {results_dir}")
    return subs[-1]


def _load(csv_path: Path) -> list[dict]:
    with csv_path.open() as f:
        raw = list(csv.DictReader(f))
    rows = []
    for r in raw:
        row = {
            "state":              r["row.state"],
            "water_density_gcc":  float(r["row.water_density_gcc"]),
            "z_flood_cm":         float(r["row.z_flood"]),
            "keff":               float(r["keff"]),
            "keff_std":           float(r["keff_std"]),
            "keff_plus_2sigma":   float(r["keff_plus_2sigma"]),
            "u235_mass_g":        float(r["u235_mass_g"]),
            "pf_achieved":        float(r["pf_achieved"]),
            "bed_height_cm":      float(r["bed_height_cm"]),
            "wall_seconds":       float(r["wall_seconds"]),
        }
        rows.append(row)
    return rows


def _series_label(state: str, density: float) -> str:
    phase = "liquid" if density >= 0.5 else "vapor"
    return f"{state}, {phase} ({density} g/cm³)"


def _series_style(state: str, density: float) -> dict:
    color = {
        ("collapsed",  1.0):   "steelblue",
        ("collapsed",  0.001): "cornflowerblue",
        ("fluidized",  1.0):   "firebrick",
        ("fluidized",  0.001): "salmon",
    }.get((state, density), "gray")
    ls = "-" if density >= 0.5 else "--"
    return {"color": color, "ls": ls}


def plot_sweep(rows: list[dict], out_png: Path) -> None:
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        key = (r["state"], r["water_density_gcc"])
        groups.setdefault(key, []).append(r)
    for key in groups:
        groups[key].sort(key=lambda r: r["z_flood_cm"])

    fig, ax = plt.subplots(figsize=(11, 6))

    overall_peak_k2s = -np.inf
    overall_peak_row = None

    for (state, density), series in sorted(groups.items()):
        z    = [r["z_flood_cm"]       for r in series]
        k    = [r["keff"]             for r in series]
        sig  = [r["keff_std"]         for r in series]
        k2s  = [r["keff_plus_2sigma"] for r in series]
        sty  = _series_style(state, density)
        lbl  = _series_label(state, density)

        ax.errorbar(z, k, yerr=[2 * s for s in sig],
                    fmt="o", capsize=3, ms=4,
                    color=sty["color"], alpha=0.7, zorder=2)
        ax.plot(z, k2s, marker="s", ms=4,
                color=sty["color"], ls=sty["ls"], lw=1.4,
                label=lbl + "  (k+2σ)", zorder=3)

        i_peak = int(np.nanargmax(k2s))
        if k2s[i_peak] > overall_peak_k2s:
            overall_peak_k2s = k2s[i_peak]
            overall_peak_row = series[i_peak]
            overall_peak_label = lbl

    # Geometry markers
    ax.axvline(_Z_CONE_TOP, color="darkorange", lw=1.0, ls=":", alpha=0.8,
               label=f"cone top z={_Z_CONE_TOP:.1f} cm")
    ax.axvline(_Z_RT_TOP, color="purple", lw=1.0, ls=":", alpha=0.6,
               label=f"retort top z={_Z_RT_TOP:.1f} cm")

    # Auto-scale y-axis to the data range with 20% padding, ignoring the 0.95 limit line
    all_k2s = [v for series in groups.values() for v in [r["keff_plus_2sigma"] for r in series]]
    all_k   = [v for series in groups.values() for v in [r["keff"] for r in series]]
    all_sig = [v for series in groups.values() for v in [r["keff_std"] for r in series]]
    y_min = min(k - 2 * s for k, s in zip(all_k, all_sig))
    y_max = max(all_k2s)
    y_range = y_max - y_min
    ax.set_ylim(max(0.0, y_min - y_range * 0.45), y_max + y_range * 0.20)

    # Subcritical limit (not labelled — off-screen at k=0.95)
    ax.axhline(_SUBCRIT, color="black", lw=1.2, ls=":")


    ax.set_xlabel("Bottom-up flood level z_flood (cm above cone throat)", fontsize=11)
    ax.set_ylabel("k-eff", fontsize=11)
    ax.set_title(
        "Step 8 — Bottom-up flood sweep  |  bare UCO kernels, 293.6 K, 95 g charge\n"
        "2 bed states × 2 water densities × 20 z_flood levels (80 cases)",
        fontsize=10,
    )
    ax.grid(True, lw=0.3, alpha=0.5)
    ax.legend(fontsize=8, loc="lower right", ncol=1)
    fig.text(0.5, 0.01, _FOOTER, ha="center", fontsize=7, style="italic", color="dimgray")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_png), dpi=150)
    plt.close(fig)
    print(f"Plot → {out_png}")

    if overall_peak_row is not None:
        r = overall_peak_row
        print(f"\nOverall peak (k+2σ={r['keff_plus_2sigma']:.4f}):")
        print(f"  state={r['state']}, density={r['water_density_gcc']} g/cm³, "
              f"z_flood={r['z_flood_cm']:.2f} cm")
        print(f"  keff={r['keff']:.4f} ± {r['keff_std']:.4f}")


def top10_table(rows: list[dict], out_txt: Path) -> None:
    ranked = sorted(rows, key=lambda r: r["keff_plus_2sigma"], reverse=True)
    top = ranked[:10]
    best = top[0] if top else None

    lines = []
    lines.append("Top-10 most reactive flood cases (k + 2σ, bare kernels)")
    lines.append("=" * 105)
    hdr = (f"{'Rank':>4}  {'State':<10}  {'Density g/cc':>12}  {'z_flood cm':>10}  "
           f"{'k-eff':>7}  {'±σ':>6}  {'k+2σ':>8}")
    lines.append(hdr)
    lines.append("-" * 105)
    for i, r in enumerate(top, 1):
        flag = " ⚠" if r["keff_plus_2sigma"] >= _SUBCRIT else ""
        lines.append(
            f"{i:>4}  {r['state']:<10}  {r['water_density_gcc']:>12.3f}  "
            f"{r['z_flood_cm']:>10.2f}  {r['keff']:>7.4f}  {r['keff_std']:>6.4f}  "
            f"{r['keff_plus_2sigma']:>8.4f}{flag}"
        )
    lines.append("")

    if best:
        lines.append("Most reactive case — full specification:")
        lines.append("-" * 60)
        lines.append(f"  State            : {best['state']}")
        lines.append(f"  Water density    : {best['water_density_gcc']} g/cm³")
        lines.append(f"  z_flood_cm       : {best['z_flood_cm']:.2f} cm")
        lines.append(f"  k-eff            : {best['keff']:.4f} ± {best['keff_std']:.4f}")
        lines.append(f"  k + 2σ           : {best['keff_plus_2sigma']:.4f}")
        lines.append(f"  U-235 mass       : {best['u235_mass_g']:.2f} g")
        lines.append(f"  pf achieved      : {best['pf_achieved']:.4f}")
        lines.append(f"  Bed height       : {best['bed_height_cm']:.3f} cm")
        lines.append(f"  Wall time        : {best['wall_seconds']/60:.1f} min")
        lines.append("")
        lines.append("Double-contingency check for most reactive case:")
        lines.append("-" * 60)
        lines.append(textwrap.dedent("""\
          For this case to reach its k-eff, the following independent parameters
          must be simultaneously off-normal:

            1. Flood water present
               Nominal: dry CVD atmosphere (no water in retort interior)
               Required: liquid water in the retort
               Initiator: injector coolant leak, sprinkler discharge, steam ingress.
               Likelihood: single failure.

            2. Bed in collapsed state
               Nominal: fluidized during CVD operation (pf ≈ 0.333)
               Required: collapsed/settled (pf ≈ 0.50, ~1.5× denser)
               Initiator: loss of process gas flow (pump failure, valve closure).
               Likelihood: independent failure of process gas system.

            3. Charge mass at accumulated inventory level
               Nominal: 95 g per loading cycle
               Initiator: improper inventory management.
               Likelihood: administrative failure.

          Reaching the reported k-eff requires at minimum failures (1) and (2)
          simultaneously. Two independent hardware failures are required — this is
          a double-contingency scenario.

          NOTE: Results are unvalidated Stage-0 screening calculations.
          NOTE: c_H_in_H2O S(α,β) omitted (empty library data); free-gas treatment
                underestimates thermal moderation → k-eff is non-conservatively low.
        """))

    out_text = "\n".join(lines)
    print(out_text)
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    out_txt.write_text(out_text)
    print(f"\nTop-10 table → {out_txt}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--submission", default=None,
                    help="Submission ID directory name; defaults to newest.")
    args = ap.parse_args()

    if args.submission:
        sub_dir = _RESULTS / args.submission
    else:
        sub_dir = _latest_submission(_RESULTS)

    csv_path = sub_dir / "summary.csv"
    if not csv_path.exists():
        raise SystemExit(f"summary.csv not found: {csv_path}")

    print(f"Loading {csv_path} ...")
    rows = _load(csv_path)
    print(f"  {len(rows)} rows loaded")

    out_png = _RESULTS / "flood_sweep.png"
    out_txt = _RESULTS / "flood_top10.txt"

    print("\nGenerating plot...")
    plot_sweep(rows, out_png)

    print("\nGenerating top-10 table...")
    top10_table(rows, out_txt)


if __name__ == "__main__":
    main()
