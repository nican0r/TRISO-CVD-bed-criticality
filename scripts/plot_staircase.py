"""Visualize the inscribed staircase (cylindrical slab) approximation of the cone.

Usage:
    python scripts/plot_staircase.py [--n_slabs N] [--output PATH]
    python scripts/plot_staircase.py --compare 4 10 20 --output staircase.png

Draws a 2-D r-z cross-section showing the true cone frustum and inscribed slabs.
Self-contained: reads params.yaml directly, no openmc dependency.
"""

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

REPO_ROOT = Path(__file__).parent.parent
PARAMS_FILE = REPO_ROOT / "params.yaml"


def load_dimensions():
    with open(PARAMS_FILE) as f:
        raw = yaml.safe_load(f)
    dim_mm = raw["dimensions"]
    # Convert mm → cm
    return {
        "r_throat": dim_mm["nozzle"]["throat_diameter"] / 2.0 / 10.0,
        "r_retort":  dim_mm["retort"]["id"] / 2.0 / 10.0,
        "z_cone_top": dim_mm["cone"]["vertical_drop"] / 10.0,
        "half_angle_deg": dim_mm["cone"]["included_angle_deg"] / 2.0,
    }


def frustum_volume(r_top, r_bot, half_angle_deg):
    tan_theta = math.tan(math.radians(half_angle_deg))
    return (math.pi / 3.0) * (r_top**3 - r_bot**3) / tan_theta


def staircase_slabs(dim, n_slabs):
    r_throat    = dim["r_throat"]
    z_cone_top  = dim["z_cone_top"]
    tan_theta   = math.tan(math.radians(dim["half_angle_deg"]))
    z_apex      = -r_throat / tan_theta
    dz          = z_cone_top / n_slabs

    slabs = []
    for i in range(n_slabs):
        z_bot = i * dz
        z_top = (i + 1) * dz
        r = tan_theta * (z_bot - z_apex)   # inscribed: radius at slab bottom
        slabs.append({"z_bot": z_bot, "z_top": z_top, "r": r, "vol": math.pi * r**2 * dz})

    total_vol  = sum(s["vol"] for s in slabs)
    v_frustum  = frustum_volume(dim["r_retort"], r_throat, dim["half_angle_deg"])
    vol_err    = (v_frustum - total_vol) / v_frustum
    return slabs, vol_err


def draw_panel(ax, dim, n_slabs):
    r_throat   = dim["r_throat"]
    r_retort   = dim["r_retort"]
    z_cone_top = dim["z_cone_top"]
    tan_theta  = math.tan(math.radians(dim["half_angle_deg"]))

    slabs, vol_err = staircase_slabs(dim, n_slabs)

    # True cone interior (shaded)
    z_fill = np.linspace(0, z_cone_top, 300)
    r_fill = tan_theta * (z_fill + r_throat / tan_theta)
    ax.fill_betweenx(z_fill, -r_fill, r_fill, color="#d0d0d0", alpha=0.5, zorder=0,
                     label="True cone interior")

    # True cone walls
    for sign in (+1, -1):
        ax.plot([sign * r_throat, sign * r_retort], [0, z_cone_top],
                "k-", lw=2, zorder=4)

    # Staircase slabs
    colors = ["#4878CF", "#6ACC65"]
    for k, slab in enumerate(slabs):
        r, zb, zt = slab["r"], slab["z_bot"], slab["z_top"]
        c = colors[k % 2]
        ax.fill_betweenx([zb, zt], [-r, -r], [r, r],
                         color=c, alpha=0.6, linewidth=0, zorder=2)
        # Outline top and sides of each slab
        ax.plot([-r, -r], [zb, zt], color="navy", lw=0.7, zorder=3)
        ax.plot([ r,  r], [zb, zt], color="navy", lw=0.7, zorder=3)
        ax.plot([-r,  r], [zt, zt], color="navy", lw=0.7, zorder=3)
        ax.plot([-r,  r], [zb, zb], color="navy", lw=0.4, ls="--", zorder=3, alpha=0.6)

    # Reference lines
    ax.axhline(0,          color="black", lw=1.0)
    ax.axhline(z_cone_top, color="black", lw=0.6, ls="--", alpha=0.5)
    ax.axvline(0,          color="black", lw=0.5, ls=":", alpha=0.4)

    # Dimension arrows
    arrow_kw = dict(arrowstyle="<->", lw=1.2)
    # r_throat: annotate with arrow pointing down-left from just outside the throat
    ax.annotate(
        f"r_throat = {r_throat:.2f} cm",
        xy=(r_throat, 0.0),
        xytext=(r_retort * 0.55, z_cone_top * 0.12),
        fontsize=7, color="firebrick", ha="left", va="center",
        arrowprops=dict(arrowstyle="->", color="firebrick", lw=1.0),
    )

    ax.annotate("", xy=(r_retort, z_cone_top + 0.35), xytext=(0, z_cone_top + 0.35),
                arrowprops=dict(**arrow_kw, color="darkgreen"))
    ax.text(r_retort / 2, z_cone_top + 0.52,
            f"r_retort = {r_retort:.2f} cm", ha="center", va="bottom",
            fontsize=7, color="darkgreen")

    ax.set_xlim(-r_retort * 1.25, r_retort * 1.25)
    ax.set_ylim(0, z_cone_top * 1.18)
    ax.set_xlabel("Radius (cm)", fontsize=10)
    ax.set_title(f"n = {n_slabs} slabs\nVol. error: {vol_err*100:.1f}%", fontsize=10)
    ax.set_aspect("equal")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n_slabs", type=int, default=10)
    parser.add_argument("--output",  default=None)
    parser.add_argument("--compare", nargs="+", type=int, default=None,
                        help="Side-by-side comparison of multiple slab counts")
    args = parser.parse_args()

    dim    = load_dimensions()
    counts = args.compare if args.compare else [args.n_slabs]

    fig, axes = plt.subplots(1, len(counts),
                             figsize=(4.5 * len(counts), 6),
                             sharey=True)
    if len(counts) == 1:
        axes = [axes]

    for ax, n in zip(axes, counts):
        draw_panel(ax, dim, n)

    axes[0].set_ylabel("z (cm)", fontsize=10)

    fig.suptitle(
        "Inscribed staircase approximation of cone frustum\n"
        f"(60° included angle, z = 0 → {dim['z_cone_top']:.1f} cm, units: cm)",
        fontsize=11, y=1.01,
    )
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"Saved → {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
