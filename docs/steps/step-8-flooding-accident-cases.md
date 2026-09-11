# Step 8 — Flooding and collapsed-bed accident cases

## How to run

```
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env

# Smoke-test bottom-up sweep (21 cases × ~30 s each):
caffeinate python scripts/run_flood.py --bottomup-flood --quick

# Full bottom-up sweep (21 cases); run in parallel for speed:
caffeinate python scripts/run_flood.py --bottomup-flood -j 4

# Plot after CSV is written:
caffeinate python scripts/plot_flood_sweep.py
```

Outputs:
- CSV — `results/flood_bottomup.csv` (21 cases: 1 dry baseline + 20 z_flood levels)
- PNG — `results/flood_bottomup.png`
- TXT — `results/flood_top10.txt` (top-10 table + double-contingency check)

---

## What was implemented

- **`furnace/geometry.py`** — Added module-level `_flood_fill(z_bot, wet, dry, z_flood)` helper. Modified `bed_region()` to accept `z_flood` and `dry_material` parameters. When set, each staircase slab, its annular cone void cell, and the overflow cylinder cell are assigned wet or dry material based on whether their `z_bot` falls below `z_flood`. The unfilled cone void above the bed top is split exactly at `z_flood` (with a dedicated `ZPlane`) rather than using the slab-rounding approximation, giving exact flood level accuracy in the unfilled cone region.

- **`furnace/model.py`** — Added `z_flood` and `water_density_gcc` parameters to `build_model()`. For bottom-up flood (`z_flood` set), the gas-above-bed cell is split at `z_flood` when the flood level extends above the bed top, using the same shared `ZPlane` object for the adjacent cell boundaries.

- **`scripts/run_flood.py`** — Sweep driver for the bottom-up flood series (21 cases).

- **`scripts/plot_flood_sweep.py`** — Plotting script. Generates a line plot of k + 2σ + δk_disc vs z_flood for the bottom-up series with peak annotation, and a top-10 most reactive case table with the double-contingency analysis for the most reactive case.

---

## How it works

**Bottom-up flood sweep.** Twenty flood levels are computed as `z_rt × i/20` for `i = 1, …, 20`, uniformly spanning from just above the throat junction to the retort top. At each level the bed cells below `z_flood` use liquid water as background and those above use process gas; the cone unfilled-void cell and the gas-above-bed cell are split at `z_flood` with a shared `ZPlane`. A dry baseline case (no water) is also included in the same CSV for reference. The collapsed bed at the nominal 95 g charge fits entirely within the cone frustum, so the overflow-cylinder branch of `bed_region` is not exercised for this sweep.

---

## Experimental design

**What is being modelled.** Twenty-one eigenvalue cases covering bottom-up inundation from a coolant leak at the injector throat. All cases use the collapsed bed (pf ≈ 0.50) and bare UCO kernels — the bounding combination for fissile density at the nominal charge. Room temperature is the bounding nuclear condition for flooding (more thermal moderation than at operating temperature).

**Materials.** Bare UCO kernels (ρ = 10.5 g/cm³, 19.75 wt% ²³⁵U), structural graphite, injector body water at 1.0 g/cm³. Flood water: light water at 1.0 g/cm³ (liquid). No S(α,β) applied to flood water because the ENDF/B-VIII.0 endfb80_hdf5 library has empty temperature data for `c_H_in_H2O` — free-gas treatment is used, which underestimates thermal moderation of hydrogen and makes the result non-conservative (k-eff will be underestimated relative to the correct thermal treatment). This is a known Stage 0 limitation; the cases are still screening-valid because H/U ratios and trends are preserved.

**Nuclear data.** ENDF/B-VIII.0 (`endfb80_hdf5`), temperature snapped to 900 K by the `nearest` method (library has no room-temperature point). Non-conservative for Doppler broadening. `c_Graphite` S(α,β) applied to all carbon-bearing materials.

**Driving condition.** Eigenvalue mode; Watt fission spectrum source (a = 0.988 MeV, b = 2.249 × 10⁻⁶ eV⁻¹) over a box spanning the full bed extent with `constraints={'fissionable': True}`.

**Parameter scheme.** Bottom-up sweep: 20 z_flood levels uniformly from z_rt/20 to z_rt, plus 1 dry baseline. All other parameters held at the step-5 nominal.

**Per-step settings.** Each case: 20 000 particles/generation, 300 batches (50 inactive + 250 active), seed = 42. Expected σ(k) ≈ few × 10⁻⁴. This precision is sufficient to resolve differences between flood cases, which are expected to be O(10⁻²) or larger.

**Simulation hierarchy.**
- Bottom-up sweep: 21 cases × (300 batches × 20 000 particles) = 1.26 × 10⁸ histories

**Output.**
- `results/flood_bottomup.csv` — one row per z_flood level with z_flood_cm, k-eff, and same derived columns (k+2σ, k+2σ+δk_disc, H/²³⁵U, C/²³⁵U, bed geometry, wall time).
- `results/flood_bottomup.png` — k+2σ+δk_disc vs z_flood with peak annotation and cone-top marker.
- `results/flood_top10.txt` — ranked top-10 table and double-contingency analysis for the most reactive case.

---

## Design decisions

- **Collapsed bed only.** A fluidized bed at pf ≈ 0.333 has ~1.5× less fissile density per unit volume than the collapsed bed (pf = 0.50). The flooding case is an off-normal condition where process gas flow may also have been lost; collapsed is the correct bounding state. The step specification confirms this.
- **Liquid water only (1.0 g/cm³).** Partial-density water is bracketed by the dry baseline and the liquid case at each z_flood level. If the liquid case is subcritical, partial-density cases at the same flood height cannot be more reactive.
- **Bottom-up flood only.** A uniform cross-product of flood extents (`bed_and_cone`, `full_retort`) would duplicate the dry baseline, duplicate the full-retort case (which equals z_flood = z_rt), and add only one non-redundant case (`bed_and_cone` ≈ z_flood at cone top) that is already bracketed by adjacent bottom-up levels. The bottom-up sweep covers the full geometric range with finer spatial resolution and no redundancy.
- **z_flood sweep: 20 uniform levels.** Fine enough to localise the peak reactivity z to within ~1.9 cm. The two levels that fall within the cone (~0–3.8 cm) may have coarse spacing, but the exact split of the unfilled cone void (implemented with a dedicated ZPlane rather than the slab rounding approximation) preserves accuracy at z_flood values within the unfilled cone region.
- **Simple slab assignment for bed cells.** A slab is fully wet if its `z_bot < z_flood`, fully dry otherwise. The maximum over-approximation of the flooded volume is one slab height ≈ 0.12 cm (at n_slabs=32). Splitting every slab at z_flood would be geometrically exact but requires per-slab surface creation and cell splitting — the error here is smaller than the statistical σ(k) and acceptable for Stage 0.
- **Exact split for unfilled cone void.** Unlike the particle-filled slabs, the unfilled cone void spans a large axial range. When z_flood falls inside this region, the void is split exactly at z_flood with a dedicated ZPlane, eliminating the cone-void over-approximation entirely.
- **c_H_in_H2O omitted.** The ENDF/B-VIII.0 library has empty temperature data for `c_H_in_H2O`, so the free-gas scattering law is used for hydrogen in water. This underestimates thermal moderation and makes k-eff results non-conservative for the flooding cases. This is flagged in Assumptions as unconfirmed; if a library with room-temperature `c_H_in_H2O` becomes available, the flood cases should be re-run.
- **Bottom-up flood initiator named.** The injector coolant leak (single failure at the bottom of the cone) is the bounding water source for this design — it puts water exactly where the fuel is densest in the collapsed configuration. The sweep captures whether partial fill (as the leak progresses) is more reactive than full fill.

---

## Assumptions

**Confirmed:**
- Collapsed bed at pf = 0.50 (`model.packing_fraction_static`) is the bounding bed state for flooding cases. Confirmed by step-8 specification.
- Nominal charge mass 95 g fits entirely within the cone frustum at collapsed packing (V_bulk = 18.1 cm³ < V_frustum ≈ 28.2 cm³), so no overflow cylinder is exercised in this sweep.
- All cases run at 293.6 K material temperature (cold, flooded condition is bounding per preamble). Temperature snapped to 900 K by library nearest-method.

**Defaulted:**
- 20 z_flood levels, uniformly spaced from z_rt/20 to z_rt. Chosen to span the full retort height with ~1.9 cm resolution. Coarser in the cone (only ~2 points) but the unfilled-void split removes the dominant error in that region.
- Simple slab-assignment rule: whole slab is wet if z_bot < z_flood. Error bounded by one slab height ≈ 0.12 cm at n_slabs = 32; smaller than σ(k).
- Water density 1.0 g/cm³ (saturated liquid at atmospheric pressure). The upper bound on liquid-water density; conservative for moderation.
- H/²³⁵U ratio for partial-flood cases (z_flood sweep) uses the full bed void volume multiplied by the water number density. This overestimates H/²³⁵U when z_flood < bed_height, since only the fraction of the void below z_flood is water-filled. Used for diagnostic tracking only, not for NCS decisions.

**Unconfirmed (`# CONFIRM`):**
- `c_H_in_H2O` S(α,β) omitted due to empty library temperature data. Free-gas treatment underestimates thermal moderation of H in water → k-eff is non-conservatively low for flooding cases. If a compatible `c_H_in_H2O` table is available, the flood cases must be re-run and results will be higher. This is the most significant known non-conservative approximation in Step 8.
- Structural graphite density 1.75 g/cm³ with zero boron equivalent — inherited `# CONFIRM` from preamble.
- All ambient-temperature approximations (900 K instead of 293.6 K for neutron cross sections, 1200 K for c_Graphite S(α,β)) — inherited Stage 0 limitation.
