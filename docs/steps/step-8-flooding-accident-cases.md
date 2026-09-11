# Step 8 — Flooding accident cases: bottom-up flood sweep

## How to run

```
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env

# Smoke-test bottom-up sweep (~30 s/case):
caffeinate python scripts/run_flood.py --bottomup-flood --quick -j 4

# Full bottom-up sweep (82 cases); run in parallel for speed:
caffeinate python scripts/run_flood.py --bottomup-flood -j 4

# Plot after CSV is written:
caffeinate python scripts/plot_flood_sweep.py
```

Outputs:
- CSV — `results/flood_bottomup.csv` (80 cases: 2 states × 20 z_flood levels × 2 densities)
- PNG — `results/flood_bottomup.png`
- TXT — `results/flood_top10.txt` (top-10 table + double-contingency check)

---

## What was implemented

- **`furnace/geometry.py`** — Added module-level `_flood_fill(z_bot, wet, dry, z_flood)` helper. Modified `bed_region()` to accept `z_flood` and `dry_material` parameters. When set, each staircase slab, its annular cone void cell, and the overflow cylinder cell are assigned wet or dry material based on whether their `z_bot` falls below `z_flood`. The unfilled cone void above the bed top is split exactly at `z_flood` (with a dedicated `ZPlane`) rather than using the slab-rounding approximation, giving exact flood level accuracy in the unfilled cone region.

- **`furnace/model.py`** — Added `flood_extent`, `z_flood`, and `water_density_gcc` parameters to `build_model()`. Material selection logic determines bed-fill and above-bed-fill materials based on the flood scenario. For bottom-up flood (`z_flood` set), the gas-above-bed cell is split at `z_flood` when the flood level extends above the bed top, using the same shared `ZPlane` object for adjacent cell boundaries.

- **`furnace/sweeps.py`** — Added three CSV columns (`flood_extent`, `z_flood_cm`, `water_density_gcc`) to `_CSV_FIELDNAMES`. Extended `run_case()` to extract `_flood_extent`, `_z_flood`, and `_water_density` from underscore-prefixed override keys and forward them to `build_model()`. Updated the `fill_mat` selection for H/²³⁵U ratio computation so that water is used when any flooding is present.

- **`scripts/run_flood.py`** — Bottom-up flood sweep driver. Sweeps 20 z_flood levels from throat to retort top for both bed states (collapsed, fluidized) and both water densities (liquid 1.0 g/cm³, vapor 0.001 g/cm³). 80 cases total; dry baselines omitted (already run in step 5). CSV → `results/flood_bottomup.csv`.

- **`scripts/plot_flood_sweep.py`** — Generates a k+2σ+δk_disc vs z_flood line plot for the bottom-up series with peak annotation, and a top-10 most reactive case table with double-contingency analysis.

---

## How it works

The bottom-up flood sweep models the injector coolant leak: water enters at the cone throat and fills upward. This is the bounding initiator for this design because the cone is where a collapsed bed concentrates its mass — partial fill from below puts moderator exactly where the fuel is densest. Fills from other directions (e.g. sprinklers from above) are less conservative.

Twenty flood levels are swept from just above the throat to the retort top for each combination of bed state and water density:

- **Bed states:** collapsed (pf = 0.50), fluidized (pf = 0.333)
- **Water densities:** liquid (1.0 g/cm³), vapor (0.001 g/cm³, near-atmospheric steam)
- **z_flood levels:** 20 uniform steps from z_rt/20 to z_rt
- **Dry baselines:** one per bed state for reference

At each z_flood, bed slab cells with `z_bot < z_flood` are assigned water as background; those above use process gas. The unfilled cone void is split exactly at z_flood with a dedicated `ZPlane`. All cases run at 293.6 K (room temperature, bounding for neutron moderation).

**Expect a reactivity peak at partial fill, not full flood.** Report the flood level at which k-eff peaks. If the peak is at partial fill, say so explicitly — draining a partial coolant leak partway does not necessarily improve the nuclear safety condition.

---

## Experimental design

**Scope.** 82 eigenvalue cases covering the bottom-up flooding scenario for two bed states and two water phases. All cases use bare UCO kernels. Room temperature is the bounding nuclear condition for flooding.

**Rationale for both bed states.** Collapsed (pf ≈ 0.50) has higher fissile density per unit bed volume. Fluidized (pf ≈ 0.333) has higher void fraction, giving approximately 2× the H/²³⁵U ratio when flooded:

| State | pf | Void fraction | H/²³⁵U ratio (relative) |
|---|---|---|---|
| Collapsed | 0.50 | 0.50 | 1.0× |
| Fluidized | 0.333 | 0.667 | ~2.0× |

Whichever state produces the higher k-eff under flooding is the bounding case; this is determined empirically by running both.

**Rationale for two water densities.** Liquid (1.0 g/cm³) represents saturated water from the coolant leak. Vapor (0.001 g/cm³) represents near-atmospheric steam — relevant if the leak flashes to steam in the hot retort environment or if steam ingress occurs on cooldown. Both are tested because the reactivity optimum may not be at maximum density; undermoderated systems can peak at a lower hydrogen density.

**Materials.** Bare UCO kernels (ρ = 10.5 g/cm³, 19.75 wt% ²³⁵U), structural graphite, injector body water. No S(α,β) applied to flood water (library limitation — see Assumptions).

**Nuclear data.** ENDF/B-VIII.0 (`endfb80_hdf5`), temperature snapped to 900 K by the `nearest` method. `c_Graphite` S(α,β) applied to all carbon-bearing materials.

**Parameter scheme.** 2 bed states × 20 z_flood levels × 2 water densities = 80 cases. Dry baselines omitted — already run in step 5.

**Per-case settings.** 20 000 particles/generation, 300 batches (50 inactive + 250 active), seed = 42. Expected σ(k) ≈ few × 10⁻⁴.

**Simulation totals.**
- 80 cases × (300 batches × 20 000 particles) ≈ 4.80 × 10⁸ histories

**Output.**
- `results/flood_bottomup.csv` — one row per case: state, z_flood_cm, water_density_gcc, k-eff, σ, k+2σ, k+2σ+δk_disc, H/²³⁵U, C/²³⁵U, bed geometry, wall time.
- `results/flood_bottomup.png` — k+2σ+δk_disc vs z_flood for both bed states and water densities with peak annotation and cone-top marker.
- `results/flood_top10.txt` — ranked top-10 table and double-contingency analysis for the most reactive case.

---

## Design decisions

- **Both collapsed and fluidized bed states.** Lower fissile density per volume is not always conservative when the system is undermoderated: the fluidized bed has ~2× the H/²³⁵U ratio when flooded (void fraction 0.667 vs 0.500), which can increase k-eff if the collapsed case sits below the moderation optimum. Both states are run and the higher k-eff governs.
- **Two water densities (liquid and vapor).** Reactivity as a function of moderator density is not monotone; the optimum H/²³⁵U ratio may be reached at partial density. Testing vapor explicitly ensures the steam-ingress and flash-vaporisation initiators are covered.
- **Bottom-up flood only.** The injector coolant leak fills from the cone throat upward, placing moderator exactly where fuel density is highest in the collapsed configuration. Uniform inundation from above or externally is not tested because it is bounded by the bottom-up scenario for the most reactive spatial configuration.
- **z_flood sweep: 20 uniform levels per state/density.** Fine enough to localise the peak reactivity z to within ~1.9 cm. The exact split of the unfilled cone void (implemented with a dedicated ZPlane) removes the dominant geometric error within the unfilled cone region.
- **Simple slab assignment for bed cells.** A slab is fully wet if its `z_bot < z_flood`, fully dry otherwise. Maximum over-approximation of the flooded volume is one slab height ≈ 0.12 cm (at n_slabs = 32); smaller than σ(k) and acceptable for Stage 0.
- **Exact split for unfilled cone void.** Unlike particle-filled slabs, the unfilled cone void spans a large axial range. When z_flood falls inside this region, the void is split exactly at z_flood with a dedicated ZPlane, eliminating the cone-void over-approximation entirely.
- **c_H_in_H2O omitted.** The ENDF/B-VIII.0 library has empty temperature data for `c_H_in_H2O`, so the free-gas scattering law is used for hydrogen in water. This underestimates thermal moderation and makes k-eff results non-conservative for the flooding cases. Results must be re-run if a compatible `c_H_in_H2O` table becomes available.

---

## Assumptions

**Confirmed:**
- Nominal charge mass 95 g fits entirely within the cone frustum at collapsed packing (V_bulk = 18.1 cm³ < V_frustum ≈ 28.2 cm³); no overflow cylinder is exercised in this sweep.
- All cases run at 293.6 K material temperature. Temperature snapped to 900 K by library nearest-method.

**Defaulted:**
- 20 z_flood levels per state/density, uniformly spaced from z_rt/20 to z_rt. Coarser in the cone (~2 points) but the unfilled-void split removes the dominant error in that region.
- Simple slab-assignment rule: whole slab is wet if z_bot < z_flood. Error bounded by one slab height ≈ 0.12 cm at n_slabs = 32.
- Vapor density 0.001 g/cm³ (near-atmospheric steam). Representative of steam at ~100 °C and 1 atm (ρ ≈ 0.0006 g/cm³); value rounded up slightly for conservatism.
- H/²³⁵U ratio for partial-flood cases uses the full bed void volume multiplied by water number density. Overestimates H/²³⁵U when z_flood < bed_height; used for diagnostic tracking only.

**Unconfirmed (`# CONFIRM`):**
- `c_H_in_H2O` S(α,β) omitted due to empty library temperature data. Free-gas treatment underestimates thermal moderation of H in water → k-eff is non-conservatively low for flooding cases. Most significant known non-conservative approximation in Step 8.
- Structural graphite density 1.75 g/cm³ with zero boron equivalent — inherited `# CONFIRM` from preamble.
- All ambient-temperature approximations (900 K instead of 293.6 K for neutron cross sections) — inherited Stage 0 limitation.
