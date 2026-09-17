---
group: Flooding accident cases
title: Bottom-up flood sweep — k-eff vs. flood level for collapsed and fluidized beds
---

## How to run

```bash
micromamba activate triso-env

# Smoke-test (~30 s/case):
caffeinate python scripts/run_flood.py --bottomup-flood --quick -j 4

# Full sweep (80 cases):
caffeinate python scripts/run_flood.py --bottomup-flood -j 4

# Plot:
caffeinate python scripts/plot_flood_sweep.py
```

Or submit to AWS Batch (see `docs/aws-batch.md`):
```bash
python scripts/batch/submit_sweep.py --sweep step8_flood --manifest manifests/step8_flood_bottomup.json --dry-run
python scripts/batch/submit_sweep.py --sweep step8_flood --manifest manifests/step8_flood_bottomup.json
python scripts/batch/pull_results.py --sweep step8_flood
```

Outputs:
- `results/flood_bottomup.csv` — one row per case
- `results/flood_bottomup.png` — k+2σ vs z_flood for all scenarios with peak annotation
- `results/flood_top10.txt` — top-10 most reactive cases

## What was implemented

- **`furnace/geometry.py`** — `exact_cone_bed` accepts `z_flood` and `dry_material`; the above-bed void cell is split at z_flood with a dedicated `ZPlane` for exact flood level boundary
- **`furnace/model.py`** — `build_model` accepts `z_flood` and `water_density_gcc`; material selection routes bed cells to wet or dry background based on flood level
- **`scripts/run_flood.py`** — sweeps 20 uniform z_flood levels from throat to retort top for 2 bed states × 2 water densities = 80 cases
- **`scripts/plot_flood_sweep.py`** — k+2σ line plot plus top-10 table

## How it works

The bottom-up flood models an injector coolant leak: water enters at the cone throat and fills upward. This is the bounding flood initiator because the collapsed bed concentrates its mass at the cone base — partial fill from below puts moderator exactly where fuel density is highest.

Sweep matrix: 2 bed states (collapsed pf=0.50, fluidized pf=0.333) × 20 z_flood levels × 2 water densities (liquid 1.0 g/cm³, vapor 0.001 g/cm³) = 80 cases. Dry baselines are the nominal-case runs.

The unfilled cone void above the bed is split exactly at z_flood with a dedicated `ZPlane`. All cases: 20 000 particles/gen, 300 batches (50 inactive + 250 active), seed 42, 293.6 K material temperature.

**Note:** `c_H_in_H2O` S(α,β) is omitted (empty library data); free-gas scattering is used for flood water. This underestimates thermal moderation — k-eff results are non-conservatively low for the flooding cases.

## Results

Canonical run: `step8_flood/1789561376_2144557fdce0`, git `030835a`, manifest `manifests/step8_flood_bottomup.json`.

**Peak k-eff by scenario (95 g nominal charge):**

| state | water density (g/cm³) | z_flood at peak (cm) | peak k-eff | peak k+2σ |
|---|---|---|---|---|
| collapsed | 1.0 (liquid) | 11.34 | 0.026771 | 0.026857 |
| fluidized | 1.0 (liquid) | 11.34 | 0.024401 | 0.024492 |
| collapsed | 0.001 (vapor) | any | ~0.021601 | ~0.021624 |
| fluidized | 0.001 (vapor) | any | ~0.016650 | ~0.016669 |

**Collapsed / liquid detail (most reactive scenario):**

| z_flood (cm) | k-eff | k+2σ |
|---|---|---|
| dry | 0.021637 | 0.021661 |
| 5.67 | 0.025635 | 0.025712 |
| 9.45 | 0.026617 | 0.026707 |
| **11.34 ← peak** | **0.026771** | **0.026857** |
| ≥ 13.23 | plateau ~0.02665–0.02671 | — |

**Key findings:**
- All cases are deeply subcritical: peak k+2σ = 0.02686, three orders of magnitude below 0.95
- Reactivity peaks at partial fill (z_flood ≈ 11 cm), not at full flood — water above the bed acts as a reflector and the gain saturates once the reflector is ~7–8 cm thick
- Plateau at z_flood ≥ 13 cm: further rise has negligible effect
- Collapsed is more reactive than fluidized under liquid flooding (+0.0024 at peak)
- Vapor flooding (0.001 g/cm³) is indistinguishable from the dry baseline at all flood levels — steam at atmospheric pressure provides negligible moderation
