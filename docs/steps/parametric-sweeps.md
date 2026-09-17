---
group: Parametric sweeps
title: Sweep driver, seed-variation validation, and charge-mass sweep
---

## Parametric driver

### How to run

```bash
micromamba activate triso-env

# Seed-check sweep (quick, ~30 s total):
caffeinate python3 scripts/run_sweep.py --seed-check --quick

# Full seed-check at production settings:
caffeinate python3 scripts/run_sweep.py --seed-check

# Three cases in parallel:
caffeinate python3 scripts/run_sweep.py --seed-check --quick --jobs 3
```

### What was implemented

- **`furnace/sweeps.py`** — `run_case(base_params, overrides, tag, run_dir, ...)` deep-copies params, applies nested override dict, builds and runs the model in an isolated directory, appends results to a crash-safe locked CSV
- **`scripts/run_sweep.py`** — CLI driver with `--seed-check`, `--mass-sweep`, `--quick`, `--force`, `--threads`, `--jobs` flags

### How it works

`run_case` is the composable unit: it receives a frozen `base_params` and a nested `overrides` dict, unfreezes params, deep-merges the overrides, re-freezes, then calls `build_model`. Underscore-prefixed keys (`_state`, `_stage`, `_background`) are extracted and passed directly to `build_model` rather than merged into params.

CSV rows are appended with `fcntl.flock` so parallel jobs writing to the same file don't corrupt it.

The **seed-check sweep** runs three cases with seeds 42, 43, 44, all other parameters at nominal. Because the seed threads through both the packing RNG and transport RNG, each case is an independent realisation. Consistent H/²³⁵U and C/²³⁵U ratios across seeds confirms the ratio computation is stable; k-eff spread bounds combined geometric + transport uncertainty.

**Quick mode** reduces charge to 5 g and settings to 1 000 particles / 12 batches — for pipeline smoke-testing only, not meaningful k-eff.

Outputs: `results/seed_check.csv`, per-case run dirs under `results/seed_check/`.

## Mass sweep

k-eff vs. charge mass for the collapsed bare-kernel bed (95 g – 3610 g).

### How to run

```bash
micromamba activate triso-env

# Quick smoke-test (~30 s/case):
caffeinate python scripts/run_sweep.py --mass-sweep --quick

# Full sweep:
caffeinate python scripts/run_sweep.py --mass-sweep -j 3

# Plot:
caffeinate python scripts/plot_mass_sweep.py
```

Or submit to AWS Batch (see `docs/aws-batch.md`):
```bash
python scripts/batch/submit_sweep.py --sweep step7a_mass --manifest manifests/step7a_mass.json --dry-run
python scripts/batch/submit_sweep.py --sweep step7a_mass --manifest manifests/step7a_mass.json
python scripts/batch/pull_results.py --sweep step7a_mass
```

Outputs: `results/mass_sweep.csv`, `results/mass_sweep.png`, per-case run dirs under `results/mass_sweep/`.

### What was implemented

- **`furnace/sweeps.py`** — `vessel_capacity_g(params, state, stage)` computes the maximum admissible charge mass from retort usable volume × packing fraction × particle density; cases above this are dropped before submission
- **`scripts/run_sweep.py`** — `--mass-sweep` flag; sweeps multipliers `[1.5×, 2×, 3×, 5×, 10×, 20×, 38×]` of the 95 g nominal charge
- **`scripts/plot_mass_sweep.py`** — log-x k-eff ± 2σ vs charge mass with the 0.95 subcritical reference line

### How it works

Each case runs the nominal model at `state='collapsed'`, `stage='bare_kernel'`, `background='gas'` with only the charge mass overridden. The collapsed bed (pf = 0.50) is the most reactive gravity-only configuration — it concentrates the same fissile mass into a smaller volume than the fluidized state.

The bed grows naturally from the cone base into the retort cylinder as mass increases. Vessel capacity at pf = 0.50 is ~3 653 g; the grid terminates at 38× (3 610 g).

**Settings:** 20 000 particles/generation, 300 batches (50 inactive + 250 active), seed 42 — matching the nominal case.

### Results

Canonical run: `step7a_mass/1789561375_788a31460970`, git `030835a`, manifest `manifests/step7a_mass_tiled.json`.

| charge_mass_g | bed_height_cm | u₂₃₅ mass (g) | k-eff | k+2σ | wall (s) |
|---|---|---|---|---|---|
| 95 | 3.21 | 17.82 | 0.021637 | 0.021661 | 207 |
| 142.5 | 3.75 | 26.73 | 0.024794 | 0.024820 | 238 |
| 190 | 4.20 | 35.66 | 0.027437 | 0.027465 | 279 |
| 285 | 5.13 | 53.48 | 0.032148 | 0.032181 | 334 |
| 475 | 6.97 | 89.12 | 0.038874 | 0.038906 | 370 |
| 950 | 11.58 | 178.22 | 0.047674 | 0.047719 | 463 |
| 1900 | 20.79 | 356.42 | 0.054148 | 0.054199 | 517 |
| 3610 | 37.38 | 677.18 | 0.057324 | 0.057379 | 549 |

**Key findings:**
- All cases deeply subcritical: max k+2σ = 0.0574 at 3 610 g, far below 0.95
- k-eff increases monotonically with mass and flattens significantly above ~500 g (doubling from 1 900 g to 3 610 g raises k by only +0.003)
- Bed overflows from cone into retort cylinder at ≥ 190 g (~2× nominal)
- pf_achieved = 0.4997–0.4999 consistently across all masses
