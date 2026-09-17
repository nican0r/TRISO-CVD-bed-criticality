---
group: Nominal case
title: Fluidized and collapsed bed states at 95 g bare UCO kernel charge
---

## How to run

```bash
micromamba activate triso-env

# Smoke test (~15 s):
python3 scripts/smoke_test.py

# Nominal fluidized and collapsed cases:
caffeinate python3 scripts/run_nominal.py
caffeinate python3 scripts/run_collapsed.py
```

Results written to `results/step5_nominal/` and `results/step5_collapsed/`.

## What was implemented

- **`furnace/model.py`** — `build_model(params, *, state, stage, background, n_inactive, n_active, n_particles, seed, charge_mass_g)` assembles and returns `(openmc.Model, BedStats)`. Includes three tallies: flux spectrum, material reaction rates, and U-235 fission spatial distribution.
- **`scripts/run_nominal.py`** — fluidized bed eigenvalue run
- **`scripts/run_collapsed.py`** — collapsed bed eigenvalue run
- **`scripts/smoke_test.py`** — pipeline check: 5 g charge, 1+1 batches, 10 000 particles (~15 s)

## How it works

`build_model()` calls `exact_cone_bed()` for the TRISO bed (state-dependent geometry from the bed-geometry group), `furnace_shell_cells()` for the retort shell, then fills the retort interior above the bed with process gas. All three regions combine into a single root universe.

The two nominal cases use identical materials and charge (95 g bare kernel, process gas background) — they differ only in bed state:

- **Fluidized** (pf = 0.333): bed top at z ≈ 3.75 cm, V_bulk ≈ 27.1 cm³ (~96% of frustum)
- **Collapsed** (pf = 0.500): bed top at z ≈ 3.21 cm, V_bulk ≈ 18.1 cm³ (~64% of frustum); upper cone is empty gas

**Settings:** 20 000 particles/generation, 50 inactive + 250 active batches, seed 42.

## Results

Canonical run: `step5_nominal/1789561375_fa29b5e83832`, git `030835a`.

| state | bed_height_cm | pf_achieved | u₂₃₅ mass (g) | k-eff | σ | k+2σ |
|---|---|---|---|---|---|---|
| fluidized | 3.751 | 0.3331 | 17.820 | 0.016647 | 1.20×10⁻⁵ | 0.016672 |
| collapsed | 3.212 | 0.4997 | 17.820 | 0.021637 | 1.20×10⁻⁵ | 0.021661 |

Δk(collapsed − fluidized) = **+4.99×10⁻³** (~70 MC-σ). Both states are deeply subcritical (k+2σ < 0.022 ≪ 0.95) at the 95 g bare-kernel nominal charge in a dry process-gas atmosphere.
