# TRISO-CVD-bed-criticality

OpenMC criticality safety (NCS) screening model of a laboratory-scale fluidized-bed chemical vapor deposition furnace used to coat TRISO fuel kernels. The model evaluates whether HALEU UCO TRISO particles accumulating inside a graphite retort during the CVD coating process can approach criticality under normal and off-normal (flooding) conditions.

## Quickstart

### Prerequisites
- Python 3.10, OpenMC (built with MPI/HDF5), and `openmc-data`-compatible libraries on `PATH`.
- Nuclear data pointed to by `OPENMC_CROSS_SECTIONS` (ENDF/B-VIII.0 at 900/1200 K is what the model expects). See [`.claude/CLAUDE.md`](.claude/CLAUDE.md) for the exact paths used in development.
- Repo params live in [`params.yaml`](params.yaml); the runners load them via `furnace.params.load_params()`.
- On macOS wrap long OpenMC runs with `caffeinate` to prevent sleep.

### Smoke test (verify the environment)
```bash
python scripts/smoke_test.py           # tiny geometry, <2 min, catches XS/geometry errors
```

### Step 5 — nominal k-eff cases
```bash
caffeinate python scripts/run_nominal.py     # fluidized bed
caffeinate python scripts/run_collapsed.py   # collapsed bed
```
Results land under `results/step5_nominal/` and `results/step5_collapsed/`.

### Step 3 — cone-slab convergence study
```bash
python scripts/run_convergence.py --local-check   # build one model locally, no run
python scripts/run_convergence.py                 # write manifest, print submit command
python scripts/run_convergence.py --submit        # write manifest and submit to Batch
```

### Steps 6–7 — parametric sweeps (local)
```bash
caffeinate python scripts/run_sweep.py --seed-check --quick        # 3-case smoke sweep
caffeinate python scripts/run_sweep.py --seed-check                # full seed check
caffeinate python scripts/run_sweep.py --mass-sweep --jobs 4       # mass sweep, 4 in parallel
```
CSV output is appended to `results/*.csv`; re-runs skip completed rows unless `--force`.

### Step 8 — flooding sweep
```bash
caffeinate python scripts/run_flood.py --quick   # 1 000 particles, 12 batches per case
caffeinate python scripts/run_flood.py           # full 80-case sweep
```

### AWS Batch sweeps
Every submission must be traceable to a manifest under [`manifests/`](manifests). Always dry-run first.
```bash
python scripts/batch/submit_sweep.py --sweep <name> --manifest manifests/<name>.json --dry-run
python scripts/batch/submit_sweep.py --sweep <name> --manifest manifests/<name>.json
python scripts/batch/pull_results.py --sweep <name>      # → results/<sweep>/<submission>/summary.csv
```
See [`docs/aws-batch.md`](docs/aws-batch.md) for the CloudFormation stack, queue, and job-definition details, and [`job-status.md`](job-status.md) for the current sweep tracker.

### Plotting results
```bash
python scripts/plot_mass_sweep.py                             # k+2σ vs mass (auto-detects CSV schema)
python scripts/plot_mass_sweep.py --csv results/step7a_mass_tiled/<sub>/summary.csv
python scripts/plot_flood_sweep.py                            # k+2σ+δk vs z_flood
python scripts/plot_staircase.py                              # cone-slab discretisation diagnostics
python scripts/plot_tiling_coverage.py                        # bed-tiling coverage check
```

### Packing cache (optional, speeds up repeated runs)
```bash
python scripts/warm_packing_cache.py            # pre-bakes fluidized/collapsed 95 g configs
python scripts/warm_packing_cache_parallel.py   # same, parallel
```

## Documentation

The model is documented as seven grouped topics under `docs/steps/`:

| Group | Summary | Details |
|-------|---------|---------|
| Preamble | Physical description, key modeling choices, non-conservatisms | [docs/steps/preamble.md](docs/steps/preamble.md) |
| Scaffold | Repo layout, `params.yaml`, `load_params()`, `check_env()` | [docs/steps/scaffold.md](docs/steps/scaffold.md) |
| Materials & furnace hardware | OpenMC material factories (kernel, coatings, graphite, gas, water, air) plus the retort/cone walls, heater annulus, vacuum gap, and water-cooled injector | [docs/steps/materials-and-hardware.md](docs/steps/materials-and-hardware.md) |
| Bed geometry | TRISO single-particle universe, random/tiled packing with disk cache, exact-cone frustum bed, state-dependent bed height, staircase-vs-exact-cone convergence study | [docs/steps/bed-geometry.md](docs/steps/bed-geometry.md) |
| Nominal case | Full `openmc.Model` assembly, tallies, fluidized/collapsed baseline runs at 95 g bare kernel | [docs/steps/nominal-case.md](docs/steps/nominal-case.md) |
| Parametric sweeps | `run_case` driver (deep-copy overrides, crash-safe CSV, parallel jobs), seed-check validation, and 1× → 38× charge-mass sweep with vessel-capacity guard and k+2σ vs mass plot | [docs/steps/parametric-sweeps.md](docs/steps/parametric-sweeps.md) |
| Flooding accident cases | Bottom-up flood sweep: 2 bed states × 20 z_flood levels × 2 water densities (liquid/vapor); peak-reactivity results and top-10 case table | [docs/steps/flooding-accident.md](docs/steps/flooding-accident.md) |
