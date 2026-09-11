# Step 6 — Parametric driver

## How to run

```
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env
```

Seed-check sweep (quick, ~30 s total):

```
caffeinate python3 scripts/run_sweep.py --seed-check --quick
```

Full seed-check at nominal 20 000-particle / 300-batch settings:

```
caffeinate python3 scripts/run_sweep.py --seed-check
```

Three cases in parallel (each gets its own run directory; CSV writes are file-locked):

```
caffeinate python3 scripts/run_sweep.py --seed-check --quick --jobs 3
```

Force re-run even if tags exist in CSV:

```
caffeinate python3 scripts/run_sweep.py --seed-check --quick --force
```

Per-case OpenMC thread count:

```
caffeinate python3 scripts/run_sweep.py --seed-check -s 4
```

---

## What was implemented

- **`furnace/sweeps.py`** — full implementation replacing the step-6 stub:
  - `run_case(base_params, overrides, tag, run_dir, *, csv_path, threads, mpi_args, quick, force) → dict` — deep-copies params, applies overrides, builds the model, runs it in an isolated directory, extracts all result quantities, appends to CSV
  - `_unfreeze(obj)` — recursively converts MappingProxyType / tuples to mutable dicts / lists
  - `_freeze(obj)` — re-wraps the modified dict back into a frozen MappingProxyType (local copy; does not call into `furnace.params`)
  - `_merge(base, patch)` — deep-merges a nested override dict into the mutable params dict in-place; underscore-prefixed keys are skipped (reserved for build_model kwargs)
  - `_compute_hc_ratios(stats, params, stage, fill_mat) → (float, float)` — H/²³⁵U and C/²³⁵U atom ratios from actual model materials and volumes; solid-phase contributions from `particle_atom_counts()`, gas/water contributions from OpenMC material atom densities
  - `_thermal_flux_fraction(sp) → float` — thermal / total flux from the `flux_spectrum_bed` tally in a statepoint
  - `_load_dk_disc(n_slabs) → float` — reads δk_disc from `results/convergence_n_slabs.csv` for a given n_slabs; 0.0 if absent
  - `_DK_DISC_N8` — module-level constant (dead code; loaded at import from CSV for n_slabs=8; superseded by `_load_dk_disc` which reads dynamically for any n_slabs)
  - `_append_to_csv(csv_path, row)` — appends one result row; creates header if file is new; uses `fcntl.flock` for cross-process crash safety on macOS/Linux
  - `_export_and_run(model, out_dir, threads, mpi_args) → Path` — exports XML and runs OpenMC with optional shared-memory / MPI parallelism; separate from `model.export_and_run` so it does not touch the running nominal case

- **`scripts/run_sweep.py`** — CLI driver replacing the step-6 stub:
  - `--seed-check` flag: runs the 3-case seed-variation validation sweep
  - `--quick` flag: overrides to 1 000 particles / 12 batches for pipeline smoke-testing
  - `--force` flag: re-runs cases even if tag exists in CSV
  - `-s/--threads N`: OpenMC shared-memory threads per case
  - `-j/--jobs N`: number of cases to run concurrently via `ProcessPoolExecutor`
  - `_worker(args)` — module-level worker function (required for process pickling)
  - `_print_csv(csv_path)` — prints a readable subset of columns after the sweep

---

## How it works

`run_case` is the single composable unit of the sweep infrastructure. It receives a frozen `base_params` object and a nested `overrides` dict (mirroring the params YAML structure), unfreezes the params into a mutable dict, deep-merges the overrides, and re-freezes before calling `build_model`. The re-frozen params are then the single source of truth for that case — the geometry, materials, and settings all read from it, so any override (enrichment, charge mass, gas composition, n_slabs, seed) is automatically propagated everywhere.

Underscore-prefixed override keys (`_state`, `_stage`, `_background`) are not merged into params; they are extracted separately and passed as `build_model` keyword arguments, since those dimensions have no equivalent params.yaml key.

The result dict captures every derived quantity immediately after the run: H/²³⁵U and C/²³⁵U ratios (the x-axes for subsequent sweep plots), thermal flux fraction, wall time, and all k-eff combinations including the discretisation bias. The CSV row is appended with a file lock before the function returns, so a later crash in the calling script cannot lose completed results.

The seed-check sweep in `run_sweep.py` runs three cases with seeds 42, 43, 44 and no other changes. Because `build_model` threads the seed through both the packing RNG and the OpenMC transport RNG, each case produces a genuinely independent realisation. Consistent H/²³⁵U and C/²³⁵U ratios across all three (the gas composition and bed volume are unchanged) confirms the ratio computation is stable; any variation in k-eff across seeds bounds the combined geometric + transport uncertainty.

---

## Design decisions

- **`_freeze` re-implemented locally, not imported from `furnace.params`** — `furnace.params._freeze` is a private symbol. Re-implementing the five-line function avoids a dependency on a private interface and keeps `sweeps.py` self-contained.

- **Underscore-prefix convention for build_model kwargs** — `state`, `stage`, and `background` have no params.yaml equivalents; they cannot be injected via the normal nested-dict merge. Reserving `_state`/`_stage`/`_background` in the overrides dict is the least-surprise API: callers can express a complete case specification in one dict without needing to know which params go to yaml and which go to `build_model`.

- **`overrides_json` as a single CSV column** — The overrides schema varies across sweep types (seed sweep overrides `model.seed`; mass sweep overrides `dimensions.bed.charge_mass_g`; a combined sweep overrides both). A fixed-width CSV that spans all possible override combinations would require adding new columns for each new sweep, breaking existing CSV readers. One JSON column gives full auditability without schema changes.

- **`fcntl.flock` for CSV writes** — When `--jobs N > 1` runs N cases in parallel processes, each process appends to the same CSV. `flock` provides an OS-level exclusive lock that prevents interleaved writes even across separate Python processes. The fallback (`except ImportError: pass`) makes the code run on Windows without crashing, at the cost of possible interleaving if parallel jobs are used there.

- **`_DK_DISC_N8` loaded at module import (dead code)** — Originally loaded δk_disc for n_slabs=8 as a fixed constant. Superseded by `_load_dk_disc(n_slabs)` which reads the CSV dynamically for the actual n_slabs configured in params. `_DK_DISC_N8` is retained for auditability but not used in any live code path.

- **Quick mode: 5 g charge mass + 1 000 particles / 2 inactive + 10 active** — The charge mass is reduced to 5 g (≈12 000 TRISO particles, matching the smoke test) so that the geometry XML stays small and OpenMC can read it in seconds rather than minutes. Without the mass reduction, quick mode still builds the full 95 g / 226 k-particle geometry, which produces a multi-GB XML file that takes 5+ minutes to parse. Transport settings (1 000 particles, 12 batches) produce meaningless k-eff statistics; the goal is only to confirm the pipeline runs end-to-end and the CSV is written correctly.

- **`_export_and_run` in `sweeps.py` rather than extending `model.export_and_run`** — `model.py` is imported by the currently running `run_nominal.py`. Adding optional `threads`/`mpi_args` kwargs would be a backward-compatible change, but modifying any module imported by a live process is risky during the nominal run. The sweep-specific helper avoids touching `model.py` entirely.

- **`ProcessPoolExecutor` for `--jobs`** — Multiple concurrent OpenMC runs each write to separate directories (`results/seed_check/seed_42/`, etc.), so there is no file contention on the XML or statepoint files. Each worker calls `load_params()` independently to avoid passing a `MappingProxyType` across the process boundary (picklability of `MappingProxyType` is CPython-version-dependent). The CSV append lock handles the only shared resource.

---

## Assumptions

**Confirmed:**
- δk_disc is loaded from `results/convergence_n_slabs.csv` by `_load_dk_disc(n_slabs)` for the configured n_slabs. At n_slabs=32, δk_disc=0.0 (no n=64 baseline; estimated missing bias ~3×10⁻⁵, negligible vs σ). Direction is conservative: staircase overestimates k, so adding the bias makes the upper bound more conservative.
- H/²³⁵U and C/²³⁵U ratios are computed bed-only: solid-phase from `particle_atom_counts()` (all TRISO layers present at the given stage), gas-phase from `fill_mat.get_nuclide_atom_densities()` applied to the void volume `V_bulk × (1 − pf_achieved)`
- Nuclide-prefix convention: hydrogen isotopes start with 'H', carbon isotopes start with 'C' — consistent with OpenMC's ENDF/B nuclide naming

**Defaulted:**
- Quick mode settings: 1 000 particles, 12 batches (2 inactive + 10 active). These are only for pipeline validation; they produce statistically unreliable k-eff values.
- Default CSV path: `results/sweep.csv` (general sweeps) and `results/seed_check.csv` (seed-check sweep). Each sweep type gets its own CSV to avoid schema conflicts between different override structures.
- `_state` = `'fluidized'`, `_stage` = `'bare_kernel'`, `_background` = `'gas'` when the corresponding `_` keys are absent from `overrides`.

**Unconfirmed (`# CONFIRM`):**
- The δk_disc bias from the convergence study (collapsed bed, water background, low particle count) is applied conservatively to all cases regardless of bed state or background. For the fluidized / gas-background cases it may overestimate the true discretisation error (see step-3 docs). This is a known Stage 0 approximation.
