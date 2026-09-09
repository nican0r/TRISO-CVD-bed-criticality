---
step: 0 — Scaffold
title: Repository layout, params.yaml, load_params(), check_env()
status: complete
---

## What was implemented

- `params.yaml` — single source of truth for all dimensions, TRISO layer parameters, material densities, gas composition, and model settings; every `# CONFIRM` item flagged.
- `furnace/params.py` — `load_params()` reads the YAML, validates required keys, converts mm→cm and µm→cm into `_cm`-suffixed keys, and returns a recursive `MappingProxyType`; `check_env()` validates the OpenMC library path and prints version/table count.
- `furnace/__init__.py` — package root; re-exports `load_params` and `check_env`.
- `furnace/materials.py`, `triso.py`, `geometry.py`, `model.py`, `sweeps.py`, `postproc.py` — empty stubs with `# TODO` markers for steps 1–9.
- `scripts/run_nominal.py`, `scripts/run_sweep.py` — empty stubs.
- `cases/` (gitignored), `results/` directories.
- `.gitignore` updated to exclude `cases/`, `results/*.{csv,png,h5}`, `statepoint.*.h5`, `summary.h5`.
- `progress.md` created at repo root.

## How it works

`load_params()` reads `params.yaml` with PyYAML, validates that all required top-level and materials keys are present, then applies two unit conversions: every linear key in the `dimensions` section that contains a dimension-related substring (diameter, height, length, etc.) gets a `_cm` companion key at `value × 0.1`; every key in the `triso` section gets a `_cm` companion at `value × 1e-4`. The entire nested structure is then wrapped in `types.MappingProxyType` recursively, making it read-only at every level. `check_env()` reads `OPENMC_CROSS_SECTIONS`, opens the XML via `openmc.data.DataLibrary.from_xml()`, and prints version and table count; `import openmc` is deferred to call time so the module is importable without OpenMC.

## Design decisions

- **Return type: recursive `MappingProxyType` (frozen dict), not nested dataclasses.** The parameter schema grows with each step; dataclasses would require updating field definitions for every new parameter. `MappingProxyType` gives the same mutation protection with zero per-field boilerplate. Tradeoff: no IDE key autocompletion for the returned object.

- **`_cm` conversion keyed on name substrings, not on all numeric values.** `_add_cm_keys()` only adds a `_cm` companion for keys whose names contain "diameter", "depth", "height", "length", "drop", "thickness", "id", or "od". This prevents `included_angle_deg`, `pressure_pa`, `temperature_k`, and `static_bulk_cc` from receiving a nonsensical `× 0.1` or `× 1e-4` result.

- **`static_bulk_cc` stored under `dimensions.bed_volume`, not `dimensions.bed`.** It is in cm³ (not mm); mixing it into a mm-dimensioned sub-dict would require a special-case guard in the converter. The separate sub-dict makes the unit difference explicit in the YAML.

- **`check_env()` defers `import openmc`.** `furnace.params` is importable in environments without OpenMC (unit tests, CI parameter checks). Importing at module level would break any `from furnace import load_params` in such environments.

- **`_cm` keys added alongside originals, not replacing them.** `params['triso']['kernel_diameter']` → 425.0 µm (natural, matches YAML); `params['triso']['kernel_diameter_cm']` → 0.0425 cm (OpenMC-ready). Both are correct; no hidden state.

- **`DataLibrary.from_xml()` used, not `from_hdf5()`.** OpenMC 0.15.3 uses `from_xml()`; `from_hdf5()` was removed. Noted here so a future upgrade check knows this was tested.

- **Seed = 42 in `model.seed`.** Arbitrary but deterministic and recorded; must appear in every result file per the preamble working agreement.

## Assumptions

### Confirmed
- All material densities and enrichment match the reference repo (`nican0r/TRISO-neutronics`) and AGR-1 specs (Demkowicz et al., NED 329 (2018) 102–111; INL/EXT-10-19476).
- Gas: 98 mol% H₂ + 2 mol% CH₃SiCl₃ at 101325 Pa / 1873 K — user-confirmed as worst-case moderation scenario; held constant across all CVD steps.
- Structural graphite boron impurity = 0 ppm (no unconfirmed poison credited).
- Nuclear data library: ENDF/B-VIII.0 processed at 900 K and 1200 K; path set in `OPENMC_CROSS_SECTIONS`. 690 nuclide/thermal tables confirmed present.

### Defaulted
- TRISO layer thicknesses (buffer 100 µm, IPyC 40 µm, SiC 35 µm, OPyC 40 µm): AGR-1 nominal from INL/EXT-10-19476 Table 3. `# CONFIRM` flags retained.
- Kernel diameter 425 µm: AGR-1 nominal from Demkowicz et al. 2018. `# CONFIRM` retained.
- Structural graphite density 1.75 g/cm³: from AGR-1 compact matrix (INL/EXT-10-19476 Table 5); may not match the furnace graphite grade. `# CONFIRM`.
- Graphite felt density 0.20 g/cm³: placeholder from preamble. `# CONFIRM`.
- Placeholder simulation settings: 1000 particles, 110 batches (10 inactive). These are for development; production runs need ≥ 50 000 particles.
- Seed = 42: arbitrary.

### Unconfirmed (`# CONFIRM`)
- `bed.static_depth` (35 mm) — wrong value shifts fissile mass distribution in cone geometry; directly affects volume fraction of cone vs. cylinder.
- `bed_volume.static_bulk_cc` (50 cm³) — sets total TRISO particle count; scales directly with k-eff.
- `bed.fluidized_height` (180 mm) — needed for the expanded-bed geometry in steps 3–4.
- `model.fuel_temperature_k` (1500 K) — library only has 900 K and 1200 K; effective XS temperature will be 1200 K via nearest-temperature lookup regardless of this value.
- `graphite_felt_thickness_mm` (null) — insulation not yet modelled; if thick and close to the fissile region, moderation contribution may matter.
- `gas.pressure_pa` (101325 Pa) — sub-atmospheric operation is possible; lower pressure → lower gas density → less H moderation → lower k-eff. Flag if sub-atm cases are needed.
