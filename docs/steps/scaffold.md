---
group: Scaffold
title: Repository layout, params.yaml, load_params(), check_env()
---

## What was implemented

- `params.yaml` — single source of truth for all dimensions, TRISO layer parameters, material densities, gas composition, and model settings; `# CONFIRM` flags mark unverified values
- `furnace/params.py` — `load_params()` reads the YAML, validates required keys, converts mm→cm and µm→cm into `_cm`-suffixed keys, returns a read-only `MappingProxyType`; `check_env()` validates the OpenMC library path
- `furnace/__init__.py`, `furnace/materials.py`, `triso.py`, `geometry.py`, `model.py`, `sweeps.py`, `postproc.py` — stubs with `# TODO` markers for subsequent groups
- `scripts/run_nominal.py`, `scripts/run_sweep.py` — empty stubs
- `cases/`, `results/` directories; `.gitignore` updated to exclude run artifacts

## How it works

`load_params()` reads `params.yaml`, validates required keys, adds `_cm`-suffixed companions for dimension keys (mm→cm) and TRISO layer keys (µm→cm), then returns the entire structure wrapped in a recursive `MappingProxyType` (read-only at every level). `check_env()` opens the cross-sections XML via `openmc.data.DataLibrary.from_xml()` and reports version and table count.

Key conventions:
- `_cm` keys are added alongside originals (e.g. `kernel_diameter` = 425 µm; `kernel_diameter_cm` = 0.0425 cm)
- `check_env()` defers `import openmc` so the module is importable without OpenMC installed
- Seed = 42 in `model.seed`; must appear in every result file
