---
group: Materials and furnace hardware
title: OpenMC material definitions and the retort/heater/injector shell
---

## Materials

### What was implemented

`furnace/materials.py` — ten material factory functions, each taking the frozen params dict:

| Function | Material | Density | S(α,β) |
|---|---|---|---|
| `uco_kernel(params)` | HALEU UCO fuel kernel | 10.5 g/cm³ | none |
| `buffer_pyc(params)` | porous carbon buffer | 1.0 g/cm³ | c_Graphite |
| `ipyc(params)` | inner PyC | 1.87 g/cm³ | c_Graphite |
| `sic(params)` | SiC pressure-retention layer | 3.20 g/cm³ | none |
| `opyc(params)` | outer PyC | 1.87 g/cm³ | c_Graphite |
| `graphite_structural(params)` | retort, heater, cone | 1.75 g/cm³ | c_Graphite |
| `process_gas(params)` | H₂/MTS CVD atmosphere | ~2.06×10⁻⁴ g/cm³ | none |
| `water(density)` | light water (flooding sweeps) | variable | none* |
| `air()` | standard dry air | ideal gas | none |

*`c_H_in_H2O` omitted — the ENDF/B-VIII.0 build has empty temperature data for this table; free-gas scattering is used instead.

### How it works

Each factory creates an `openmc.Material`, sets density and composition via `add_element()` with atom-fraction percentages, attaches S(α,β) tables where applicable, and sets temperature. All numeric values come from the frozen params dict — no constants are hardcoded.

**All temperatures fixed at 1200 K** (ENDF/B-VIII.0 has no room-temperature tabulation in this build). The bounding case for moderation is cold, but the library forces 1200 K; this is a known non-conservatism accepted for Stage 0.

**Process gas density** is computed from the ideal gas law at 293.6 K and 101325 Pa (ρ ≈ 2.06×10⁻⁴ g/cm³) for the SiC CVD atmosphere (98 mol% H₂ + 2 mol% MTS).

**Water factory** takes density as an argument to support flooding density sweeps without modifying params.

### Notable assumptions

- UCO: U(C₀.₅O₀.₄) atom fractions, 19.75 wt% U-235; S(α,β) not applied (UCO is not graphite)
- SiC: no S(α,β) applied (`c_SiC` exists in ENDF/B but not credited at this stage)
- Structural graphite boron impurity: 0 ppm (NCS convention)
- U-234 content in HALEU: OpenMC's built-in U-234/U-235 = 0.008 mass ratio approximation (may be inaccurate for cascade-enriched HALEU)

Run the material check with:
```
micromamba run -n triso-env python3 -m furnace.materials
```

## Furnace hardware: retort wall, heater, vacuum gap, water-cooled injector

### What was implemented

`furnace/geometry.py` — five surface constructors and one cell-builder:

- `retort_outer_cylinder(params)` — ZCylinder at r = retort OD / 2 (3.3 cm)
- `outer_cone_surface(params)` — ZCone for the outer wall of the graphite cone (same 30° half-angle as inner cone, apex offset to match outer retort radius at z = z_cone_top)
- `heater_surfaces(params)` — (inner_cyl, outer_cyl, z_bot, z_top) for the heater annulus; bottom flush with retort base
- `injector_surfaces(params)` — five ZCylinders (gas bore, body inner, coolant inner, coolant outer, body outer) plus top/bottom ZPlanes
- `outer_boundary_surfaces(params)` — vacuum boundary surfaces (radial ZCylinder at heater OD, top ZPlane at retort top, bottom ZPlane at injector bottom)
- `furnace_shell_cells(params, graphite, water, process_gas)` — assembles 13 cells for all structural and void regions outboard of the bed interior; returns a dict including `boundary_surfaces`

`params.yaml` — added `dimensions.injector` section (six keys: gas_throat, body_id, coolant_id, coolant_od, body_od, cooled_length). Graphite felt insulation removed from model.

### How it works

`furnace_shell_cells()` builds two nested ZCone surfaces (same 30° half-angle, apices offset for inner/outer retort radius) and three cylindrical surfaces defining 13 non-overlapping cells from r = r_retort_inner outward to the vacuum boundary, and from z = −cooled_length (injector bottom) to z = z_retort_top.

Boundary surfaces are returned from `furnace_shell_cells()` so the model assembly can use them to bound the bed interior cells without creating duplicate surfaces.

**Injector:** four concentric annular cells (gas bore / inner graphite / water coolant / outer graphite) below the cone base. Coolant water at 1.0 g/cm³. The injector gas bore and throat are modelled as void/graphite — no process gas below the retort (conservative: adding H₂ below the bed would raise k-eff slightly).

**Outer vacuum boundary** is placed at the heater OD (no insulation layer — conservative choice).

The function does not cover the bed interior (r < r_retort_inner); that is assembled exclusively in `bed_region()`. The nominal-case model combines the two: `all_cells = bed['cells'] + shell['all']`.
