---
step: 4 — Retort wall, heater element, vacuum gap, and water-cooled injector
title: Full furnace shell geometry — structural graphite, vacuum envelope, and injector coolant
status: complete
---

## What was implemented

- `furnace/geometry.py` — five new surface constructors and one cell-builder function:
  - `retort_outer_cylinder(params)` — ZCylinder at r = retort OD / 2 (3.3 cm)
  - `outer_cone_surface(params)` — ZCone for the outer wall of the graphite cone; apex derived so the surface radius equals r_retort_outer at z = z_cone_top
  - `heater_surfaces(params)` — (inner_cyl, outer_cyl, z_bot, z_top) for the heater annulus; bottom flush with retort base
  - `injector_surfaces(params)` — all five ZCylinders (gas bore, body inner, coolant inner, coolant outer, body outer) plus top/bottom ZPlanes
  - `outer_boundary_surfaces(params)` — three vacuum boundary surfaces (radial ZCylinder at heater OD, top ZPlane at retort top, bottom ZPlane at injector bottom)
  - `furnace_shell_cells(params, graphite, water, process_gas)` — assembles 13 cells for all structural and void regions outboard of the bed interior; returns a dict keyed by component plus `boundary_surfaces` and `all`

- `furnace/materials.py` — removed `graphite_felt_insulation()` function and its reference in `__main__`

- `params.yaml` — three changes:
  - Added `dimensions.injector` section (six keys: `gas_throat`, `body_id`, `coolant_id`, `coolant_od`, `body_od`, `cooled_length`); `load_params()` auto-generates `_cm` suffixed keys
  - Removed `materials.graphite_felt_density_gcc` — felt insulation dropped from the model
  - Removed `materials.graphite_felt_thickness_mm` — same

- `scripts/check_geometry.py` — new script: builds the full shell geometry with a simplified process-gas bed interior, exports XMLs, and runs OpenMC in plot mode to detect overlapping cells or undefined regions; produces two XZ PNGs

## How it works

`furnace_shell_cells()` constructs two nested ZCone surfaces (inner and outer, same 30° half-angle, apices offset to match the inner and outer retort radii at z = z_cone_top) and three cylindrical surfaces (retort outer, heater inner, heater outer). These, combined with six Z-planes at the key axial transitions, define thirteen non-overlapping cells that cover everything from r = r_retort_inner outward to the outer vacuum boundary, and from z = −cooled_length (injector bottom) to z = z_retort_top. The function returns the boundary surface objects in `boundary_surfaces` so the model assembly (step 5) can use them to bound the bed interior cells without creating duplicate surfaces at the same position.

The injector is modelled as four concentric annular cells (gas bore / inner graphite / water coolant / outer graphite) stacked at z ∈ [−6.0, −3.0] cm below the cone base. Gas is assumed to be entirely inside the retort; the injector bore and the throat are both void/graphite with no process-gas fill. The water coolant annulus is the only non-graphite, non-void region in the nozzle zone. The existing simplified throat (r < 3 mm, full-graphite) is superseded by a throat cell in `furnace_shell_cells()`.

## Experimental design

Not applicable — this step adds structural geometry only. No transport calculation is run.

## Design decisions

- **Graphite felt insulation removed at user request.** The outer vacuum boundary is placed at the heater OD (r = 4.7 cm). This is conservative relative to adding an insulation reflector.

- **Heater bottom flush with retort base (z = z_cone_top = 3.8 cm).** Position confirmed by user; no axial offset parameter needed.

- **Outer cone surface derived analytically, not parameterised.** The outer cone has the same 30° half-angle as the inner cone. Its apex is z0_outer = z_cone_top − r_retort_outer / tan(30°) ≈ −1.916 cm, giving r = r_retort_outer at z = z_cone_top. This makes the cone wall exactly 8 mm thick at the retort base — matching the retort cylinder wall — and tapers to ~8 mm at the throat junction.

- **Injector inner graphite is one combined cell (r = 0.3 to 0.8 cm).** The drawing provides `body_id` (r = 0.6 cm) as a machining reference inside the solid graphite sleeve; it is not a material boundary. The combined cell is physically correct and avoids a spurious surface.

- **Nozzle and injector bore are void; no process gas modelled below the retort.** All CVD process gas is assumed to be inside the retort. The throat (r < 3 mm, solid graphite) and the injector gas bore (r < 3 mm, void) model the structural material and empty passage but not flowing gas. This is conservative: replacing void with process gas would add a small amount of hydrogen moderator directly under the bed, raising k-eff slightly.

- **`furnace_shell_cells()` does not cover bed interior (r < r_retort_inner).** That region is the exclusive responsibility of `bed_region()`. The cell-builder function is designed to be combined with bed results in the step 5 model assembly: `all_cells = bed['cells'] + bed['cone_void_cells'] + shell['all']`.

- **Boundary surfaces returned from `furnace_shell_cells()`.** The three vacuum boundary surfaces (radial, top, bottom) are created inside the function and returned in `boundary_surfaces`. The model assembly uses these same objects to bound the bed interior cells, preventing duplicate parallel surfaces at the same axial/radial position.

- **Injector height = h_cooled − h_throat = 6.0 − 3.0 = 3.0 cm.** User confirmed `cooled_length` = 60 mm is the total cooled extent below the cone base (throat + injector body), giving an injector body of the same axial height as the throat.

- **Coolant water at 1.0 g/cm³ (liquid water).** The injector is water-cooled in normal operation; this is a fixed condition, not a swept variable. Water density from the existing `water()` factory in materials.py.

## Assumptions

### Confirmed
- Heater bottom flush with retort base (z = z_cone_top).
- Graphite felt insulation removed from model.
- Injector body material: graphite_structural.
- Cone wall has the same 8 mm wall thickness as the retort cylinder wall (same ZCone half-angle applies to outer surface).
- Injector dimensions from fabrication drawing (six values, all confirmed by user).
- Total cooled length (60 mm) spans throat (30 mm) + injector body (30 mm).

### Defaulted
- **Throat cell: full graphite (r < 3 mm).** Carries over the Step 3 simplification. The gas bore through the throat is process gas in reality but modelled as graphite in this cell. Effect on k-eff is negligible (small volume, no fissile material). Marked `# TODO:` above.
- **Injector gas bore is void.** All process gas is assumed to be inside the retort; no gas material is assigned below the cone base.
- **Coolant water density 1.0 g/cm³.** Liquid water at ambient temperature and normal operating pressure. Sub-cooling is expected at the injector; 1.0 g/cm³ is the correct bound.

### Unconfirmed (`# CONFIRM`)
- `graphite_structural_density_gcc: 1.75` — retort wall, cone wall, heater, throat, and injector graphite all use this density. Placeholder from Step 1; must be confirmed against material certificate for the specific nuclear graphite used.
- `nozzle_height: 30.0` mm — throat axial height; placeholder from Step 3.
- All Step 1–3 `# CONFIRM` flags still apply (kernel diameter, layer thicknesses, etc.).
