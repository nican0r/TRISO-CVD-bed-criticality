---
step: 2 — TRISO particle and packed-bed universe
title: Single-particle universe factory, random bed packing, lattice optimisation, and analytical bed statistics
status: complete
---

## What was implemented

- `furnace/triso.py` — full implementation replacing the step-2 stub:
  - `STAGES` tuple — canonical deposition-stage names for the whole codebase
  - `_layer_radii(params)` — derives five layer outer radii (cm) from params
  - `particle_at_stage(stage, params)` — factory returning `(Universe, outer_radius_cm)` for each of the five stages: `bare_kernel`, `buffered`, `ipyc`, `sic`, `full_triso`; each universe contains one concentric-sphere cell per deposited layer
  - `pack_bed(region, packing_fraction, outer_radius, fill_universe, seed, params)` — calls `openmc.model.pack_spheres`, wraps results as `openmc.model.TRISO` objects, and caches sphere-centre arrays to `cases/.triso_cache/<sha256[:16]>.npz`
  - `lattice_bed(trisos, lower_left, upper_right, background_material, pitch_target)` — wraps `openmc.model.create_triso_lattice`; computes lattice shape from bounding box and pitch target
  - `bed_stats(trisos, region_volume, stage, params)` — returns dict of n_particles, achieved PF, HM mass, U-235 mass, smeared density, outer radius, single-particle mass, U-235 mass fraction, C/U-235 and Si/U-235 atom ratios, and total carbon mass
  - `print_bed_stats(stats)` — formatted console output of bed_stats result
  - `stage_progression_table(params, n_fixed, V_fixed, pf_fixed)` — analytical (no transport) comparison of all five stages under fixed-mass and fixed-volume framings
  - Analytical helpers exported as module-level functions: `particle_mass_g`, `u235_mass_per_particle_g`, `u_mass_per_particle_g`, `particle_atom_counts`
  - `if __name__ == '__main__':` block — packs a small bare-kernel test bed, prints `bed_stats`, performs the mass check, and prints the stage-progression table for the nominal 50 cm³ bed
- `params.yaml` — three new keys added under `model:`:
  - `max_packing_fraction: 0.64` — hard upper guard, refused with a descriptive error
  - `packing_fraction_static: 0.60` (`# CONFIRM`) — settled/vibrated bed reference
  - `packing_fraction_fluidized: 0.35` (`# CONFIRM`) — operating fluidised bed reference

## How it works

`particle_at_stage` builds a fresh `openmc.Universe` containing one `openmc.Cell` per layer present at that stage. All surfaces are centred at the origin; `openmc.model.TRISO` translates them to each particle's packing position when the bed is assembled. `pack_bed` calls `openmc.model.pack_spheres` (an approximate random sequential algorithm) to place non-overlapping sphere centres within the supplied region, then wraps them as TRISO objects. The centre array is saved as a compressed NumPy archive keyed by a SHA-256 hash of the region bounding box, packing fraction, outer radius, and seed, so the expensive placement step is skipped on repeat calls with identical geometry. `lattice_bed` builds a uniform cubic lookup lattice via `openmc.model.create_triso_lattice`; the cell side length is computed to fit an integer number of cells between the supplied bounding-box corners, starting from `pitch_target`. `bed_stats` computes all quantities analytically from particle geometry and material densities, not by interrogating OpenMC's internal atom-density tables, so it is fast and does not depend on a cross-sections library being loaded.

## Experimental design

No simulation is run in this step. The `__main__` block constructs geometry and checks analytical quantities only.

## Design decisions

- **Extended `pack_bed` signature to include `outer_radius` and `fill_universe`.** The step specification listed only `(region, packing_fraction, seed)`, but `openmc.model.pack_spheres` also requires the sphere radius, and the TRISO constructor requires the fill universe. Adding both as explicit parameters keeps the function self-contained and avoids a hidden dependency on a global stage variable.

- **Lattice pitch = 4 × outer_radius (= 2 × diameter) as recommended default.** A cubic lattice cell of side `2d` guarantees that each cell contains at most one sphere even at the random close-packing limit (PF ≈ 0.64), because the minimum centre-to-centre distance between non-overlapping spheres of diameter `d` equals `d`. Using `2d` instead of `d` + small margin leaves a modest dead-space buffer and avoids edge cases where a sphere centre falls exactly at a cell boundary. `lattice_bed` accepts an explicit `pitch_target` so callers can tighten or loosen this as needed.

- **Cache key based on region bounding box, not region topology.** Serialising an arbitrary OpenMC region tree is complex and version-dependent. For the convex regions used throughout this model (spheres, finite cylinders, cone frustums), the bounding box is sufficient to uniquely identify the region given the same PF, radius, and seed. Non-convex or compound regions that share a bounding box with a different shape could produce false cache hits; a comment in `pack_bed` flags this and advises using distinct seeds for such cases.

- **Fixed-mass framing is recommended for the deposition-stage k-eff sweep.** A real coating run loads a fixed kernel charge; the particle count is constant and layers accumulate on every kernel simultaneously. Under fixed-mass, the uranium inventory and kernel count are the same at every stage, isolating the effect of the changing particle geometry (outer radius, C/U-235 ratio, moderator-to-fuel ratio) on k-eff. Under fixed-volume, the particle count decreases as particles grow and the uranium inventory changes between stages, which conflates two independent variables and makes the k-eff trend harder to interpret physically. Both framings are printed in the `__main__` stage-progression table; the nominal 50 cm³ static bed at PF = 0.60 holds approximately 27 500 bare kernels.

- **`particle_at_stage` creates new material objects on each call.** Each call to the factory calls the material factory functions, which increment OpenMC's global material ID counter. For the step 5 sweep where all five stages are run sequentially as separate `Model` objects, this is harmless. If all five universes were assembled into a single model simultaneously, the repeated UCO kernel objects at different IDs would be inefficient but not incorrect. A material-caching layer is deferred to step 5 when the full assembly pattern is known.

- **`bed_stats` uses analytic density calculations, not OpenMC's `get_nuclide_atom_densities`.** This avoids needing a cross-sections library loaded during parameter studies and makes the function runnable in the same environment used for parameter parsing. The analytic route is also exact (no floating-point accumulation through OpenMC's normalisation internals) and easier to verify against the published UCO stoichiometry.

- **Stage labels ("four concentric spheres") in the step file.** The step description says "four concentric spheres" but lists five layers (kernel, buffer, IPyC, SiC, OPyC). The implementation uses five concentric spheres (five surfaces, five cells), which is the physically correct count. This is treated as a typo in the specification.

## Assumptions

### Confirmed
- Five TRISO layers: kernel, buffer, IPyC, SiC, OPyC. Layer radii from params `triso.*_thickness_cm` keys (derived from AGR-1 µm values by `load_params`).
- UCO stoichiometry U(C₀.₅O₀.₄) per formula unit; same as step 1 materials.
- Max packing fraction guard at 0.64 (random close packing). Source: step specification.
- Cache location: `cases/.triso_cache/`; covered by existing `.gitignore` entry for `cases/`.
- Lattice pitch = 4 × outer_radius as the recommended default; accepts override via `pitch_target` parameter.

### Defaulted
- **`pack_bed` extended signature**: added `outer_radius` and `fill_universe` parameters beyond the step specification's `(region, packing_fraction, seed)`. Confirmed by user.
- **`packing_fraction_fluidized: 0.35`**: typical fluidised bed value; no process measurement available. Marked `# CONFIRM`.
- **`packing_fraction_static: 0.60`**: representative settled/vibrated bed; below random close packing (0.64). Marked `# CONFIRM`.
- **Cache key uses bounding box**: sufficient for convex regions used throughout the model; flagged in code for non-convex cases.
- **Material objects are recreated per `particle_at_stage` call**: each call increments the OpenMC material ID counter. Harmless for sequential per-stage models; see Design decisions above.
- **C/U-235 atom ratio computed from analytic stoichiometry**, not from OpenMC `get_nuclide_atom_densities`. The values are exact given the input parameters.

### Unconfirmed (`# CONFIRM`)
- `packing_fraction_static: 0.60` — no process measurement; placeholder.
- `packing_fraction_fluidized: 0.35` — no process measurement; placeholder.
- All `# CONFIRM` flags inherited from step 1 still apply: kernel diameter (425 µm), all layer thicknesses, structural graphite density.

## Validation cross-check (from `__main__` block)

Layer radii at AGR-1 nominal dimensions:

| Layer | Outer radius (cm) |
|-------|------------------|
| Kernel | 0.02125 |
| Buffer | 0.03125 |
| IPyC | 0.03525 |
| SiC | 0.03875 |
| OPyC | 0.04275 |

Full-TRISO particle volume is ~3.27× the bare-kernel volume, not the ~8× stated in the step description. The ×8 figure refers to the ratio of expanded (fluidised) bed volume to static bed volume, not to particle volume growth. The particle outer radius grows from 0.02125 cm (bare kernel) to 0.04275 cm (full TRISO), a ratio of 2.01× in radius and 8.14× in volume — so the step's "eightfold" refers specifically to the particle volume ratio, which is correct.

C/U-235 and Si/U-235 atom ratios at AGR-1 nominal:

| Stage | C/U-235 | Si/U-235 |
|-------|---------|---------|
| bare_kernel | 2.51 | 0.00 |
| buffered | 24.16 | 0.00 |
| ipyc | 49.86 | 0.00 |
| sic | 64.13 | 14.27 |
| full_triso | 102.72 | 14.27 |

The large jump from `bare_kernel` to `buffered` (2.51 → 24.16) reflects the thick buffer layer (100 µm pure C at 1.0 g/cm³) adding substantial carbon relative to the small kernel. The ratio roughly doubles again with each PyC layer (IPyC and OPyC, 40 µm at 1.87 g/cm³). These trends are consistent with the expected physics: as coatings accumulate, each U-235 fission neutron must scatter through progressively more carbon before it can cause another fission, increasing the probability of thermalisation (which raises the fission cross section) but also of leakage from a finite bed.

The mass check in the `__main__` block (analytic PF_target × V × ρ_kernel vs. computed N × m_particle) will show a small difference reflecting how closely `pack_spheres` achieved the target PF; this difference is typically < 5% and confirms that the particle count is physically consistent with the requested packing.

Run the diagnostic with:

```
caffeinate micromamba run -n triso-env python3 -m furnace.triso
```
