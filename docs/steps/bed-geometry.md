---
group: Bed geometry
title: TRISO particle universe, random packing, and exact-cone bed geometry
---

## TRISO particle and packed-bed universe

### What was implemented

`furnace/triso.py`:

- `STAGES` — canonical deposition-stage names: `bare_kernel`, `buffered`, `ipyc`, `sic`, `full_triso`
- `particle_at_stage(stage, params)` — returns `(Universe, outer_radius_cm)` with one concentric-sphere cell per deposited layer
- `pack_bed(region, packing_fraction, outer_radius, fill_universe, seed, params)` — calls `openmc.model.pack_spheres`, wraps results as `openmc.model.TRISO` objects, caches sphere centres to `cases/.triso_cache/<sha256[:16]>.npz`
- `tiled_bed(pf, outer_radius, fill_universe, seed, tile_size_cm, params)` — packs one tile cube once, wraps it in a `Universe`, broadcasts across a `RectLattice`; used in production to keep memory O(particles per tile) regardless of charge mass
- `lattice_bed(trisos, lower_left, upper_right, background_material, pitch_target)` — wraps `openmc.model.create_triso_lattice`
- `bed_stats(trisos, region_volume, stage, params)` — returns particle count, achieved PF, HM mass, U-235 mass, smeared density, and C/U-235 and Si/U-235 atom ratios (all analytic, no OpenMC library needed)

`params.yaml` additions:
- `max_packing_fraction: 0.64` — hard upper guard
- `packing_fraction_static: 0.50` — collapsed bed
- `packing_fraction_fluidized: 0.333` — derived from pf_static / bed_expansion_ratio (1.5)

### How it works

`particle_at_stage` builds a fresh `openmc.Universe` with one cell per deposited layer. `pack_bed` places non-overlapping sphere centres via `openmc.model.pack_spheres` and caches results using a SHA-256 key on the bounding box, packing fraction, radius, and seed — skipping the expensive placement step on repeat calls.

`tiled_bed` is the production packing path: it packs one ~0.5 cm cube tile once and places the same universe reference at every cell of a `RectLattice`. Unique cell count is O(particles per tile) rather than O(total particles), enabling the mass sweep up to 3610 g without OOM.

### Validation results (AGR-1 nominal dimensions)

Layer outer radii:

| Layer | r (cm) |
|---|---|
| Kernel | 0.02125 |
| Buffer | 0.03125 |
| IPyC | 0.03525 |
| SiC | 0.03875 |
| OPyC | 0.04275 |

C/U-235 atom ratios by stage:

| Stage | C/U-235 | Si/U-235 |
|---|---|---|
| bare_kernel | 2.51 | 0.00 |
| buffered | 24.16 | 0.00 |
| ipyc | 49.86 | 0.00 |
| sic | 64.13 | 14.27 |
| full_triso | 102.72 | 14.27 |

Run the diagnostic with:
```
micromamba run -n triso-env python3 -m furnace.triso
```

## Exact-cone bed geometry with state-dependent bed height

### What was implemented

`furnace/geometry.py`:

- `exact_cone_bed(params, state, stage, background_material, charge_mass_g, seed, z_flood, dry_material, *, use_tiled_bed, tile_size_cm)` — assembles the truncated-frustum bed cell at a state-dependent `z_bed_top`, packs it at `pf_state`, adds an above-bed void cell, and handles overflow into the retort cylinder
- `frustum_volume`, `_frustum_volume_below_z`, `_solve_bed_top_z` — analytic frustum volume and bisection solver for the bed-top plane
- `retort_inner_cylinder`, `cone_wall`, `nozzle_throat` — surrounding surfaces

This is the production geometry used by the nominal case, mass sweep, and flooding sweeps.

### How it works

The bed is a single OpenMC cell whose region is the truncated frustum (`-cone_surf & +zp_cone_bot & -zp_bed_top`) filled by a broadcast tile lattice.

**State-dependent bed top:** Given `V_solid = charge_mass_g / ρ_particle` and the state's packing fraction, the required bed volume is `V_bulk = V_solid / pf_state`. The code bisects for `z_bed_top` such that the frustum volume below it equals `V_bulk`. The bed cell is truncated there and the tile is packed at exactly `pf_state`. Fissile inventory is conserved by construction.

- **Collapsed** (pf = 0.50): shorter, denser column at the cone base; upper cone is empty gas
- **Fluidized** (pf = 0.333): fills more of the cone at lower density

When `V_bulk > V_frustum`, the excess spills into a retort-cylinder overflow cell (exercised in the mass sweep).

**Tiled fill:** `tiled_bed()` packs one ~0.5 cm cube tile once (cached), wraps it in a Universe, and broadcasts it across a `RectLattice`. Unique cell count is O(particles per tile), not O(charge mass).

### Convergence study: staircase vs. exact cone

A prior inscribed-staircase geometry (n equal-height cylinders approximating the cone) was tested and retired. Results at 95 g bare_kernel, full-flood, collapsed:

| config | mean k-eff | bias vs exact-cone | |Δk|/σ_seed |
|---|---|---|---|
| exact_cone | 0.02417 | — | — |
| n=32 slabs | 0.02648 | +2.31×10⁻³ | 22.4 |
| n=16 slabs | 0.02712 | +2.95×10⁻³ | 51.5 |
| n=8 slabs | 0.02745 | +3.28×10⁻³ | 33.8 |

Every staircase configuration carries a statistically significant +2–3×10⁻³ k-eff bias. No finite n_slabs meets the convergence criterion. **Exact-cone tiled geometry is the production choice.**

Data: `results/step3_convergence/` and `results/step3_convergence_exactcone_tiled/`.
