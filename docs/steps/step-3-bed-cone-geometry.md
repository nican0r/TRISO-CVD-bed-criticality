---
step: 3 — Bed and cone geometry
title: Exact-cone tiled bed geometry, state-dependent bed top, two bed states (fluidized/collapsed)
status: complete — exact-cone tiled adopted as production 2026-09-16; state-dependent bed top added 2026-09-16 (supersedes mass-conserving tile pf); staircase retired
---

## What is implemented (production)

The bed is a single OpenMC cell whose region is the truncated frustum
(`-cone_surf & +zp_cone_bot & -zp_bed_top`) filled by a broadcast tile
lattice. The bed top plane `z_bed_top` is state-dependent: it is chosen so
that the frustum volume below it equals `V_bulk = V_solid / pf_state`. The
tile inside is packed at `pf_state` directly. Together this gives a
physically distinct geometry for the two bed states — the collapsed state
occupies a shorter, denser column at the base of the cone; the fluidized
state occupies more (or all) of the cone at lower density. Fissile inventory
is conserved: `n_particles = charge_mass_g / m_particle` in both states.

Unique-cell count is O(particles per tile), independent of charge mass, so
the same geometry runs at 5 g and at 3610 g without memory scaling.

Entry point: `build_model(..., use_exact_cone=True, use_tiled_bed=True, tile_size_cm=0.5)`.
This is the geometry used by steps 3, 5, 7 and 8.

- `furnace/geometry.py`:
  - `exact_cone_bed(params, state, stage, background_material, charge_mass_g,
    seed, z_flood, dry_material, *, use_tiled_bed, tile_size_cm)` — assembles
    the truncated-frustum cell at `z_bed_top`, packs it at `pf_state`, adds
    an above-bed void cell inside the cone, and handles overflow into the
    retort cylinder when `V_bulk > V_frustum`.
  - `frustum_volume(r_top, r_bot, half_angle_deg)` — analytic frustum volume.
  - `_frustum_volume_below_z(z, r_bot, tan_theta)` — volume from cone base
    up to height z (used by the bed-top solver).
  - `_solve_bed_top_z(V_target, r_bot, tan_theta, z_cap)` — bisection for the
    z-plane at which the frustum-below volume equals `V_target`.
  - `retort_inner_cylinder(params)`, `cone_wall(params)`,
    `nozzle_throat(params)` — surrounding surfaces.
  - `bed_region(...)` — legacy staircase implementation, retained for
    reproducibility of pre-2026-09-16 results only. See "Historical:
    staircase approach (retired)" at the bottom.
- `furnace/triso.py::tiled_bed` — packs one 0.5 cm tile at a caller-specified
  pf, wraps it in a `Universe`, broadcasts it across a `RectLattice` covering
  the requested bounding box.
- `params.yaml` — nozzle_height, charge_mass_g, static_bulk_cc,
  bed_expansion_ratio, packing_fraction_fluidized (see Design decisions).

## How it works

1. **Cell region is the truncated cone.**
   `frustum_region = -cone_surf & +zp_cone_bot & -zp_bed_top`. The `ZCone`
   is the actual radial boundary; `zp_bed_top` is a state-dependent ZPlane
   that bounds the bed axially. No inscribed staircase, no annular voids,
   no per-slab discretisation.

2. **State-dependent bed top.** Given `V_solid = charge_mass_g / ρ_particle`
   and the target bulk density `pf_state`, the required bed volume is
   `V_bulk = V_solid / pf_state`. When `V_bulk < V_frustum`, the code
   bisects for the height `z_bed_top ∈ (0, z_cone_top]` at which the
   below-plane frustum volume equals `V_bulk`. The bed cell is truncated
   there and the tile is packed at exactly `pf_state`. When
   `V_bulk ≥ V_frustum`, `z_bed_top = z_cone_top` and the excess spills
   into a retort-cylinder overflow cell handled in the same function.
   Fissile inventory `n = charge_mass_g / m_particle` is conserved by
   construction in both branches, and the geometry now differs between
   the two bed states even at fixed charge mass.

3. **Above-bed void inside the cone.** When `z_bed_top < z_cone_top`
   (nominal collapsed case at 95 g), a single `background_material`
   cell fills the empty upper frustum (`-cone_surf & +zp_bed_top &
   -zp_cone_top_plane`). If `z_flood` falls inside this region it is
   split at that plane so the wet/dry boundary is exact.

4. **Tiled fill.** `tiled_bed()` packs one ~0.5 cm cubic tile once (cached,
   keyed on stage / pf / seed), wraps it in a `Universe`, and every position
   of an outer `RectLattice` covering `(-r_retort, -r_retort, 0) →
   (r_retort, r_retort, z_bed_top)` holds *the same universe object*. Unique
   cell count is set by particles per tile (~1600 at bare_kernel, pf ≈ 0.5),
   not by bed volume. The lattice bounding box shrinks with `z_bed_top`,
   so the collapsed case runs slightly fewer lattice traversals than the
   fluidized case.

5. **Boundary clipping.** TRISO particles whose centres sit inside the cone
   but whose outer shells poke past `cone_surf` (or `zp_bed_top`) are
   geometrically clipped by the cell region. The inside portion remains a
   TRISO sub-cell; the outside portion is replaced by whatever material the
   surrounding cell holds (graphite retort wall on the sides,
   `background_material` on top). Affected fraction < 1 % of fissile mass;
   measured effect on k-eff is at MC-noise level (|Δk|/σ_seed ≈ 1.7 at 5 g
   in the tiled-vs-random validation).

6. **No transport outside the packed region.** The cell region is the outer
   boundary. When a neutron crosses `cone_surf` or `zp_bed_top` going
   outward it exits `frustum_cell` entirely and enters whichever cell
   surrounds it — the tile universe is never queried at positions outside
   the truncated frustum. The `RectLattice.outer` fallback exists as a
   safety net but is unreachable from inside the bed cell.

## Design decisions

- **Exact frustum, not inscribed staircase.** Using `ZCone` as the cell
  boundary removes the +2.3–3.3 × 10⁻³ k-eff bias that the staircase
  discretisation carried (see convergence study below). Enabled by the tiled
  `RectLattice`; the direct-packing alternative (one `openmc.Cell` per
  particle) OOMs at 14 GiB Batch containers above ~200 k particles.

- **State-dependent bed top (supersedes mass-conserving tile pf).** Packing
  the entire frustum at `pf_state` over-produces particles when `V_bulk <
  V_frustum`; the previous fix (commit `eec253d`) scaled `pf_tile` by
  `target_n × v_sphere / v_frustum`, which conserved mass but made the tile
  density depend only on `charge_mass_g`, not on `pf_state`. Collapsed and
  fluidized at 95 g therefore produced *identical* geometry and identical
  k-eff. Fix: instead of scaling the density inside a fixed region, bisect
  for the axial plane `z_bed_top` at which the frustum volume equals
  `V_bulk = V_solid / pf_state`, truncate the cell region there, and pack
  the tile at `pf_state` directly. Mass is conserved by construction in
  both the truncated-cone branch and the overflow branch, and the two bed
  states now differ in bed height as well as pf. Verified 2026-09-16: at 95 g
  bare_kernel, collapsed z_bed_top = 3.212 cm (pf = 0.500), fluidized
  z_bed_top = 3.751 cm (pf = 0.333), both with the same 224 960 particles;
  smoke-test transport gives Δk ≈ +4.4 × 10⁻³ (collapsed − fluidized), a
  physically-sensible ~7σ signal in the correct direction for an
  undermoderated bed.

- **Overflow into the retort cylinder.** When `V_bulk > V_frustum` a second
  tile lattice fills a retort-cylinder cell from `z_cone_top` up to the
  height that accommodates the excess. Same tile universe reused, so memory
  stays bounded across the mass sweep.

- **Fluidized bed occupies cone (to z_bed_top) + cylinder overflow when
  needed, not cylinder only.** Uniform fluidization distributes particles
  throughout the available volume at the fluidized packing fraction.
  Placing all particles above the cone top would model a spouted-bed
  regime, not the operating regime of this furnace.

- **Collapsed bed leaves the upper cone empty.** When the packed volume is
  smaller than the frustum, the geometry now respects that: the truncated
  cell stops at `z_bed_top` and the volume above it is filled with
  background gas (or split at `z_flood`). Modelling collapsed as
  "still-fills-the-whole-cone-just-at-lower-effective-density" was
  physically incorrect — a settled bed has a real, finite top surface,
  not a diffuse distribution up to the cone lip.

- **pf_fluidized derived from pf_collapsed / bed_expansion_ratio.** Setting
  the two packing fractions independently would let them drift out of
  consistency with the expansion ratio; deriving pf_fluidized inside the bed
  builder from the ratio ensures physical consistency.

- **background_material is a required parameter.** The choice of background
  (process_gas / water / air) is a model-assembly decision that changes
  between scenarios. Making it explicit forces the caller to state which
  flooding/operating condition is being modelled and prevents a silent
  default that could produce unconservative results.

- **charge_mass_g as an optional kwarg (default from params).** Allows
  parameter sweeps over charge mass without mutating the frozen params
  mapping. Canonical value (95 g) lives in `params.yaml`.

- **Nominal 95 g charge stays within the cone at both packing states, but
  fills a different fraction of it.** At PF = 0.50 (collapsed), bare
  kernels give V_bulk = 18.10 cm³ ≈ 64 % of V_frustum (28.29 cm³) → bed top
  at z ≈ 3.21 cm (of 3.80 cm). At PF = 0.333 (fluidized), V_bulk = 27.14
  cm³ ≈ 96 % of V_frustum → bed top at z ≈ 3.75 cm. Overflow into the
  retort cylinder is exercised by the step-7 mass sweep (up to 3610 g).

- **Convergence study uses vacuum at retort inner cylinder.** Including the
  graphite retort wall would raise absolute k-eff via reflection but not
  change the *relative* bias between geometries. Vacuum gives a more
  conservative absolute value and avoids modelling the wall thickness in
  step 3.

## Convergence study — establishing the choice

**Purpose.** Quantify the k-eff bias of the staircase discretisation versus a
true-frustum reference, using real packed geometries (not homogenised
surrogates), and decide whether any finite `n_slabs` is acceptable for
production.

**Design.** 5 seeds × {exact-cone, n_slabs = 32, 16, 8}, collapsed,
bare_kernel, 95 g charge, full-flood water (z_flood = 999 cm),
20 000 particles/gen × 300 batches. 20 Batch jobs total.

**Sweeps.** `step3_convergence` (staircase, seeds 42–46, git `aab5e58`) and
`step3_convergence_exactcone_tiled` (exact-cone tiled, seeds 42–46, git
`eec253d`).

**Note.** The exact-cone results below were produced with the pre-2026-09-16
mass-conserving-pf implementation, which filled the full frustum at a
reduced tile density. Under the current state-dependent-bed-top geometry
the collapsed bed occupies only the bottom ~64 % of the frustum, so the
absolute exact_cone reference k-eff at 95 g full-flood collapsed will shift
(the fissile mass is the same but the wet-cone surface area and top-of-bed
leakage geometry change). The direction and magnitude of the
staircase-vs-exact bias is a property of the discretisation, not the
reference geometry, so the +2–3 × 10⁻³ finding still holds qualitatively;
the absolute numbers should be re-measured after the geometry change if
step 3 is reopened.

**Results.**

| config     | mean k-eff | across-seed σ | mean MC σ | pf achieved | pf deviation | bias vs exact-cone | \|Δk\|/σ_seed |
|------------|------------|---------------|-----------|-------------|--------------|--------------------|---------------|
| exact_cone | 0.024170   | 3×10⁻⁵        | 5×10⁻⁵    | 0.4997      | −0.0003      | reference          | —             |
| n=32       | 0.026475   | 9.9×10⁻⁵      | 4×10⁻⁵    | 0.4786      | −0.0214      | +2.31×10⁻³         | 22.4          |
| n=16       | 0.027116   | 4.9×10⁻⁵      | 4×10⁻⁵    | 0.4978      | −0.0022      | +2.95×10⁻³         | 51.5          |
| n=8        | 0.027453   | 9.3×10⁻⁵      | 5×10⁻⁵    | 0.5000      | −0.0000      | +3.28×10⁻³         | 33.8          |

Data: `results/step3_convergence/1789411575_c1f98c934540/summary.csv`,
`results/step3_convergence_exactcone_tiled/1789510193_1b02b68a59d0/summary.csv`.

**Findings.**

- **Every staircase config fails the δk_disc < σ_seed criterion.**
  |Δk|/σ_seed ranges from 22 to 52. n_slabs = 32 (the previous production
  choice) does not converge with real packed particles.
- **Bias direction is correct.** Staircase overestimates k-eff under full
  flood, matching the a-priori prediction (cone→cylinder relocation reduces
  moderation efficiency and raises leakage).
- **Bias is non-monotone in n.** n=32 has slightly smaller raw bias than
  n=16 because its pf deficit (0.479 vs 0.500) removes fissile mass — a
  packer artefact, not convergence.
- **Slab-thickness floor caps usable n.** The `pack_bed` RSP guard requires
  slab_h ≥ 5 × particle diameter. At bare_kernel (d ≈ 0.0425 cm) this gives
  n_slabs ≤ ~18 for this cone; coated stages tighten the ceiling further.
  n=32 (2.8 d) violates the floor, which is what produces the pf deficit and
  the seed-to-seed u235 mass variance (17.03–17.19 g at n=32 vs uniform for
  n ≤ 16). **No finite n_slabs is both accurate and safe across the sweep
  matrix.**
- **Only new artefact from the tiled path is periodic tile structure** at
  0.5 cm scale — measured contribution to k-eff is ~10⁻⁵ (MC-noise level),
  one to two orders of magnitude below the discretisation bias it replaces.

**Conclusion.** Exact-cone tiled adopted as production geometry. Staircase
retired.

## Historical: staircase approach (retired)

Before 2026-09-16 the bed was modelled as `n_slabs` equal-height inscribed
cylinders stacked from the cone base, each with radius equal to the cone
surface at its bottom face. Both bed states filled slabs bottom-up at the
state-appropriate packing fraction, with overflow into the retort cylinder
above.

**Why it existed.** Random `pack_bed` placed each TRISO as a distinct
`openmc.Cell`. At 95 g bare_kernel that is ~225 k unique cells; OpenMC's
"Preparing distributed cell instances" phase materialises each with a
per-instance transformation table, and the working set exceeded the 14 GiB
Batch cap at cell-instance preparation time. Slicing the cone into slabs kept
any single lattice small enough to survive. `n_slabs = 32` was selected by an
earlier homogenised study (below) as the finest safe discretisation.

**Why it was retired.** The tiled-bed refactor (`8478711`, `76f0e5c`) reduced
unique-cell count from O(N_particles) to O(particles_per_tile), removing the
memory constraint that motivated the staircase. Follow-up fixes `eec253d`
(pf-scaling so the tile respects `charge_mass_g`) and `1dbb223` (skip dead
`pack_spheres` in the tiled path) made whole-cone packing cheap at any charge
mass. The packed-particle convergence study above then showed the staircase
carries a +2–3 × 10⁻³ k-eff bias that no finite n_slabs eliminates while also
respecting the packer's slab-thickness floor, so the staircase became
strictly worse than the direct exact-cone geometry.

**Legacy design decisions (still apply to `bed_region` calls if run for
reproducibility):**
- Inscribed staircase (radius at slab bottom), not midpoint or circumscribed
  — required to prevent overlapping-cell geometry with the graphite cone
  wall.
- Each staircase slab uses a distinct seed (`base_seed + slab_index`) to
  avoid spurious spatial periodicity.
- Last partial slab at n=32 / 95 g / bare_kernel is near-monolayer
  (h ≈ 0.064 cm ≈ 1.5 d), trips the RSP guard in `pack_bed`, and produces
  the pf deficit visible in the convergence study.

**Superseded homogenised study.** Kept for context on why `n_slabs = 32` was
the pre-2026 default. Uses homogenised (smeared) bed materials instead of
explicit TRISO packing; conflates geometry and packing effects and cannot
quantify the real-geometry bias.

| n_slabs | V_error | ΔV/V_bulk | k_eff (homogenised) | σ | δk_disc |
|---------|---------|-----------|---------------------|---|---------|
| 4  | 30.9% | 48.4% | 0.02638 | 4×10⁻⁵ | 2.87×10⁻⁴ |
| 8  | 16.4% | 25.6% | 0.02609 | 4×10⁻⁵ | 1.56×10⁻⁴ |
| 16 |  8.7% | 13.6% | 0.02594 | 4×10⁻⁵ | 5.32×10⁻⁵ |
| 32 |  4.7% |  7.4% | 0.02588 | 4×10⁻⁵ | 3.83×10⁻⁵ |
| 64 |  2.7% |  4.3% | 0.02585 | 4×10⁻⁵ | — |

## Assumptions

### Confirmed
- Coordinate system: z = 0 at nozzle-to-cone junction; z increases upward.
- Cone half-angle = 30° (from `included_angle_deg = 60°`) → r² = tan²(30°) = 1/3.
- Nozzle throat: graphite only; particles never enter it.
- Fluidized bed occupies the cone (at pf_fluidized) plus overflow into the
  retort cylinder above (uniform fluidization, not spouted-bed).
- Background material is a caller-supplied parameter (not a params default).
- Convergence criterion: bias vs exact-cone reference < across-seed σ.
  **Result of packed-particle study: no staircase configuration meets this
  criterion; exact-cone tiled adopted as production.**

### Defaulted
- **Throat height 30.0 mm**: no fabrication drawing; a reasonable value for
  a ø6 mm orifice throat. Marked `# CONFIRM`.
- **bed_expansion_ratio = 1.5**: confirmed against process measurements for
  near-minimum-fluidization operation; yields pf_fluidized = 0.50 / 1.5 =
  0.333. Marked `# CONFIRM`.
- **charge_mass_g = 95.0 g**: representative starting charge; matches
  approximate kernel count from step 2. Marked `# CONFIRM`.
- **static_bulk_cc = 18.10 cm³**: derived from
  charge_mass_g / rho_eff / pf_static = 95 / 10.5 / 0.50.

### Unconfirmed (`# CONFIRM`)
- `nozzle_height: 30.0` mm — no fabrication drawing; placeholder.
- `bed_expansion_ratio: 1.5` — confirmed against process data for this
  operating condition; should be verified against full-range fluidization
  measurements.
- `charge_mass_g: 95.0` g — representative starting charge; must be confirmed
  against process procedure.
- All `# CONFIRM` flags from steps 1 and 2 still apply: kernel diameter, all
  layer thicknesses, structural graphite density.
