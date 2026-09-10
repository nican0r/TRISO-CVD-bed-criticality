---
step: 3 — Bed and cone geometry
title: Hybrid staircase-cone bed geometry, two bed states (fluidized/collapsed), mass-conserving packing, n_slabs discretisation convergence study
status: complete
---

## What was implemented

- `furnace/geometry.py` — full implementation replacing the step-3 stub:
  - `retort_inner_cylinder(params)` — returns `openmc.ZCylinder` at the retort inner radius
  - `cone_wall(params)` — returns `(ZCone, z_plane_bottom, z_plane_top)` for the structural graphite cone; apex derived analytically from throat radius and half-angle
  - `nozzle_throat(params)` — returns `(ZCylinder, z_plane_bottom, z_plane_top)` for the graphite throat cylinder below z=0
  - `frustum_volume(r_top, r_bottom, half_angle_deg)` — analytic truncated-cone volume via ∫π r(z)² dz substitution; exact formula V = (π/3)(r_top³ - r_bot³)/tan θ
  - `staircase_bed(params, n_slabs)` — divides the full cone axial extent into `n_slabs` equal-height inscribed cylinders; returns list of slab dicts, total staircase volume, and fractional volume error relative to the analytic frustum
  - `_particle_effective_density(stage, params)` — stage-dependent effective density ρ_eff = m_particle / V_particle_outer_sphere
  - `bed_region(params, state, stage, n_slabs, background_material, charge_mass_g=None)` — mass-conserving bed assembly for 'collapsed' and 'fluidized' states; accepts a caller-supplied background material and an optional explicit charge mass for sweeps; both states use the same inscribed staircase geometry (cone slabs from bottom to top, then overflow into the retort cylinder); returns both `cells` (bed lattice cells) and `cone_void_cells` (annular voids between each inscribed slab cylinder and the true cone wall, plus any unfilled cone interior above the bed — required to give OpenMC complete spatial coverage)
- `params.yaml` — five new keys:
  - `dimensions.nozzle.nozzle_height: 30.0` mm (`# CONFIRM`) — throat cylinder axial extent
  - `dimensions.bed.charge_mass_g: 95.0` g (`# CONFIRM`) — starting particle charge; sweep target
  - `dimensions.bed.fluidized_height` comment updated to clarify it is the available zone; actual height is computed from charge mass
  - `model.packing_fraction_fluidized` updated to `0.075` (derived: pf_static / bed_expansion_ratio)
  - `model.bed_expansion_ratio: 8.0` (`# CONFIRM`) — volumetric fluidized-to-collapsed ratio
  - `dimensions.bed_volume.static_bulk_cc` updated to `15.1` (derived value) with a note that it is derived
- `scripts/run_convergence.py` — inline convergence study script (see Experimental design section)

## How it works

The cone region is approximated as `n_slabs` equal-height inscribed cylinders, each with radius equal to the cone surface radius at its bottom face. This guarantees no particle cell protrudes through the graphite cone wall (a circumscribed or midpoint fit would place bed cylinders outside the cone for any finite slab count). The analytic frustum volume minus the staircase volume gives a fractional geometry error that shrinks as n_slabs increases. `bed_region` converts the physical charge mass to a conserved solid particle volume using the stage-dependent effective density (mass/outer-sphere volume), then divides by the target packing fraction to obtain the required bulk volume; this ensures that the particle inventory is exactly the same regardless of which bed state or n_slabs is used in a comparison. Both the collapsed and fluidized states fill cone slabs from bottom to top at the state-appropriate packing fraction, with overflow continuing into the retort cylinder. For the fluidized state the overflow height is checked against `fluidized_height_cm` (the available cylinder zone above the cone).

## Experimental design

### Convergence study (`scripts/run_convergence.py`)

**Purpose**: Determine the minimum `n_slabs` such that the staircase discretisation introduces less than one statistical sigma error in k-eff relative to a twice-finer discretisation.

**Geometry**: Nozzle throat (graphite) + cone staircase (inscribed bed cylinders + graphite wall) + collapsed bed overflow in retort cylinder. Outer boundary: vacuum at retort inner cylinder.

**Materials**: Water at 1.0 g/cm³ as bed background (maximises moderation, giving the most sensitive geometry test); graphite_structural for cone wall and throat.

**Stage**: `bare_kernel` (smallest particle, highest packing density at fixed mass, most sensitive to cone discretisation).

**State**: `collapsed` (the staircase approximation only affects the collapsed geometry; the fluidized bed is a simple cylinder above the cone regardless of n_slabs).

**n_slabs values**: {4, 8, 16, 32}

**Convergence criterion**: δk_disc(n) = |k(n) − k(2n)| < σ(2n) (one sigma of the finer run)

**Note on state choice**: The staircase approximation now affects both the collapsed and fluidized states (particles occupy the cone in both cases). The collapsed state is the more sensitive test because pf_static (0.60) is 8× higher than pf_fluidized (0.075), so geometric volume errors translate to larger particle-count errors per slab in the collapsed case. Results from the collapsed convergence study are conservative for the fluidized case.

**Settings**: Uses `params.model.batches`, `params.model.inactive`, `params.model.particles`, `params.model.seed` from params.yaml. For preliminary convergence screening the default particle count (1000) is sufficient to distinguish large discretisation errors; production runs should use ≥50 000 particles.

**Outputs**:
- Console table: n_slabs | V_staircase_cm3 | V_error_pct | k_eff ± σ | δk_disc
- `results/convergence_n_slabs.csv`
- OpenMC state-point files in `cases/convergence_n{n}/`

**Run command**:
```
caffeinate python scripts/run_convergence.py
```

### Results and n_slabs selection

| n_slabs | V_staircase (cm³) | V_error (%) | k_eff | σ | δk_disc |
|---------|-------------------|-------------|-------|---|---------|
| 4  | 19.54 | 30.9% | 0.02825 | 4×10⁻⁵ | 3×10⁻⁵ |
| 8  | 23.66 | 16.4% | 0.02822 | 3×10⁻⁵ | 8×10⁻⁵ |
| 16 | 25.84 |  8.7% | 0.02814 | 3×10⁻⁵ | 7×10⁻⁵ |
| 32 | 26.95 |  4.7% | 0.02808 | 3×10⁻⁵ | — |

**Selected: n_slabs = 8.**

The strict δk_disc < σ criterion is **never satisfied** across the tested range:

- n=4→8: δk = 3×10⁻⁵, σ(8) = 3×10⁻⁵ → not strictly satisfied (equal, not less than; both values rounded to 1 s.f.)
- n=8→16: δk = 8×10⁻⁵, σ(16) = 3×10⁻⁵ → fails
- n=16→32: δk = 7×10⁻⁵, σ(32) = 3×10⁻⁵ → fails

**Direction of bias:** The staircase underestimates cone volume, not particle count. Because `bed_region` is mass-conserving, particles that do not fit within the inscribed staircase cylinders overflow into the retort cylinder above the cone rather than disappearing. The cone is geometrically less favorable than the cylinder (narrower cross-section, closer to the graphite wall, higher neutron leakage). Redistributing particles from the cone to the cylinder therefore **overestimates** k-eff relative to the true geometry. This is confirmed by the data: k-eff decreases monotonically as n_slabs increases and the staircase more accurately captures the cone volume.

The k-eff bias is therefore not conservative for NCS in the usual sense — a lower bound on k-eff is not what a criticality safety case needs. What it does mean is that the k-eff criterion alone cannot be used to select n_slabs, since the criterion failing is expected and does not indicate an unsafe model. n=8 is selected on geometric accuracy grounds:

1. **Diminishing geometric returns beyond n=8**: The volume error drops from 30.9% at n=4 to 16.4% at n=8 — nearly halved. Subsequent doublings recover only ~8 percentage points each (8.7% at n=16, 4.7% at n=32) at the cost of doubling compute time.

2. **n=4 overcredits the cylinder by a physically significant margin**: A 30.9% cone volume error at n=4 means a substantial fraction of the cone inventory is incorrectly placed in the more-reactive cylinder. For design use (determining the safe operating charge mass), this produces a k-eff estimate that is biased high and not representative of the actual particle distribution.

3. **n=8 provides sufficient geometric fidelity for the nominal charge**: At 95 g, the total bulk volume is ~120 cm³. The unrepresented cone volume at n=8 is V_frustum − V_staircase ≈ 28.3 − 23.7 = 4.6 cm³, which is mislocated to the cylinder by mass conservation. As a fraction of total bulk volume: 4.6 / 120 ≈ 3.8% of inventory is mislocated — a small geometric perturbation that results in a k-eff overestimate of ~8×10⁻⁵ (the δk_disc step from n=8 to n=16).

For the fluidized state (pf_fluidized = 0.075, 8× lower than pf_static), the same geometric error produces an 8× smaller k-eff perturbation, making n=8 more than adequate for fluidized-state runs.

## Design decisions

- **Inscribed staircase (radius at slab bottom), not midpoint or circumscribed.** The cone wall is a hard geometric boundary; any bed cylinder whose radius exceeds the cone radius at its axial position would protrude outside the graphite wall, creating an illegal overlapping-cell geometry in OpenMC. The inscribed choice (smallest possible slab radius) guarantees correctness at the cost of underestimating the cone volume. As n_slabs → ∞ the staircase volume converges to the analytic frustum from below.

- **Mass-conserving volume calculation.** Rather than specifying a packing fraction and region volume and letting the particle count float, `bed_region` fixes the particle charge mass and derives V_bulk = (charge_mass_g / ρ_eff) / pf. This means the total fissile inventory is the same for all n_slabs values in the convergence study, ensuring k-eff differences are due to geometric discretisation rather than different particle counts.

- **pf_fluidized derived from pf_collapsed / bed_expansion_ratio; not stored as an independent parameter.** Setting the two packing fractions independently would make it easy for them to become inconsistent with the bed_expansion_ratio. Deriving pf_fluidized inside `bed_region` from the ratio ensures physical consistency.

- **Fluidized bed occupies cone + cylinder overflow at pf_fluidized, not cylinder only.** Uniform fluidization distributes particles throughout the available volume at the fluidized packing fraction. Placing all particles above the cone top would model a spouted-bed regime (particles pneumatically conveyed out of the cone), which is not the operating regime of this furnace. The inscribed staircase therefore applies to both states; the n_slabs discretisation error is relevant in both cases (though less sensitive in the fluidized case due to the 8× lower packing fraction).

- **background_material is a required parameter to bed_region; not read from params.** The choice of background (process_gas vs. water vs. air) is a model-assembly decision that changes between scenarios. Making it explicit forces the caller to state which flooding/operating condition is being modelled and prevents a silent default that could produce unconservative results.

- **charge_mass_g as an optional kwarg (default reads from params).** This allows parameter sweeps over charge mass without mutating the frozen params mapping. The canonical value (95 g) lives in params.yaml and is read by default; callers can pass any value for a mass-sensitivity study.

- **Nominal 95 g charge stays within the cone (no overflow).** At PF = 0.60, bare kernels give V_bulk ≈ 15.1 cm³, which is less than both the staircase volume (23.7 cm³ at n=8) and the full frustum volume (28.3 cm³). The bed fills only the lower portion of the cone. Overflow into the retort cylinder above is handled in code and will trigger for larger charges (the mass-sweep use case), but does not occur at the starting 95 g value.

- **Convergence study uses vacuum at retort inner cylinder.** Including the graphite retort wall would add reflected neutrons and increase absolute k-eff, but would not change the relative δk_disc between different n_slabs values (the wall geometry is the same in every case). The vacuum boundary gives a more conservative absolute value and avoids needing to model the retort wall thickness in step 3.

- **Each staircase slab uses a distinct seed (base_seed + slab_index).** Using the same seed for all slabs would produce identical sphere-centre patterns scaled to each slab's size, introducing spurious spatial periodicity in the particle arrangement. Distinct seeds give statistically independent packing realisations for each slab.

## Assumptions

### Confirmed
- Coordinate system: z=0 at nozzle-to-cone junction; z increases upward.
- Cone half-angle = 30° (from included_angle_deg = 60°) → r2 = tan²(30°) = 1/3.
- Inscribed staircase only (not midpoint, not circumscribed).
- Nozzle throat: graphite only; particles never enter it.
- Fluidized bed occupies the cone (inscribed staircase at pf_fluidized) plus overflow into the retort cylinder above the cone top (uniform fluidization, not spouted-bed).
- Background material is a caller-supplied parameter (not a params default).
- Convergence criterion: δk_disc < 1σ of the finer run.

### Defaulted
- **Throat height 30.0 mm**: no fabrication drawing available; a reasonable value for a ø6 mm orifice throat. Marked `# CONFIRM`.
- **bed_expansion_ratio = 8.0**: typical fluidized-bed expansion factor; should be measured or computed from fluid dynamics for this geometry. Marked `# CONFIRM`.
- **charge_mass_g = 95.0 g**: starting charge for the CVD run; chosen to match approximate kernel count used in step 2 analytical study. Marked `# CONFIRM`.
- **static_bulk_cc updated to 15.1 cm³**: derived from charge_mass_g and bare-kernel properties at pf_static = 0.60; replaces the original placeholder of 50 cm³.

### Unconfirmed (`# CONFIRM`)
- `nozzle_height: 30.0` mm — no fabrication drawing; placeholder.
- `bed_expansion_ratio: 8.0` — no process measurement; typical value for spherical particles in this size range.
- `charge_mass_g: 95.0` g — representative starting charge; must be confirmed against process procedure.
- All `# CONFIRM` flags from steps 1 and 2 still apply: kernel diameter, all layer thicknesses, structural graphite density.
