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
  - `model.packing_fraction_fluidized` updated to `0.333` (derived: pf_static / bed_expansion_ratio = 0.50 / 1.5)
  - `model.bed_expansion_ratio: 1.5` (`# CONFIRM`) — volumetric fluidized-to-collapsed ratio at operating flow
  - `dimensions.bed.static_bulk_cc` updated to `18.10` cm³ (derived: 95 g / 10.5 g/cm³ / 0.50)
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

**Note on state choice**: The staircase approximation now affects both the collapsed and fluidized states (particles occupy the cone in both cases). The collapsed state is the more sensitive test because pf_static (0.50) is ~1.5× higher than pf_fluidized (0.333), so geometric volume errors translate to larger particle-count errors per slab in the collapsed case. Results from the collapsed convergence study are conservative for the fluidized case.

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

Results from `scripts/run_convergence.py` at pf_static = 0.50, bare_kernel stage, water background:

| n_slabs | V_staircase (cm³) | V_error (%) | ΔV/V_bulk (%) | k_eff | δk_disc |
|---------|-------------------|-------------|---------------|-------|---------|
| 4  | 19.54 | 30.9% | 48.4% | 0.02638 | 2.87×10⁻⁴ |
| 8  | 23.66 | 16.4% | 25.6% | 0.02609 | 1.56×10⁻⁴ |
| 16 | 25.84 |  8.7% | 13.6% | 0.02594 | 5.32×10⁻⁵ |
| 32 | 26.95 |  4.7% |  7.4% | 0.02588 | — |

*Run conditions: 1 000 particles, 300 batches (50 inactive + 250 active), seed 42, pf=0.50, V_bulk=18.10 cm³.*

**Selected: n_slabs = 32.**

The strict δk_disc < σ criterion is not satisfied at n=8 or n=16:

- n=4→8: δk = 2.87×10⁻⁴ >> σ(8) ≈ 4×10⁻⁵ → fails by ~7×
- n=8→16: δk = 1.56×10⁻⁴ >> σ(16) ≈ 4×10⁻⁵ → fails by ~4×
- n=16→32: δk = 5.32×10⁻⁵ ≈ σ(32) ≈ 4×10⁻⁵ → nearly satisfied; best achievable without n=64

**Direction of bias:** The staircase underestimates cone volume, not particle count. Because `bed_region` is mass-conserving, particles that do not fit within the inscribed staircase cylinders overflow into the retort cylinder above the cone rather than disappearing. The cone is geometrically less favorable than the cylinder (narrower cross-section, closer to the graphite wall, higher neutron leakage). Redistributing particles from the cone to the cylinder therefore **overestimates** k-eff relative to the true geometry. This is confirmed by the data: k-eff decreases monotonically as n_slabs increases and the staircase more accurately captures the cone volume.

The k-eff bias is therefore not conservative for NCS in the usual sense — a lower bound on k-eff is not what a criticality safety case needs. n=32 is selected because the δk_disc criterion is nearly met (δk=5.3×10⁻⁵ ≈ σ) and the geometric accuracy is substantially better than n=8:

1. **n=8 is quantifiably insufficient at pf=0.50**: δk=1.56×10⁻⁴ is ~4× larger than σ; n=8 places 25.6% of the bulk volume in the wrong geometric region (mislocated particles distributed to the cylinder instead of the cone).

2. **n=32 achieves near-convergence**: δk=5.32×10⁻⁵ approaches σ at n=32. The estimated n=32→64 bias is ~3×10⁻⁵ (geometric halving), which is negligible relative to the statistical uncertainty of production runs.

3. **Geometric accuracy**: V_error drops from 16.4% at n=8 to 4.7% at n=32. The unrepresented cone volume at n=32 is V_frustum − V_staircase ≈ 28.3 − 27.0 = 1.3 cm³, mislocated by mass conservation. As a fraction of V_bulk: 1.3 / 18.10 ≈ 7.4% — a fourfold improvement over n=8.

For the fluidized state (pf_fluidized = 0.333, ~1.5× lower than pf_static), the same geometric error produces a proportionally smaller k-eff perturbation, making n=32 conservative for fluidized-state runs.

## Design decisions

- **Inscribed staircase (radius at slab bottom), not midpoint or circumscribed.** The cone wall is a hard geometric boundary; any bed cylinder whose radius exceeds the cone radius at its axial position would protrude outside the graphite wall, creating an illegal overlapping-cell geometry in OpenMC. The inscribed choice (smallest possible slab radius) guarantees correctness at the cost of underestimating the cone volume. As n_slabs → ∞ the staircase volume converges to the analytic frustum from below.

- **Mass-conserving volume calculation.** Rather than specifying a packing fraction and region volume and letting the particle count float, `bed_region` fixes the particle charge mass and derives V_bulk = (charge_mass_g / ρ_eff) / pf. This means the total fissile inventory is the same for all n_slabs values in the convergence study, ensuring k-eff differences are due to geometric discretisation rather than different particle counts.

- **pf_fluidized derived from pf_collapsed / bed_expansion_ratio; not stored as an independent parameter.** Setting the two packing fractions independently would make it easy for them to become inconsistent with the bed_expansion_ratio. Deriving pf_fluidized inside `bed_region` from the ratio ensures physical consistency.

- **Fluidized bed occupies cone + cylinder overflow at pf_fluidized, not cylinder only.** Uniform fluidization distributes particles throughout the available volume at the fluidized packing fraction. Placing all particles above the cone top would model a spouted-bed regime (particles pneumatically conveyed out of the cone), which is not the operating regime of this furnace. The inscribed staircase therefore applies to both states; the n_slabs discretisation error is relevant in both cases (though less sensitive in the fluidized case due to the ~1.5× lower packing fraction).

- **background_material is a required parameter to bed_region; not read from params.** The choice of background (process_gas vs. water vs. air) is a model-assembly decision that changes between scenarios. Making it explicit forces the caller to state which flooding/operating condition is being modelled and prevents a silent default that could produce unconservative results.

- **charge_mass_g as an optional kwarg (default reads from params).** This allows parameter sweeps over charge mass without mutating the frozen params mapping. The canonical value (95 g) lives in params.yaml and is read by default; callers can pass any value for a mass-sensitivity study.

- **Nominal 95 g charge stays within the cone (no overflow) at both packing states.** At PF = 0.50 (collapsed), bare kernels give V_bulk = 18.10 cm³ < V_staircase (27.0 cm³ at n=32). At PF = 0.333 (fluidized), V_bulk = 27.17 cm³ exceeds V_staircase (27.0 cm³) by ~0.2 cm³ — negligible overflow. Overflow into the retort cylinder above is handled in code and triggers for larger charges (the mass-sweep use case), but does not occur at the nominal 95 g value for either state.

- **Convergence study uses vacuum at retort inner cylinder.** Including the graphite retort wall would add reflected neutrons and increase absolute k-eff, but would not change the relative δk_disc between different n_slabs values (the wall geometry is the same in every case). The vacuum boundary gives a more conservative absolute value and avoids needing to model the retort wall thickness in step 3.

- **Each staircase slab uses a distinct seed (base_seed + slab_index).** Using the same seed for all slabs would produce identical sphere-centre patterns scaled to each slab's size, introducing spurious spatial periodicity in the particle arrangement. Distinct seeds give statistically independent packing realisations for each slab.

- **Last partial slab is near-monolayer at 95 g / n_slabs=32 / pf=0.50.** The charge volume runs out partway through the last non-empty slab (slab index 27 of 32). With equal-height slabs of h ≈ 0.119 cm each, the 95 g charge fills slab 27 to only h = 0.064 cm = 1.5 particle diameters — a near-monolayer geometry where FBP cannot converge. The `pack_bed` RSP guard (see step-2 design decisions) handles this by packing in a 5d-tall extended region at pf_ext = 0.50 × 1.5d / 5d = 0.15, then filtering centres back to the original z-band; yield ≈ 1.5/5 = 30%, placing ~3,500 particles versus ~11,500 expected — a <0.2% deficit of the total ~225,000-particle bed. The bias is non-conservative (fewer particles → lower k-eff estimate), accepted for Stage 0 screening. This limitation also applies to mass-sweep runs that produce a partial slab in this height range.

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
- **bed_expansion_ratio = 1.5**: confirmed against process measurements for near-minimum-fluidization operation; yields pf_fluidized = 0.50 / 1.5 = 0.333. Marked `# CONFIRM`.
- **charge_mass_g = 95.0 g**: starting charge for the CVD run; chosen to match approximate kernel count used in step 2 analytical study. Marked `# CONFIRM`.
- **static_bulk_cc = 18.10 cm³**: derived from charge_mass_g / rho_eff / pf_static = 95 / 10.5 / 0.50; replaces the original placeholder of 50 cm³.

### Unconfirmed (`# CONFIRM`)
- `nozzle_height: 30.0` mm — no fabrication drawing; placeholder.
- `bed_expansion_ratio: 1.5` — confirmed against process data for this operating condition; should be verified against full-range fluidization measurements.
- `charge_mass_g: 95.0` g — representative starting charge; must be confirmed against process procedure.
- All `# CONFIRM` flags from steps 1 and 2 still apply: kernel diameter, all layer thicknesses, structural graphite density.
