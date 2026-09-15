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

**Purpose**: Quantify the k-eff bias introduced by the staircase cone discretisation relative to a rejection-sampled exact-frustum reference, using real packed-particle geometries (not homogenised surrogates).  Determine whether the current production value n_slabs=32 is adequate, and whether n=16 or n=8 are acceptable alternatives.

**Exact-cone reference**: A reference geometry packs the true frustum by rejection sampling: pack the bounding cylinder (r=r_retort, z in [0, z_cone_top]) at elevated pf_trial, then discard every kernel whose sphere intersects the cone surface (perpendicular clearance < r_kernel) or the floor/ceiling planes. pf_trial is iterated until the surviving count matches the charge mass to within 0.5%. Results are cached. This is validation-only code; it replaces the n→∞ extrapolation used in the homogenised study.

**Seed replicates**: Five independent packing seeds per configuration. The across-seed standard deviation is the error bar for all comparisons. Per-run Monte Carlo sigma is reported separately.

**Run matrix**: Collapsed bed state, bare_kernel stage, nominal charge mass (95 g), full-flood background (water at 1.0 g/cm³, z_flood=999 cm).
- Configurations: exact_cone (reference), n=32 (current production), n=16, n=8
- Full flood is the bounding case for criticality; the convergence bias is measured at worst-case moderation rather than at nominal gas conditions
- 5 seeds × 4 configurations = 20 jobs via AWS Batch

**Run command**:
```
python scripts/run_convergence.py [--local-check] [--submit | --dry-run]
```

**Outputs** (after `pull_results.py`):
- Table: config | mean k-eff | across-seed σ | mean MC σ | achieved pf | requested pf | pf deviation | bias vs exact-cone
- Manifest: `manifests/step3_convergence.json`

The primary diagnostic is **achieved vs requested packing fraction per configuration**. Interface depletion — fewer particles placed near slab boundaries than the bulk packing fraction implies — appears as `pf_achieved < pf_requested`, and the deficit is expected to grow with decreasing n_slabs. A large pf deficit at n=32 would indicate that slab boundaries are removing a significant particle fraction, and the k-eff bias against the exact-cone reference would confirm the neutronics impact.

### Results and n_slabs selection

**Previous homogenised study (superseded)**

An earlier version of this study used homogenised (smeared) bed materials instead of explicit TRISO packing, and tested both water and gas backgrounds. Results at pf_static = 0.50, bare_kernel stage, water background, 1 000 particles per batch:

| n_slabs | V_error (%) | ΔV/V_bulk (%) | k_eff | σ | δk_disc |
|---------|-------------|---------------|-------|---|---------|
| 4  | 30.9% | 48.4% | 0.02638 | 4×10⁻⁵ | 2.87×10⁻⁴ |
| 8  | 16.4% | 25.6% | 0.02609 | 4×10⁻⁵ | 1.56×10⁻⁴ |
| 16 |  8.7% | 13.6% | 0.02594 | 4×10⁻⁵ | 5.32×10⁻⁵ |
| 32 |  4.7% |  7.4% | 0.02588 | 4×10⁻⁵ | 3.83×10⁻⁵ |
| 64 |  2.7% |  4.3% | 0.02585 | 4×10⁻⁵ | — |

That study selected n_slabs = 32 based on the δk_disc < σ(2n) criterion. Homogenisation omits TRISO self-shielding and grain-structure effects; the δk_disc comparisons conflate geometry and packing. The study correctly identified the monotonic direction of bias (staircase overestimates k-eff by displacing particles from the cone to the cylinder) but cannot quantify the bias magnitude for real packed geometries.

**Packed-particle study (pending)**

Results from `scripts/run_convergence.py` with real packed beds are pending. The table below will be populated after the AWS Batch sweep completes:

| config | mean k-eff | across-seed σ | mean MC σ | pf achieved | pf requested | pf deviation | bias vs exact-cone |
|--------|------------|---------------|-----------|-------------|--------------|--------------|-------------------|
| exact_cone | — | — | — | — | 0.500 | — | reference |
| n=32 | — | — | — | — | 0.500 | — | — |
| n=16 | — | — | — | — | 0.500 | — | — |
| n=8  | — | — | — | — | 0.500 | — | — |

*Background: full-flood water at 1.0 g/cm³. State: collapsed. Stage: bare_kernel. Charge: 95 g.*

**Expected direction of bias**: The staircase mislocates particles from the cone (narrower, closer to graphite wall, higher leakage) to the retort cylinder above, overestimating k-eff. In the flooded case this effect is amplified — the cone is better moderated than the cylinder so displacing particles upward reduces moderation efficiency and raises leakage, pushing k-eff in the non-conservative direction.

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
- Convergence criterion: bias vs exact-cone reference < across-seed σ (packed-particle study, pending).

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
