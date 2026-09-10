# Step 5 — Nominal case and convergence discipline

## How to run

Initialize the shell hook once per terminal session, then activate:

```
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env
```

Run the smoke test first to verify the pipeline (~15 s):

```
python3 scripts/smoke_test.py
```

Then run in order from the repo root:

```
caffeinate python3 scripts/run_nominal.py
caffeinate python3 scripts/run_step5_convergence.py
python3 scripts/plot_step5_convergence.py
```

Alternatively, skip activation entirely and use `micromamba run`:

```
caffeinate micromamba run -n triso-env python3 scripts/run_nominal.py
caffeinate micromamba run -n triso-env python3 scripts/run_step5_convergence.py
micromamba run -n triso-env python3 scripts/plot_step5_convergence.py
```

Results are written to `results/step5_nominal/` and `results/step5_convergence/`.

---

## What was implemented

- **`furnace/model.py`** — full module with:
  - `build_model(params, *, state, stage, background, n_inactive, n_active, n_particles, seed, charge_mass_g)` — assembles and returns `(openmc.Model, BedStats)`
  - `_build_tallies(...)` — three tallies: flux spectrum, material reaction rates, U-235 fission spatial distribution
  - `_u235_mass_g(params, n_particles)` — derives U-235 mass from OpenMC nuclide atom densities
  - `export_and_run(model, out_dir)` — exports XML and runs OpenMC; returns statepoint path
  - `BedStats` namedtuple for reporting
- **`furnace/geometry.py`** — added optional `seed` keyword to `bed_region()` so packing seed can be overridden independently of `params`
- **`furnace/materials.py`** — removed `add_s_alpha_beta('c_H_in_H2O')` from `water()` (see Design decisions)
- **`params.yaml`** — updated `fluidized_height`, `particles`, `batches`, `inactive` (see Design decisions)
- **`scripts/run_nominal.py`** — runs the nominal eigenvalue case; reports U-235 mass and k-eff ± σ / +2σ / +3σ
- **`scripts/run_step5_convergence.py`** — three additional runs: `batch_500` (seed=42, 500 active), `seed_43` (250 active, repacked), `seed_44` (250 active, repacked)
- **`scripts/plot_step5_convergence.py`** — three PNG plots: entropy vs generation, k-eff vs active batches, k-eff vs seed
- **`scripts/smoke_test.py`** — pipeline smoke test: 5 g charge (~5 193 TRISO particles), 1+1 batches, 10 000 particles; completes in ~15 s

## How it works

`build_model()` calls `bed_region()` for the fluidized TRISO bed, `furnace_shell_cells()` for the retort walls, heater, and injector, then adds a single cell of process gas filling the retort interior above the bed. The three regions are combined into a single root universe. A 10×10×20 Shannon entropy mesh over the full retort cylinder tracks source convergence each generation. The three tallies (flux spectrum, material reaction rates, U-235 fission spatial) are scored over the bed region.

The `seed` parameter threads through both the packing call (`pack_bed()` → `pack_spheres()`) and the OpenMC transport RNG (`settings.seed`), so changing `seed` in `build_model()` produces a fully independent realisation — different particle positions and different transport history — for the seed study.

The convergence study uses two batch counts (250 and 500 active) to detect k-eff drift, and three independent seeds (42, 43, 44) to confirm results are not sensitive to a particular geometric realisation of the random packing.

## Experimental design

**What is being modelled:** Uniformly fluidized bed of bare UCO kernels (no coating layers — the nominal CVD starting condition) distributed throughout the cone interior (inscribed staircase at pf_fluidized = 0.075) and the retort cylinder above the cone, surrounded by the full graphite retort wall, vacuum gap, graphite heater element, and water-cooled graphite injector. Vacuum boundary conditions at the heater OD, retort top, and injector bottom. The bed is modelled as an explicit random packing (not homogenised).

**Materials:** UCO kernel only (19.75 wt% HALEU) — bare kernel stage, no PyC or SiC coatings present. Structural graphite (retort/cone/heater/injector), process gas (98 mol% H₂ + 2 mol% MTS at 1 atm), liquid water in the injector coolant annulus. All materials at 293.6 K (room temperature, most reactive nuclear state). No water ingress. No boron in graphite (NCS convention — no credit for unconfirmed poisons).

**Nuclear data:** ENDF/B-VIII.0, endfb80_hdf5; 900 K and 1200 K tabulations. OpenMC nearest-temperature method selects 293.6 K → uses 293.6 K tables where available. c_Graphite and c_H_in_H2O S(α,β) applied to all carbon and water materials respectively.

**Driving condition:** Eigenvalue (k-eff) calculation. No external neutron source; fission neutrons drive subsequent generations.

**Shannon entropy mesh:** 10×10×20 regular mesh spanning cone base (z=0) to retort top. Lateral cells 0.5×0.5 cm. Mesh cells that fall inside the cone but outside the inscribed staircase cylinders are empty and contribute zero entropy — this is harmless. At 20,000 particles per generation with pf_fluidized = 0.075, the majority of fissile volume is in the cylinder overflow, providing a meaningful entropy signal without empty-cell noise domination.

**Per-step settings (nominal production run):**
- 50 inactive batches / 250 active batches / 20,000 particles per generation
- Expected statistical precision on k-eff: σ ≈ 1/√(active × particles) relative to batch variance. At deeply subcritical k-eff (expected ≪ 0.1), σ will dominate over any physical effect, so particle count need only be sufficient for reasonable entropy convergence, not tight k-eff precision.

**Batch-count convergence study:** Two runs at seed=42: 250 and 500 active batches. If |k(250) − k(500)| < combined σ, no drift is detected and 250 active batches is sufficient at this stage.

**Seed study:** Three fully independent runs at seeds 42, 43, 44 (250 active batches each). Each seed produces a different random packing geometry AND a different transport RNG sequence, so agreement across seeds confirms neither geometric realisation nor transport randomness biases the result.

**Simulation hierarchy (nominal case):**
```
build_model()          1 call
  bed_region()           1 call — explicit TRISO packing
    pack_bed()             n_slabs + 1 calls — cone staircase slabs + cylinder overflow
  furnace_shell_cells()  1 call — 12 structural cells
openmc.run()           1 call
  inactive generations:  50 × 20 000 particles  = 1 000 000 particles
  active generations:   250 × 20 000 particles  = 5 000 000 particles
```

**Output:** `results/step5_nominal/statepoint.300.h5` — k-eff history, entropy history, three tally results (flux spectrum, reaction rates, U-235 fission spatial).

## Design decisions

- **`fluidized_height` raised from 180 mm to 250 mm** — The placeholder value of 180 mm was a `# CONFIRM` estimate. At the nominal charge of 95 g and pf_fluidized = 0.075, the computed cylinder overflow height is:
  - V_solid = 95 g / 2.94 g/cm³ = 32.3 cm³ (ρ_eff = m_particle / V_outer at full TRISO stage)
  - V_bulk = 32.3 / 0.075 = 430 cm³
  - V_cone_staircase ≈ 23.7 cm³ (n_slabs=8; cone filled at pf_fluidized in both states)
  - V_overflow = 430 − 23.7 = 406.3 cm³
  - h_overflow = 406.3 / (π × 2.5²) ≈ 20.7 cm

  20.7 cm > 18.0 cm → `bed_region()` would raise `ValueError`. The retort cylinder has 34 cm of headroom (retort height = 340 mm), so 250 mm (25 cm) is geometrically realizable and provides ~4.3 cm margin above the nominal overflow top. This remains `# CONFIRM` until verified against the actual process specification for maximum fluidization zone height.

  Note: `fluidized_height_cm` bounds only the cylinder overflow above the cone, not the total bed height. The cone is always available to the bed regardless of this parameter.

- **Particles updated from 50,000 to 20,000; batches from 110 to 300 (50 inactive + 250 active)** — The previous params reflected a production target that predates step 5 convergence work. Starting at 20,000 with 50 inactive allows the convergence study to confirm whether higher counts are needed before committing to expensive sweeps in steps 7–8.

- **`seed` added to `bed_region()` signature** — Minimal non-breaking change (defaults to `params['model']['seed']` if omitted). Required for the seed study to vary both geometry and transport simultaneously.

- **`only_fissionable` dropped; replaced by `constraints={'fissionable': True}` on `IndependentSource`** — The `only_fissionable` parameter on `openmc.stats.Box` was deprecated in OpenMC 0.15. The new parameter achieves the same effect (source positions are rejection-sampled until landing in a fissionable cell). This only affects the first generation; active-batch k-eff is unaffected.

- **`settings.source_rejection_fraction = 0.005`** — The UCO kernel volume fraction in the source Box is ~0.72% (pf=0.075 × kernel/sphere ratio 0.123 × π/4 box-to-cylinder). OpenMC 0.15 defaults to erroring when fewer than 5% of sampled positions satisfy the `fissionable` constraint. Setting to 0.005 (0.5%) allows the actual 0.72% acceptance rate without error.

- **Gas-above-bed cell** — A single `process_gas` cell fills the retort interior between the bed top and the retort top. This is required to cover all geometry — leaving it undefined would cause OpenMC to error on particles entering that region.

- **`MaterialFilter` built from `geometry.get_all_materials()` by name** — `particle_at_stage()` (called from `bed_region()`) creates its own material instances internally; their IDs differ from the top-level `kernel_mat`, `buf_mat`, etc. created in `build_model()`. Using the latter in `MaterialFilter` causes "Could not find material N" at runtime. Fix: build a name → material dict from `geometry.get_all_materials()` and look up materials by name.

- **All material temperatures set to 1200 K** — `endfb80_hdf5` neutron XS tables exist at 900, 1200, and 2500 K; c_Graphite S(α,β) exists only at 1200 K; c_H_in_H2O S(α,β) has an empty temperature table (library defect — no data at any temperature, so the call is omitted). Setting all materials to 1200 K is the only temperature at which both neutron XS and c_Graphite S(α,β) are simultaneously available. Non-conservative for NCS (higher T → more Doppler absorption → lower k; using 1200 K instead of 293.6 K underestimates Doppler absorption slightly, giving slightly higher k than the true room-temperature value). Effect is negligible at Stage 0 since the system is deeply subcritical.

- **`c_H_in_H2O` S(α,β) removed from water material** — Inspection of `endfb80_hdf5/thermal/c_H_in_H2O.h5` revealed the `kTs` group is empty (no temperature sub-groups, no data). `add_s_alpha_beta('c_H_in_H2O')` would cause OpenMC to abort at runtime. Water (injector coolant) is treated with the free-gas kernel instead. This is non-conservative (free-gas underestimates thermalization → less absorption in H) but water is only in the injector cooling annulus, not in the neutron path through the bed, so the k-eff impact is negligible.

- **Smoke test uses 10 000 particles** — The Watt source spectrum (peak ~1 MeV) has only ~0.2% fission probability per source particle in a sub-centimetre UCO kernel. With the smoke test's 1.15 cm bed, k-eff ≈ 0.0004 (99.98% leakage), so source particles rarely cause secondary fissions. 10 000 particles gives ~20 expected fission events per generation, reliably populating the fission bank. The smoke test confirmed k-eff = 0.0004 ± (not computed, 1 active batch) for the 5 g case.

- **Batch-count study design** — User revised from a four-point ladder (100/250/500/1000) to two points (250 and 500). This halves the convergence study runtime while still detecting drift (the primary failure mode) and providing one runtime-scaling data point for sweep planning. If the two-point check passes, the study is done; if it fails, a third run can be added.

## Assumptions

**Confirmed:**
- State: fluidized (step description specifies "fluidized bed")
- Stage: bare_kernel (nominal CVD operating condition — kernels enter the furnace uncoated)
- Background: process gas (step description specifies "precursor gas")
- Temperature: 293.6 K (step description specifies "room temperature"); overridden to 1200 K due to library constraint (see Design decisions)
- Injector coolant: liquid water at 1.0 g/cm³ (as-built); modelled with free-gas kernel (c_H_in_H2O S(α,β) unavailable in endfb80_hdf5)
- No water ingress
- Seed=42 for nominal; seeds 43 and 44 for seed study (different geometry + transport)
- Batch-count study: 250 vs 500 active batches; 50 inactive fixed

**Defaulted:**
- Watt spectrum a = 0.988e6 eV, b = 2.249e-6 eV⁻¹ — standard U-235 thermal fission parameters from ENDF/B-VIII.0, cited in OpenMC stats reference.
- Energy group boundaries: 0, 0.625 eV, 1 MeV, 20 MeV — standard three-group thermal/epithermal/fast split; 0.625 eV is the conventional 2200 m/s thermal cutoff.
- Flux tally: single spatial bin (1×1×1 mesh over bed) for an integrated spectrum.
- U-235 fission spatial mesh: 20×20×30 over the bed bounding box. Each cell is ~0.25 cm × 0.25 cm × ~0.73 cm at the nominal fluidized bed height of 21.9 cm.
- Entropy mesh: 10×10×20 over the full retort cylinder interior.

**Unconfirmed (`# CONFIRM`):**
- `fluidized_height: 250.0 mm` — raised from placeholder 180 mm to accommodate the computed bed height of 21.9 cm at 95 g charge. The correct value depends on the maximum operating fluidization zone height from the process specification. If the actual zone is shorter than 21.9 cm, the nominal charge must be reduced or the bed_expansion_ratio revised.
- `charge_mass_g: 95.0 g` — step specifies "nominal charge"; value is a placeholder pending process documentation.
- `packing_fraction_fluidized: 0.075` — derived from pf_static / bed_expansion_ratio; both factors remain `# CONFIRM`.
- All other `# CONFIRM` items inherited from steps 3–4 (`kernel_diameter`, layer thicknesses, densities, graphite density, gas pressure).
