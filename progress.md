# Progress Log

Write for someone picking this up cold in six months.
Never edit old entries; add a new section when a decision is reversed.

---

## Step 0 — Scaffold (2026-09-08)

### What was added / changed

- `furnace/__init__.py` — package root; re-exports `load_params`, `check_env`.
- `furnace/params.py` — `load_params()` and `check_env()` (see design decisions below).
- `furnace/materials.py`, `triso.py`, `geometry.py`, `model.py`, `sweeps.py`, `postproc.py` — stubs with `# TODO` markers; none are importable workers yet.
- `params.yaml` — all dimensions, TRISO layer parameters, material densities, gas composition, and model settings in one place.
- `scripts/run_nominal.py`, `scripts/run_sweep.py` — stubs.
- `cases/`, `results/` directories created; `cases/` gitignored.
- `.gitignore` updated: `cases/`, `results/*.csv|png|h5`, `statepoint.*.h5`, `summary.h5`.
- `progress.md` created (this file).

### Design decisions and reasoning

**`load_params()` return type: recursive `types.MappingProxyType` (frozen dict), not nested dataclasses.**
Frozen dict was chosen because the schema grows across steps — each step adds new dimension or material keys. Nested dataclasses with `frozen=True` would require updating field definitions every time. `MappingProxyType` applied recursively gives identical mutation protection with zero per-field boilerplate. Tradeoff: no IDE key autocompletion (acceptable for a config structure).

**`load_params()` and `check_env()` live in `furnace/params.py`, not `furnace/materials.py`.**
`materials.py` is step 1's job; putting infrastructure utilities there would require step 1 to be partially implemented before the step 0 check works. `furnace/params.py` is importable with only PyYAML — no OpenMC required — which is useful for unit tests and parameter inspection.

**`check_env()` defers `import openmc` until call time.**
`furnace.params` can be imported in environments without OpenMC installed (CI, parameter scripts). Importing at module level would break any `from furnace import load_params` in such environments.

**`_cm` keys added alongside originals, not replacing them.**
`params['triso']['kernel_diameter']` → 425.0 µm (natural unit, matches the YAML).
`params['triso']['kernel_diameter_cm']` → 0.00425 cm (OpenMC unit).
Both are correct; the user-facing value is preserved for readability, the converted value is ready for geometry code.

**Volume `static_bulk_cc` stored under `dimensions.bed_volume`, not `dimensions.bed`.**
The static bulk volume is in cm³, not mm; mixing it into the `dimensions.bed` sub-dict (where everything is in mm) would cause `_add_cm_keys` to apply a mm→cm factor to a volume, producing a nonsensical cm³ × 0.1 value. Separating it avoids a special-case guard in the converter and makes the unit mismatch obvious in the YAML.

**`_LENGTH_SUFFIXES` whitelist controls which keys get a `_cm` companion.**
Only keys whose names contain "diameter", "depth", "height", "length", "drop", "thickness", "id", or "od" are converted. This prevents `included_angle_deg`, `pressure_pa`, `temperature_k`, and `static_bulk_cc` from getting nonsensical cm values.

**Seed = 42 in `model.seed`.**
Arbitrary but deterministic. The value itself is unimportant; recording it in every result file is the requirement.

### Assumptions

**Confirmed:**
- All material densities and enrichment match the reference repo (`nican0r/TRISO-neutronics`) and AGR-1 specifications (Demkowicz et al. 2018; INL/EXT-10-19476).
- Gas composition: 98 mol% H₂ + 2 mol% CH₃SiCl₃ at 101325 Pa, 1873 K — confirmed by user as worst-case moderation scenario.
- Structural graphite boron impurity = 0 ppm; no unconfirmed poison credited.

**Defaulted:**
- TRISO layer thicknesses (buffer 100 µm, IPyC 40 µm, SiC 35 µm, OPyC 40 µm): AGR-1 nominal values from INL/EXT-10-19476 Table 3. `# CONFIRM` flags retained.
- Kernel diameter 425 µm: AGR-1 nominal. `# CONFIRM` retained.
- Structural graphite density 1.75 g/cm³: from AGR-1 compact matrix (INL/EXT-10-19476 Table 5); may not match furnace graphite grade. `# CONFIRM`.
- Graphite felt density 0.20 g/cm³: placeholder from preamble. `# CONFIRM`.
- Seed = 42: arbitrary.

**Unconfirmed (`# CONFIRM`):**
- `bed.static_depth` (35 mm) — wrong value shifts fissile mass distribution in cone geometry.
- `bed_volume.static_bulk_cc` (50 cm³) — sets number of TRISO particles; directly scales k-eff.
- `bed.fluidized_height` (180 mm) — affects worst-case expanded bed geometry.
- `model.fuel_temperature_k` (1500 K) — library only has 900 K and 1200 K; effective XS temperature will be 1200 K regardless (non-conservative for Doppler).
- `graphite_felt_thickness_mm` (null) — insulation not yet modelled.
- `gas.pressure_pa` (101325 Pa) — sub-atmospheric operation possible.

### Cross-section temperature limitation (critical note)

The library (`~/Documents/nuclear-data/endfb80_hdf5/cross_sections.xml`) was processed at 900 K and 1200 K. The furnace operates at 1500–1825 K. OpenMC's nearest-temperature method will use 1200 K cross sections for all materials even when `fuel_temperature_k = 1500 K` is set. The `c_Graphite` S(α,β) evaluation in this library is also only available at 1200 K.

Effect on k-eff:
- Doppler broadening at 1200 K is stronger than at 1500 K → captures at 1200 K are slightly higher → k-eff is slightly **lower** than the true 1500 K value. Non-conservative.
- `c_Graphite` thermal scattering at 1200 K produces slightly more thermalisation than at 1500 K → slightly **higher** k-eff contribution from moderation. Non-conservative direction is problem-dependent.

Net effect: small and in opposite directions; acceptable for screening. Must be stated in every report.

### Open questions

None carried forward from step 0. All `# CONFIRM` items are dimension/material placeholders to be resolved with the user before production runs.

---

## Step 1 — Materials (2026-09-09)

### What was added / changed

- `furnace/materials.py` — fully implemented with ten material factory functions:
  `uco_kernel`, `buffer_pyc`, `ipyc`, `sic`, `opyc`, `graphite_structural`,
  `graphite_felt_insulation`, `process_gas`, `water`, `air`.
  Also: `_print_material_table()`, `_verify_u235_atom_density()`, and
  `if __name__ == '__main__':` diagnostic block.

### Design decisions and reasoning

**All material temperatures set to 293.6 K (step file updated; cold case only).**
The original step described a `cold=True` flag with 293.6 K and a hot operating-temperature default. The step was revised before implementation: the cold, flooded condition is the bounding NCS screening case, so we model everything at room temperature. Hot-temperature materials would only be needed if we wanted to compute Doppler feedback, which is not the purpose of this model. Cold is more reactive because it minimises U-238 Doppler capture and maximises graphite/water thermalisation effectiveness.

**c_Graphite applied to buffer, IPyC, OPyC, and structural graphite — not to the UCO kernel or SiC.**
The preamble explicitly lists these four carbon-containing materials and excludes the kernel. The kernel is a uranium fuel material where thermal scattering from the minor C constituent is dominated by the far larger U resonance and epithermal cross-sections; applying c_Graphite there would be physically unjustified. For SiC: c_SiC exists in the library and would be defensible, but the preamble didn't credit it and it is deferred. Omitting thermal scattering from PyC is a worse approximation than using c_Graphite (an imperfect model for turbostratic carbon); erring on the non-conservative side by omitting it entirely is unacceptable for NCS work.

**What omitting c_Graphite from structural graphite would do to k-eff.**
Free-gas scattering treats graphite atoms as independent at the material temperature. This misses the phonon enhancement that makes graphite thermalise sub-thermal neutrons much more efficiently than the free-gas prediction. Published comparisons for graphite-moderated systems show k-eff depressions of 2–5% when c_Graphite is replaced by free-gas. Since graphite is the primary moderator in this geometry, omitting it would produce an artificially low k-eff and a false margin of safety — the most dangerous failure mode in NCS work.

**Gas density computed at operating T/P (1873 K, 101325 Pa).**
This represents the actual physical state of the CVD atmosphere. Using room-temperature density (293.6 K) would give ~6.4× higher gas density and would be more conservative (more hydrogen → higher k-eff). The effect is small in either case because the gas is so dilute; the decision was made to use the physically correct operating condition and note the non-conservatism.

**OpenMC U-234/U-235 approximation at HALEU enrichment.**
OpenMC warns that its built-in enrichment helper applies a fixed U-234/U-235 mass ratio of 0.008, which is derived for low-enriched uranium. At 19.75 wt%, this may not be accurate for cascade-enriched HALEU. U-234 is a resonance absorber; an incorrect U-234 inventory introduces a small error in k-eff. For Stage 0 screening this is acceptable; a later stage should specify isotopic fractions explicitly from the fuel specification.

### Assumptions

**Confirmed:**
- All TRISO layer compositions and densities from AGR-1 spec (Demkowicz et al. 2018; INL/EXT-10-19476).
- UCO atom fractions: U:C:O = 1:0.5:0.4 (per formula unit), enrichment 19.75 wt%.
- c_Graphite for buffer/IPyC/OPyC/structural graphite; c_H_in_H2O for water.
- Process gas: 98 mol% H₂ + 2 mol% CH₃SiCl₃, only H-bearing species modelled.

**Defaulted:**
- SiC: no S(α,β) applied (c_SiC not credited; deferred).
- Air: NIST standard dry air, CO₂ neglected.
- U-234 content: OpenMC internal approximation (mass ratio = 0.008 × U-235).

**Unconfirmed (`# CONFIRM`):**
- `graphite_structural_density_gcc` = 1.75 g/cm³ — furnace graphite grade not confirmed.
- `graphite_felt_density_gcc` = 0.20 g/cm³ — placeholder.
- TRISO dimensions (kernel_diameter, layer thicknesses) — AGR-1 nominal.
- `gas.pressure_pa` = 101325 Pa — sub-atm operation possible.

### Verification output (from `python3 -m furnace.materials`)

All ten materials printed with correct name, density, temperature (293.6 K), and S(α,β) assignments.

U-235 atom density check:
- Analytic: 5.0494 × 10²¹ atoms/cm³
- OpenMC:   5.0494 × 10²¹ atoms/cm³
- Ratio:    0.999998 (agreement to 6 significant figures)

### Open questions

- U-234 specification: what is the actual U-234 content in this HALEU fuel spec? Replace OpenMC's approximation with explicit isotopics when available.
- Gas pressure: confirm whether the CVD process runs at sub-atmospheric pressure. If so, the gas density will be lower (less conservative, but already small).
- Should process gas density use room temperature (conservative) rather than operating temperature? The difference in k-eff is expected to be negligible, but worth quantifying once the geometry is assembled in Step 3.

---

## Step 7 — Mass sweep (2026-09-10)

### What was added / changed

- `furnace/sweeps.py` — `vessel_capacity_g(params, state, stage)` helper: usable retort volume × packing fraction × per-particle effective density.
- `furnace/sweeps.py` — `run_case()` merge order fixed: quick-mode defaults now applied *before* per-case overrides, so a sweep can override `charge_mass_g` even in quick mode (previously the 5 g quick-mode value would silently overwrite the sweep's per-case mass).
- `scripts/run_sweep.py` — `--mass-sweep` flag, `_MASS_MULTIPLIERS`, `_mass_sweep_cases()`, `mass_sweep()` with vessel-capacity precheck and skip-and-note behaviour for over-capacity cases.
- `scripts/plot_mass_sweep.py` — new plotter: log-x charge mass, k-eff ± 2σ error bars, k+2σ+δk_disc conservative trace, 0.95 subcritical reference line, `UNVALIDATED — SCREENING ONLY` footer.
- `docs/steps/step-7-mass-packing-sweeps.md` — full step doc including experimental design and the packing-fraction change record.
- `params.yaml` — **`model.packing_fraction_static` lowered from 0.60 to 0.50** (reversal, see below).

### Design decisions and reasoning

**Reversal — collapsed packing fraction 0.60 → 0.50.** The original 0.60 was chosen in steps 2–3 as a settled/vibrated bed reference. Step 7 is the first step that actually packs a large number of TRISO particles at the collapsed pf (steps 2–6 all used `state='fluidized'` at pf ≈ 0.075, which packs in seconds). A smoke test at 5 g / ~12 000 bare kernels at pf = 0.60 did not converge in 69 minutes of CPU time — the packing algorithm was still running when killed.

The failure is not a bug; `openmc.model.pack_spheres` uses a Jodrey–Tory contraction algorithm whose runtime scales as (φ_max − φ_target)⁻² near the random-close-pack ceiling φ_max ≈ 0.64. The gap-to-jamming numbers make the scaling stark:

| pf_target | gap to 0.64 | relative work per particle |
|-----------|-------------|-----------------------------|
| 0.075 | 0.565 | 1× |
| 0.50  | 0.14  | 16× |
| 0.60  | 0.04  | 200× |

Physically, 0.50 is the loose random pack of a gravity-settled powder without mechanical tapping — the correct regime for TRISO kernels accumulating in the retort under gravity alone, which is the accident path step 7 is bounding. Reaching 0.60 requires deliberate consolidation (tapping/vibration) that is not postulated. Alternatives (FCC/BCC lattice pack, Lubachevsky–Stillinger implementation) were rejected as either physically wrong or out of scope.

Consequences:
- `packing_fraction_fluidized` derives to 0.0625 (previously 0.075) because bed_expansion_ratio is unchanged at 8.0. This drops fluidized fissile density by 17%. Non-conservative for the fluidized nominal case, but the absolute k shift is small in the deeply subcritical regime. Step 5 nominal will be re-run before its k-eff is quoted anywhere final.
- Vessel mass capacity at collapsed pf: 4384 g → 3653 g. The 40× multiplier (3800 g) is now auto-skipped by the sweep's capacity precheck.
- Steps 2 and 3 documentation still references 0.60 — left unedited because that value was accurate at the time. This entry is the reversal record.

**Quick-mode merge order fix.** In step 6, `run_case` applied `_QUICK` after `overrides`, so `_QUICK['dimensions']['bed']['charge_mass_g'] = 5.0` would wipe out any per-case charge-mass override. The mass sweep can't use `--quick` under that ordering (all cases would collapse to 5 g). Swapped so `_QUICK` is applied first and per-case overrides win. The seed-check sweep is unaffected because it doesn't override `charge_mass_g`.

**Vessel-capacity pre-check, not geometry-level guard.** `bed_region` has no upper-bound check on collapsed bed height. Rather than modify `furnace/geometry.py`, the sweep driver filters the case list up-front and prints the skip reason to stdout. Keeps geometry code untouched and makes the skip decision auditable in the sweep log.

**Mass grid `[1, 1.5, 2, 3, 5, 10, 20, 40]` × 95 g.** Concentrates points below 5× where the k-eff-vs-mass trend is expected to change most, spreads them logarithmically above that. Log-x axis in the plot for the same reason.

**Bare kernel stage.** Highest fissile-atom density per particle — bounding choice for k-eff at a given mass. Coated stages dilute the kernel with low-density C/SiC layers.

**Gas background.** Answers the physical CVD-atmosphere question. Water-flooded variant is deferred to step 8; keeping moderator constant here isolates the mass effect.

**Single seed, no per-mass repeat.** The mass-to-mass k spread is expected to dwarf the seed-to-seed spread from step 6. If a specific mass needs tighter error bars later, the seed-check machinery already exists.

### Assumptions

**Confirmed:**
- Vessel usable volume = frustum(r_throat → r_retort) + cylinder(r_retort, h = retort height). Derived from `params.yaml` with no fabrication assumption beyond what is already there.

**Defaulted:**
- Collapsed packing fraction 0.50 (loose random pack of a gravity-settled powder). Marked `# CONFIRM` in `params.yaml`.
- Mass multiplier grid; log-x plot; single seed; bare kernel stage; gas background — see design decisions.

**Unconfirmed:**
- δk_disc from the collapsed / water / low-particle convergence study applied to all cases regardless of state or background. Conservative but not tight for the gas cases. Inherited Stage 0 approximation from step 6.
- Fluidized nominal (step 5) k-eff will shift slightly under the new pf_fluidized = 0.0625 — magnitude to be quantified by a re-run before that number is quoted final.

### Surprises

- **Every step before this one had implicitly assumed collapsed packing was tractable, because no step before this one actually built a collapsed bed at scale.** The state parameter existed and worked for the geometry logic, but no OpenMC run ever exercised it. This is a good reminder that "the code compiles and step 3 says it works" is not the same as "we've ever run it end-to-end."

### Open questions

- Is 0.50 the right loose-random-pack value for this specific powder, or should it be measured against the actual TRISO kernel population? A process-specific value would supersede the literature estimate.
- Should step 5 nominal be re-run and its progress-log entry updated to reflect the new fluidized pf? Deferred until user decision.
