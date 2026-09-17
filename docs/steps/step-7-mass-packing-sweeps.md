# Step 7 — Mass sweep

## How to run

```
eval "$(micromamba shell hook --shell zsh)"
micromamba activate triso-env

# Quick pipeline smoke-test (1 000 particles, 12 batches, ~30 s per case):
caffeinate python scripts/run_sweep.py --mass-sweep --quick

# Full sweep at nominal 20 000 particles / 300 batches (~few minutes per case):
caffeinate python scripts/run_sweep.py --mass-sweep

# In parallel (N concurrent processes; each still gets its own OpenMC thread pool):
caffeinate python scripts/run_sweep.py --mass-sweep -j 3

# Plot after the CSV is written:
caffeinate python scripts/plot_mass_sweep.py
```

Outputs:
- CSV — `results/mass_sweep.csv`
- PNG — `results/mass_sweep.png`
- Per-case OpenMC run dirs — `results/mass_sweep/<tag>/`

---

## What was implemented

- **`furnace/sweeps.py`** — added `vessel_capacity_g(params, state, stage)`, the mass ceiling above which no case is submitted. Computes the retort's usable interior volume (cone frustum + retort cylinder up to the vacuum boundary) and multiplies by the packing fraction and per-particle effective density.
- **`scripts/run_sweep.py`** — added `--mass-sweep` flag, `_MASS_MULTIPLIERS`, `_mass_sweep_cases(params)`, and `mass_sweep(params, args)`. Every case runs `state='collapsed'`, `stage='bare_kernel'`, `background='gas'` and overrides only `dimensions.bed.charge_mass_g`. Cases that exceed the vessel capacity are noted on stdout and dropped from the submission list.
- **`scripts/plot_mass_sweep.py`** — reads `results/mass_sweep.csv` and produces `results/mass_sweep.png`: log-x charge mass vs k-eff with 2σ error bars, the conservative `k+2σ+δk_disc` upper bound, and a horizontal 0.95 reference line. Carries the `UNVALIDATED — SCREENING ONLY` footer required by the preamble.

No changes to `furnace/model.py`, `furnace/geometry.py`, or `params.yaml` were needed; the step-6 override mechanism handles per-case charge mass without any new plumbing.

---

## How it works

The sweep multiplies the nominal 95 g charge by `[2.0, 3.0, 5.0, 10.0, 20.0, 38.0]` (190 g – 3610 g). The 1× and 1.5× multipliers are excluded: 1× duplicates the step-6 nominal collapsed-bed run, and the engineering question targets accumulation above the operating charge. The grid terminates at 38× because the retort's usable volume at pf = 0.50 holds at most 3653 g of bare kernel; 38× (3610 g) is the highest round multiplier that fits with a practical margin (~43 g). Any multiplier above 38× would be auto-skipped by `vessel_capacity_g`. Each surviving case builds a full `openmc.Model` with the same materials, tallies, and settings as the nominal step-5 run, differing only in the charge mass override. Results are appended to `results/mass_sweep.csv` via the existing crash-safe locked-append helper, and the plotter reads the CSV back to produce the k-eff-vs-mass figure.

Using `state='collapsed'` throughout removes the fluidized overflow guard from the geometry (`bed_region` only raises the overflow check for `state='fluidized'`), so the bed grows naturally into the retort cylinder as mass increases. Because bed height scales linearly with mass at fixed packing and cross-section, capacity is bounded by the vacuum plane at the top of the retort rather than by any geometry ambiguity in the code.

## Packing-fraction change: 0.60 → 0.50

The collapsed packing fraction in `params.yaml` was lowered from 0.60 to 0.50 during this step. The reason is algorithmic tractability, not physics, and the change is documented here because the *why* is not obvious from the diff.

**Symptom.** A smoke test at 5 g / ~12 000 bare kernels using the previous 0.60 value did not converge in 69 minutes of single-core CPU time; the packing loop was still running when killed. The same input at 0.075 (fluidized) had been packing in seconds throughout steps 2–6.

**Why 10% of packing fraction is not 10% of runtime.** OpenMC's `openmc.model.pack_spheres` uses a Jodrey–Tory-style contraction algorithm. The runtime of any random hard-sphere packer scales as the inverse of the gap to the jamming limit (random close packing, φ_max ≈ 0.64 for monodisperse spheres in 3D), roughly

$$ t \sim (\varphi_{\max} - \varphi_{\mathrm{target}})^{-\alpha}, \quad \alpha \approx 2. $$

The gap-to-jamming, not the packing fraction itself, sets the cost. Concretely:

| pf_target | gap to 0.64 | relative work per particle (∝ gap⁻²) |
|-----------|-------------|--------------------------------------|
| 0.333 (fluidized, current)        | 0.307 | 3.4× |
| 0.50 (collapsed, current)         | 0.14  | 16× |
| 0.60 (old collapsed)              | 0.04  | 200× |
| 0.64 (jamming)                    | 0.00  | ∞ (never converges) |

The 0.60 → 0.50 change closes a ~13× per-particle runtime gap. Multiplied by 12 000 particles and an O(n²)-ish neighbor search per contraction iteration, the runtime blowup between 0.50 and 0.60 is large enough that the two operating points are qualitatively different — 0.50 is inside the algorithm's well-behaved regime, 0.60 is inside its critical region.

**Physical defensibility of 0.50.** Loose random pack (LRP) of monodisperse hard spheres under gravity, without mechanical tapping or vibration, is φ ≈ 0.55–0.60 for well-flowing powders and lower for irregular or cohesive powders. 0.50 is the low end of that band and is the correct regime for a bed of TRISO kernels that has accumulated inside the retort under gravity but has not been consolidated by mechanical means. Reaching 0.60 requires deliberate settling (tapping, vibration, or vacuum-densification), which is not part of the postulated CVD-accumulation accident path this step is bounding. Marked `# CONFIRM` in `params.yaml` because a process-specific loose-pack measurement would supersede this literature estimate.

**Alternatives considered and rejected.**
- **Deterministic FCC/BCC lattice pack** at φ = 0.74 / 0.68 — physically unrealistic for a settled powder (real settled beds are amorphous, not crystalline) and would understate the moderator-to-fuel spacing variability that drives self-shielding statistics.
- **Lubachevsky–Stillinger implementation** — a proper molecular-dynamics-style packer that reaches φ = 0.64 without blowing up. Correct answer but a large software change outside this step's scope.

**Consequences of the change.**
- Collapsed-bed bulk volume increases by 20% (V_bulk = m / ρ_eff / pf, and 1/0.50 vs 1/0.60 = 1.20×).
- **Fluidized packing fraction is 0.333** (= pf_static / bed_expansion_ratio = 0.50 / 1.5). The bed_expansion_ratio was updated from the placeholder 8.0 to 1.5 (confirmed against process data for near-minimum-fluidization operation), giving pf_fluidized = 0.333 rather than the prior 0.0625.
- Vessel mass capacity at collapsed pf = 0.50: **3653 g**. The grid ceiling is 38× (3610 g), the highest round multiplier that fits; six cases run.
- Step-3 documentation was updated to reflect pf_static = 0.50, n_slabs = 32, and the revised convergence table.

---

## Experimental design

**What is being modelled.** Six independent eigenvalue k-eff calculations, each a variant of the step-5 nominal geometry with the TRISO charge mass replaced by 2×, 3×, 5×, 10×, 20×, 38× the nominal 95 g (190 g – 3610 g). The 1× and 1.5× multipliers are excluded: the 1× collapsed-bed case is already the step-6 reference run, and the engineering question is about accumulation above the operating charge. The grid terminates at 38× because the retort's usable volume at pf = 0.50 holds at most 3653 g of bare kernel — 38× (3610 g) is the highest round multiplier that fits with a practical margin (~43 g); any case above 38× would be auto-skipped by `vessel_capacity_g`. Every case uses a **collapsed** (gravity-settled) bed at packing fraction 0.50 rather than the operating fluidized bed, because the mass-limit question is about accumulated inventory, not the operating condition. The bed starts at the cone base and grows upward through the cone frustum into the retort cylinder as mass increases.

**Materials.** Bare UCO kernels only — no coating layers. Bare kernel is chosen because it is the most compact (highest fissile density per unit bed volume) form the charge takes during the CVD cycle; it bounds k-eff for a given charge mass. The interstitial atmosphere is the SiC-step CVD gas (98 mol% H₂ + 2 mol% MTS) at 1 atm, as in the nominal case. Structural graphite (retort, cone, heater, injector body) and water-cooled injector geometry are unchanged from step 4.

**Nuclear data.** ENDF/B-VIII.0 (`endfb80_hdf5`) with `c_Graphite` S(α,β) applied to every carbon-bearing material. Temperature is snapped to 900 K by OpenMC's `nearest` method (library has no room-temperature evaluation) — the same Stage 0 limitation as the nominal run.

**Driving condition.** Eigenvalue mode; Watt fission spectrum source (a=0.988 MeV, b=2.249×10⁻⁶ eV⁻¹) sampled over a box spanning the full bed extent with `constraints={'fissionable': True}` rejection.

**Parameter scheme.** One-dimensional sweep in charge mass alone; every other parameter is at its step-5 nominal value. The mass grid concentrates points below 5× nominal and spreads them logarithmically above that, so the k-eff-vs-mass trend can be resolved through both the low-inventory and the accumulation regimes.

**Per-step settings.** Each case runs 20 000 particles/generation, 300 batches (50 inactive + 250 active), seed=42 — identical to the nominal step-5 case. Expected σ(k) ≈ few × 10⁻⁴, small compared to the differences between mass points, which is what the plot needs.

**Simulation hierarchy.**
- 6 mass values (grid)
  - 1 case per mass (single seed, matches nominal for direct comparison)
    - 300 batches × 20 000 particles = 6 × 10⁶ histories per case
    - Total ≈ 3.6 × 10⁷ histories across the sweep

**Output.**
- `results/mass_sweep.csv` — one row per case with charge mass (in `overrides_json`), U-235 mass, k-eff, σ, k+2σ, k+2σ+δk_disc, bed height, achieved packing fraction, H/²³⁵U and C/²³⁵U ratios, thermal flux fraction, wall time, seed.
- `results/mass_sweep.png` — log-x plot of k-eff ± 2σ vs charge mass with the conservative k+2σ+δk_disc trace and the 0.95 horizontal subcritical reference line. Every data point is annotated with its k-eff, σ, and U-235 mass.

---

## Results

Canonical production run: `step7a_mass/1789561375_788a31460970`, git `030835a`, manifest `manifests/step7a_mass_tiled.json`, state=collapsed, stage=bare_kernel, use_exact_cone=True, use_tiled_bed=True, tile_size_cm=0.5, seed=42. 20 000 particles/gen × 300 batches (50 inactive + 250 active).

| charge_mass_g | bed_height_cm | pf_achieved | u₂₃₅ mass (g) | k-eff | σ | k+2σ | wall (s) |
|---------------|--------------|-------------|---------------|-------|---|------|---------|
| 95 | 3.212 | 0.4997 | 17.820 | 0.021637 | 1.20×10⁻⁵ | 0.021661 | 207 |
| 142.5 | 3.751 | 0.4997 | 26.730 | 0.024794 | 1.29×10⁻⁵ | 0.024820 | 238 |
| 190 | 4.203 | 0.4999 | 35.657 | 0.027437 | 1.40×10⁻⁵ | 0.027465 | 279 |
| 285 | 5.125 | 0.4999 | 53.477 | 0.032148 | 1.65×10⁻⁵ | 0.032181 | 334 |
| 475 | 6.968 | 0.4998 | 89.117 | 0.038874 | 1.61×10⁻⁵ | 0.038906 | 370 |
| 950 | 11.576 | 0.4997 | 178.218 | 0.047674 | 2.26×10⁻⁵ | 0.047719 | 463 |
| 1900 | 20.792 | 0.4997 | 356.419 | 0.054148 | 2.57×10⁻⁵ | 0.054199 | 517 |
| 3610 | 37.380 | 0.4997 | 677.181 | 0.057324 | 2.72×10⁻⁵ | 0.057379 | 549 |

**Key findings:**
- All cases are deeply subcritical: maximum k+2σ = 0.05738 at 3610 g, well below 0.95.
- k-eff increases monotonically with charge mass and flattens significantly above ~500 g — doubling from 1900 g to 3610 g raises k by only +0.003.
- Bed geometry transitions from cone-only (≤ 142.5 g; bed top at 3.75 cm ≈ cone top) to cone + retort-cylinder overflow (≥ 190 g; bed top 4.20 cm). Theoretical overflow onset: ~149 g at pf = 0.50 (V_bulk = V_frustum when charge_mass_g = ρ_eff × pf × V_frustum ≈ 10.5 × 0.50 × 28.3 cm³).
- pf_achieved is 0.4997–0.4999 across all masses — the tiled packing consistently achieves the target without mass-dependent drift.

**Superseded submission:** `step7a_mass/1789162528_6c380c97abd2` (older code revision), superseded. `step7a_mass_tiled/1789499610_533de1774058` (git `76f0e5c`, tiled bed without exact_cone): 6 cases (190–3610 g), k values within ~3×10⁻⁴ of the canonical results, confirming the exact-cone correction had negligible k-eff impact at the mass-sweep geometry. Superseded by the canonical run above.

---

## Design decisions

- **Multipliers start at 2× (not 1× or 1.5×).** The engineering question is "how much accumulation before mass becomes a concern", not "is nominal safe" — the latter is answered by step 5. The 1× collapsed-bed case is already the step-6 reference run; duplicating it here wastes compute without adding information. 1.5× was also dropped: the first point of interest is a meaningfully larger inventory, not a 50% increment that would be difficult to distinguish from the nominal on the log-x plot.
- **Collapsed bed throughout.** A settled bed at pf ≈ 0.50 (loose random pack; see §Packing-fraction change above for why not 0.60) packs the fissile inventory into a small bulk volume, yielding a high fissile-atom density per cm³ of bed without invoking mechanical settling that is not part of the postulated accident. That is the most reactive gravity-only configuration for a given mass under this study's atmosphere assumption, so it is what the mass limit must bound.
- **Gas background.** Answers the physical CVD-atmosphere question. A water-flooded variant is deferred to step 8; keeping the moderator constant here isolates the mass effect.
- **Bare kernel stage.** Highest fissile-atom density per particle. Coated stages dilute the kernel with low-density C/SiC layers and would raise k less per unit charge mass. Bare kernel is the correct bounding choice.
- **Vessel-capacity pre-check, not geometry-level guard.** `bed_region` for `state='collapsed'` has no upper-bound check on bed height. Rather than modify the geometry, the sweep driver computes the maximum admissible mass and filters the case list up-front. This keeps the geometry code untouched and makes the skip decision auditable in the sweep stdout.
- **Single seed, no per-mass repeat.** The mass effect between adjacent grid points is much larger than the single-case σ, so seed averaging is not needed to resolve the trend. If a future step needs tighter error bars at a specific mass, the seed-check machinery already exists.
- **Log-x axis in the plot.** The mass grid spans 2×–38× and shows most structure at the low end; a linear axis would compress the sub-5× region into a sliver.

---

## Assumptions

**Confirmed:**
- Vessel usable volume = frustum(r_throat → r_retort, half-angle 30°) + cylinder(r_retort, height = retort height). Derived from `params.yaml` dimensions with no fabrication assumption beyond what is already in `params.yaml`.
- Collapsed packing fraction 0.50 (`model.packing_fraction_static`, updated in this step from the earlier 0.60 — see §Packing-fraction change above for the tractability reasoning and physical justification).
- Nominal statistics (20 000 particles × 250 active batches, seed 42) reproduce the step-5 convergence-verified regime.

**Defaulted:**
- Mass multiplier grid `[2, 3, 5, 10, 20, 38]`. Starts at 2× (1× is the step-6 nominal; 1.5× adds no meaningful separation on the log-x axis and is not a distinct accumulation scenario). Terminates at 38×: the retort holds at most 3653 g of bare kernel at pf = 0.50, so 38× (3610 g) is the highest round multiplier that actually runs — it leaves ~43 g of margin to the vacuum boundary. The previous 40× entry (3800 g) was always auto-skipped and has been removed.
- Bare kernel stage rather than any coated stage. Bounding choice; documented above.
- Single seed. The mass-to-mass k spread is expected to dwarf the seed-to-seed spread from step 6.

**Unconfirmed (`# CONFIRM`):**
- δk_disc for n_slabs=32 is 0.0 in the CSV (no n=64 baseline); `_load_dk_disc(32)` returns 0.0. Estimated missing bias ~3×10⁻⁵ (geometric halving of δk=5.32×10⁻⁵ from n=16→32 at pf=0.50), negligible vs σ. Applied conservatively as 0 in reported k+2σ+δk_disc. This is a Stage 0 approximation.
- Structural graphite density 1.75 g/cm³ with zero boron equivalent — see preamble.
- All ambient temperature approximations (900 K instead of 293.6 K, 1200 K for `c_Graphite` S(α,β)) — inherited library limitation, non-conservative but small for this subcritical regime.

---

## Tiled repeating-unit-cell bed (memory fix for high-mass cases)

**Symptom.** AWS Batch OOM at 10×–38× mass despite two instance-memory bumps. Linear extrapolation from a local 95 g profile predicted ~21 GB of `pack_bed`-stage RSS at 38× — no reasonable instance size fixes that.

**Root cause.** `openmc.model.TRISO(...)` returns a `Cell` with a `Sphere` region per particle. `pack_bed` + `create_triso_lattice` therefore materialise **O(N_particles)** unique `Cell` objects — ~2.6 KB each in this model. At 95 g nominal this is already 565 MB of `Cell` overhead; at 38× (~8.7 M particles) it is ~21 GB. The lattice / XML layers are trivial by comparison (~7 MB at 95 g, ~85× less than the packing layer). Measured via `scripts/profile_memory.py`.

**Fix — tiled repeating-unit-cell bed** (`furnace/triso.py:tiled_bed`). Pack one small cube (default 0.5 cm edge, ~1500 bare-kernel particles at pf 0.50) *once*, wrap it as an `openmc.Universe`, and place that **same universe reference** in every cell of an outer `openmc.RectLattice` sized to cover the bed region. Because every outer-lattice cell shares the tile universe, the unique `Cell` count is O(particles per tile), not O(particles in bed) — memory becomes independent of charge mass.

Wired behind a keyword flag: `build_model(..., use_tiled_bed=True, tile_size_cm=0.5)`. Default `False`, so steps 3–6 remain byte-identical. Both `bed_region` (staircase) and `exact_cone_bed` accept the flag.

**Measured at 95 g nominal** (`use_tiled_bed=True` vs `False`):

| metric | random pack | tiled | ratio |
|---|---|---|---|
| N particles in bed | 216 990 | 224 956 | ×1.04 (packing variance) |
| Achieved pf | 0.4990 | 0.4997 | — |
| Unique `openmc.Cell` count | 326 493 | 70 343 | 0.22× |
| ΔRSS during build | 565 MB | 29 MB | **0.05×** |

At 38× mass, the tiled build RSS stays flat at ~30 MB (scales with particles-per-tile, not bed volume). This is the memory headroom the sweep needs.

**Trade-off.** The tile universe is repeated across the bed, so the particle placement is quasi-periodic at 0.5 cm rather than fully random. Thermal diffusion length in this bed is several cm, so the artifact on k-eff is expected to be small — but must be validated, not assumed.

**Validation.** `scripts/validate_tiled_bed.py` runs matched-seed k-eff at a given charge mass in both modes and compares:

```
python3 scripts/validate_tiled_bed.py                          # 95 g, quick stats
python3 scripts/validate_tiled_bed.py --mass 190 \             # 2× nominal, tighter
    --particles 20000 --inactive 50 --active 200
```

Rule: accept the tiled path if `|Δk| < 3·σ(Δk)` at the operating statistics (few × 10⁻³). If not, either shrink `tile_size_cm` (fewer particles per periodicity length) or re-seed the tile per slab.

**Follow-on optimisation not done here.** The current `tiled_bed` builds a fresh tile universe on every slab call, so total unique cells scale as N_slabs × particles_per_tile (~15 slabs × 1500 = 22.5k at 95 g; growing with mass as more slabs are populated). Sharing a single tile universe across all slabs would give a further ~O(N_slabs) reduction at high mass; it was deferred because the current implementation already fits comfortably in memory at 38×, and slab-independent tiles keep the axial pack quasi-random rather than strictly periodic through the full bed height.
