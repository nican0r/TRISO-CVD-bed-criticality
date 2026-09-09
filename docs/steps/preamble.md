---
step: preamble
title: Requirements, working agreement, and physical description
status: complete — no code produced
---

## What was implemented

- Physical system described and requirements locked in (see below).
- Working agreement established governing all subsequent steps.
- All material compositions confirmed against reference repo (`nican0r/TRISO-neutronics`, `src/triso/materials.py`).
- All ambiguities surfaced and resolved before step 0 begins.
- No files created; repository layout is step 0's job.

## How it works

This model is a criticality safety (NCS) screening calculation for a laboratory-scale fluidized-bed CVD coater used to coat TRISO fuel kernels. The physical scenario is: HALEU UCO TRISO particles packed at varying densities inside a graphite retort, with CVD process gas filling the interstitial space, surrounded by a graphite heating element. The criticality risk arises from simultaneous fissile mass accumulation, graphite moderation, and the converging cone that tends to concentrate particles. All results are unvalidated until benchmarked against ICSBEP evaluations.

## Design decisions

- **UCO composition confirmed as atom ratios U:C:O = 1.0:0.5:0.4**, 19.75 wt% U-235 enrichment, density 10.5 g/cm³. Source: Demkowicz et al., Nucl. Eng. Des. 329 (2018) 102–111; INL/EXT-10-19476.
- **`c_Graphite` S(α,β) applied to buffer, IPyC, OPyC, and structural graphite.** PyC is turbostratic — `c_Graphite` is an acknowledged approximation, but omitting thermal scattering entirely is non-conservative and a worse approximation. Noted here; do not silently remove this.
- **No `c_Graphite` on UCO kernel carbon.** The carbon in UCO is chemically bound in a ceramic matrix, not graphitic. Free gas scattering kernel applied instead. This errs slightly in the conservative direction for moderation.
- **Cross-section temperature cap at 1200 K.** The ENDF/B-VIII.0 library (`~/Documents/nuclear-data/endfb80_hdf5/cross_sections.xml`) was processed at 900 K and 1200 K only. The furnace lowest operating temperature is 1500 K. OpenMC's `nearest` temperature method will use 1200 K data for all materials when `temperature = 1500 K` is set. This means Doppler broadening is underestimated relative to the actual operating condition. For an NCS evaluation this is non-conservative with respect to Doppler feedback, but since this is a screening model (not a safety limit determination) it is acceptable for Stage 0. This limitation must appear on every reported k-eff.
- **`c_Graphite` thermal scattering capped at 1200 K.** The library contains only a 1200 K evaluation of `c_Graphite`. No 900 K or 1500 K points exist in this build. At 1500 K the code will silently use the 1200 K evaluation; moderation via thermal scattering is thus slightly overestimated (lower temperature → less phonon broadening, more effective thermalization). Non-conservative for moderation. Noted; revisit if a higher-temperature evaluation becomes available.
- **Model temperature parameter = lowest operating temperature (1500 K).** Placed in `params.yaml` as `fuel_temperature_K: 1500  # CONFIRM`. Effective cross-section temperature will be 1200 K due to library limit above.
- **Interstitial gas = SiC coating atmosphere: 98 mol% H₂ + 2 mol% CH₃SiCl₃ (MTS), 101325 Pa, 1873 K (1600 °C).** This is the most hydrogen-rich gas composition in the CVD cycle (the SiC deposition step). It is used throughout the entire model — not just during the SiC step — to bound worst-case moderation: any other CVD atmosphere (e.g. acetylene/propylene for carbon layers) contains less hydrogen and would yield a lower k-eff. Using the H₂-dominant atmosphere everywhere is a deliberate conservative choice for this NCS screening calculation. Pressure confirmed at 1 atm (101325 Pa); note in `params.yaml` with `# CONFIRM` on pressure so it can be revisited if sub-atmospheric operation is ever used.
- **50 cm³ charge = bulk volume** (particles + interstitial void). Determines number of TRISO particles given packing fraction and kernel/coating dimensions.
- **Nozzle simplified** to a single axial opening of equivalent area. The 6 × ø1.5 mm orifices in the ø6 mm throat add geometry complexity with negligible effect on k-eff (small volume, far from fissile region). A `# TODO` marker will indicate the simplification.
- **Vacuum boundary conditions** on all external surfaces. Standard for NCS evaluations.
- **Boron impurity in graphite = zero.** No unverified neutron poison is credited. If boron content is eventually characterized and verified, it can be added as a parametric case.
- **Structural graphite density placeholder = 1.75 g/cm³ (`# CONFIRM`).** From AGR-1 compact matrix (INL/EXT-10-19476 Table 5); may not match the actual furnace graphite grade.
- **All results carry footer: `UNVALIDATED — SCREENING ONLY`** until a validation basis against ICSBEP evaluations is established (to be specified in step 9).

## Open questions carried forward

- Process gas pressure confirmed at 101325 Pa (`# CONFIRM` retained in `params.yaml` for sub-atmospheric operation review).
- Kernel diameter (425 µm AGR-1 nominal), coating thicknesses (buffer 100 µm, IPyC 40 µm, SiC 35 µm, OPyC 40 µm) — all `# CONFIRM`.
- Insulation: graphite felt, placeholder 0.2 g/cm³, thickness unknown — `# CONFIRM`.
- Gas injector coolant channel geometry — `# CONFIRM`.
- Actual furnace graphite grade and density — `# CONFIRM`.
