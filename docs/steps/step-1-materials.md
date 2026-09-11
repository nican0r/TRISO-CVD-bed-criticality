---
step: 1 — Materials
title: All OpenMC material definitions for the CVD furnace NCS model
status: complete
---

## What was implemented

- `furnace/materials.py` — ten material factory functions, each taking the frozen params dict:
  - `uco_kernel(params)` — HALEU UCO fuel kernel (10.5 g/cm³, 19.75 wt% U-235)
  - `buffer_pyc(params)` — porous carbon buffer (1.0 g/cm³, c_Graphite)
  - `ipyc(params)` — inner PyC (1.87 g/cm³, c_Graphite)
  - `sic(params)` — SiC pressure-retention layer (3.20 g/cm³, 1:1 Si:C)
  - `opyc(params)` — outer PyC (1.87 g/cm³, c_Graphite)
  - `graphite_structural(params)` — retort, heater, cone (1.75 g/cm³, c_Graphite)
  - `graphite_felt_insulation(params)` — insulation blanket (0.20 g/cm³, c_Graphite)
  - `process_gas(params)` — H₂/MTS CVD atmosphere; density from ideal gas at 293.6 K, 101325 Pa
  - `water(density)` — light water factory for flooding sweeps (c_H_in_H2O)
  - `air()` — standard dry air for the vented unflooded condition
- `_print_material_table()` and `_verify_u235_atom_density()` helper functions for the diagnostic check block.
- `if __name__ == '__main__':` block that runs the full material table and U-235 density verification.

## How it works

Each factory creates an `openmc.Material`, sets density and composition via `add_element()` with explicit atom-fraction percent type, attaches S(α,β) tables where appropriate, and sets the temperature to 293.6 K. The params dict (a frozen MappingProxyType from `load_params()`) is the sole source of numeric values; no constants are hardcoded in the material functions. The process gas density is computed from the ideal gas law at 293.6 K and the pressure stored in `params['gas']`, consistent with the room-temperature bounding case. The water factory takes density as an argument to support continuous flooding density sweeps without touching params.

## Experimental design

No simulation is run in this step. Materials are constructed and checked analytically only.

## Design decisions

- **All temperatures fixed at 293.6 K.** Per the Step 1 update: the cold, flooded condition is the bounding screening case, so all material cross-sections are evaluated at the ENDF/B-VIII.0 room-temperature grid point (293.6 K ≈ 20 °C). The cold temperature bounds k-eff from above because (1) Doppler broadening of U-238 resonances is minimised at low temperature, reducing resonance capture and increasing k-eff, and (2) thermal scattering in graphite and water is more effective at low temperature, thermalising neutrons more efficiently and increasing the fission cross-section. Both effects push k-eff upward at cold temperature.

- **UCO C and O composition: atom fractions, 0.5 and 0.4 per formula unit.** U(C₀.₅O₀.₄) is the standard AGR-1 UCO stoichiometry expressed as atom ratios (U:C:O = 1:0.5:0.4). Using `percent_type='ao'` with these values correctly represents the crystal chemistry. The alternative (weight fractions) would produce a chemically incorrect composition with too much oxygen relative to carbon.

- **No S(α,β) on the UCO kernel.** The preamble explicitly lists "buffer, IPyC, OPyC, and structural graphite" as the materials requiring c_Graphite — the kernel is not in this list. Physically, the carbon in UCO is a minor constituent embedded in a fuel matrix dominated by uranium; the thermal scattering lattice dynamics of graphite are not applicable to this material form. Attaching c_Graphite to the kernel would be physically unjustified and would overcredit its moderating contribution.

- **No S(α,β) on SiC.** The preamble list omits SiC. The ENDF/B-VIII.0 library contains `c_SiC` for beta-SiC; applying it would be physically defensible since the SiC layer is crystalline. However, the preamble does not credit it and the SiC volume fraction per particle is small. Deferred to a later stage.  `# TODO: evaluate c_SiC impact on k-eff in a later stage.`

- **c_Graphite applied to PyC layers.** PyC (turbostratic carbon) is not graphite, but `c_Graphite` is a better approximation than a free-gas treatment. Omitting thermal scattering from PyC underestimates thermalisation in the coating layers, which errs non-conservatively (lower k-eff than reality). The preamble explicitly requires this approximation and directs that it be noted.

- **What happens to k-eff if c_Graphite is omitted from structural graphite.** Without `c_Graphite`, OpenMC applies free-gas scattering to carbon, which treats graphite atoms as independent scatterers at the material temperature. This misses the crystal phonon modes that make graphite such an effective thermaliser at thermal neutron energies (the 1/v and sub-thermal enhancement). The practical effect is that the free-gas model under-thermalises neutrons, hardening the spectrum and reducing the fission cross-section weighted flux. Published comparisons for graphite-moderated systems show k-eff differences of 2–5% between free-gas and c_Graphite treatments; in a geometry where graphite is the primary moderator, this is not conservative and would give a false subcritical margin.

- **Gas density computed at 293.6 K, 101325 Pa.** The CVD gas composition (98 mol% H₂ + 2 mol% MTS) is held at the SiC-step process conditions as the worst-case hydrogen-rich atmosphere, but the density used in the model is the room-temperature ideal-gas value: ρ ≈ 2.06×10⁻⁴ g/cm³. This is consistent with the 293.6 K cross-section temperature applied to all materials in the room-temperature bounding case.

- **Air composition: N₂, O₂, Ar only (NIST standard dry air); CO₂ neglected.** CO₂ is 0.036 mol% of dry air; its contribution to neutron interaction is negligible. Density computed from ideal gas at 293.6 K, 101325 Pa.

- **U-234 assumption in OpenMC's enrichment helper.** OpenMC warns that at 19.75 wt% enrichment, it applies a fixed U-234/U-235 mass ratio of 0.008 — a correlation valid only for low-enriched uranium. For HALEU from an enrichment cascade, the actual U-234 content depends on the feed and tails assay and is not necessarily 0.8% of U-235 by mass. U-234 is a resonance absorber; the error in k-eff from this approximation is small (U-234 is a minor constituent) but non-zero. To resolve: specify the U-234, U-235, U-236, U-238 isotopic fractions explicitly from the fuel specification when available.

## Assumptions

### Confirmed
- UCO kernel: U(C₀.₅O₀.₄) atom fractions, 19.75 wt% U-235, 10.5 g/cm³. Source: [AGR1]/[INL].
- Buffer: pure C, 1.0 g/cm³; IPyC/OPyC: pure C, 1.87 g/cm³; SiC: 1:1 Si:C, 3.20 g/cm³. Source: [INL] Table 3.
- c_Graphite S(α,β) applied to buffer, IPyC, OPyC, and structural graphite. c_H_in_H2O applied to water.
- Process gas: 98 mol% H₂ + 2 mol% CH₃SiCl₃ (SiC CVD step composition); density at 293.6 K, 101325 Pa (ρ ≈ 2.06×10⁻⁴ g/cm³).
- All material temperatures: 293.6 K (room temperature, cold-bounding case).
- Structural graphite boron impurity: 0 ppm (NCS convention, no unconfirmed poisons credited).

### Defaulted
- SiC has no S(α,β) applied. Could add `c_SiC` but not credited at this stage.
- Air composition from NIST standard dry air (N₂/O₂/Ar); CO₂ neglected.
- U-234 content in HALEU: OpenMC's built-in U-234/U-235 = 0.008 mass ratio approximation used. May be inaccurate for cascade-enriched HALEU.

### Unconfirmed (`# CONFIRM`)
- `graphite_structural_density_gcc` = 1.75 g/cm³ — placeholder; actual furnace graphite grade may differ.
- `graphite_felt_density_gcc` = 0.20 g/cm³ — placeholder; thickness also unknown.
- `kernel_diameter`, all layer thicknesses — AGR-1 nominal values; fabrication drawings may specify different values.
- `gas.pressure_pa` = 101325 Pa — sub-atmospheric operation is possible; lower pressure → lower gas density → marginally lower k-eff (non-conservative direction if pressures are lower than atmospheric).

[AGR1]: Demkowicz et al., Nucl. Eng. Des. 329 (2018) 102-111  
[INL]: INL/EXT-10-19476, AGR-1 fuel specification
