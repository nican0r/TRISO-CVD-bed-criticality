"""TRISO layer and structural material definitions for the CVD furnace NCS model.

All materials are set to 293.6 K (ENDF room-temperature evaluation point).
The cold, flooded condition is the bounding screening case for this system.

Sources:
  [AGR1] Demkowicz et al., Nucl. Eng. Des. 329 (2018) 102-111
  [INL]  INL/EXT-10-19476, AGR-1 fuel specification
"""

from __future__ import annotations

import types

import openmc

# ENDF/B-VIII.0 room-temperature evaluation; nearest available point to 20 °C.
_ROOM_TEMP_K = 293.6

# SI gas constant (J/mol/K)
_R = 8.314462618


def uco_kernel(params: types.MappingProxyType) -> openmc.Material:
    """UCO kernel: U(C0.5 O0.4) by atom fraction, 19.75 wt% U-235 enrichment.

    Composition fixed per [AGR1]/[INL]: atom fractions 1:0.5:0.4 (U:C:O).
    No c_Graphite S(α,β): carbon is a minor constituent in a fuel matrix
    dominated by uranium interaction cross-sections.
    """
    m = openmc.Material(name='uco_kernel')
    m.set_density('g/cm3', params['materials']['kernel_density_gcc'])
    enrich = params['materials']['enrichment_wt_pct']
    m.add_element('U', 1.0, percent_type='ao', enrichment=enrich, enrichment_type='wo')
    m.add_element('C', 0.5, percent_type='ao')
    m.add_element('O', 0.4, percent_type='ao')
    m.temperature = _ROOM_TEMP_K
    return m


def buffer_pyc(params: types.MappingProxyType) -> openmc.Material:
    """Porous carbon buffer layer. [INL] Table 3, 1.0 g/cm3."""
    m = openmc.Material(name='buffer_pyc')
    m.set_density('g/cm3', params['materials']['buffer_density_gcc'])
    m.add_element('C', 1.0, percent_type='ao')
    m.add_s_alpha_beta('c_Graphite')
    m.temperature = _ROOM_TEMP_K
    return m


def ipyc(params: types.MappingProxyType) -> openmc.Material:
    """Inner pyrolytic carbon layer. [INL] Table 3, 1.87 g/cm3.

    PyC is turbostratic; c_Graphite is an approximation but omitting thermal
    scattering entirely errs non-conservatively and is a worse approximation.
    """
    m = openmc.Material(name='ipyc')
    m.set_density('g/cm3', params['materials']['ipyc_density_gcc'])
    m.add_element('C', 1.0, percent_type='ao')
    m.add_s_alpha_beta('c_Graphite')
    m.temperature = _ROOM_TEMP_K
    return m


def sic(params: types.MappingProxyType) -> openmc.Material:
    """SiC pressure-retention layer: stoichiometric 1:1 Si:C. [INL] Table 3, 3.20 g/cm3.

    No S(α,β) applied: the preamble explicit material list (buffer/IPyC/OPyC/graphite)
    excludes SiC, and c_SiC is not credited here pending review of which
    thermal scattering table is appropriate for this CVD-deposited beta-SiC.
    # TODO: evaluate c_SiC from ENDF/B-VIII.0 for a later stage.
    """
    m = openmc.Material(name='sic')
    m.set_density('g/cm3', params['materials']['sic_density_gcc'])
    m.add_element('Si', 1.0, percent_type='ao')
    m.add_element('C', 1.0, percent_type='ao')
    m.temperature = _ROOM_TEMP_K
    return m


def opyc(params: types.MappingProxyType) -> openmc.Material:
    """Outer pyrolytic carbon layer. [INL] Table 3, 1.87 g/cm3. Same as IPyC."""
    m = openmc.Material(name='opyc')
    m.set_density('g/cm3', params['materials']['opyc_density_gcc'])
    m.add_element('C', 1.0, percent_type='ao')
    m.add_s_alpha_beta('c_Graphite')
    m.temperature = _ROOM_TEMP_K
    return m


def graphite_structural(params: types.MappingProxyType) -> openmc.Material:
    """Nuclear-grade structural graphite: retort, heater element, cone.

    Boron content = 0 ppm per NCS convention; no unconfirmed neutron poisons credited.
    Density is a placeholder — see params.yaml # CONFIRM.

    Without c_Graphite, OpenMC treats carbon scattering as free-gas, which
    underestimates thermal thermalisation in the Bragg regime, reducing the
    computed neutron moderation and lowering k-eff by several percent.
    Omitting it errs non-conservatively.
    """
    m = openmc.Material(name='graphite_structural')
    m.set_density('g/cm3', params['materials']['graphite_structural_density_gcc'])
    m.add_element('C', 1.0, percent_type='ao')
    m.add_s_alpha_beta('c_Graphite')
    m.temperature = _ROOM_TEMP_K
    return m


def process_gas(params: types.MappingProxyType) -> openmc.Material:
    """CVD process gas: H2 (98 mol%) + MTS/CH3SiCl3 (2 mol%).

    Only hydrogen-bearing species are modelled. Inert carrier gas (Ar) is absent
    from params mole_fractions and treated as void — it contributes nothing to
    moderation at these densities.

    Density is computed from the ideal gas law at 293.6 K (room temperature),
    101325 Pa, consistent with the room-temperature bounding case used for all
    material cross-section temperatures. This gives ρ ≈ 2.06×10⁻⁴ g/cm³.

    MTS = CH3SiCl3: M = 12.011 + 3*1.008 + 28.086 + 3*35.453 = 149.480 g/mol
    """
    g = params['gas']
    x_H2 = g['mole_fractions']['H2']
    x_MTS = g['mole_fractions']['CH3SiCl3']

    M_H2 = 2 * 1.008                                           # g/mol
    M_MTS = 12.011 + 3 * 1.008 + 28.086 + 3 * 35.453          # g/mol  (149.480)
    M_avg = x_H2 * M_H2 + x_MTS * M_MTS

    # Ideal gas: ρ [g/m3] = P [Pa] * M [g/mol] / (R [J/mol/K] * T [K]); ×1e-6 → g/cm3
    P = g['pressure_pa']
    rho_gcc = P * M_avg / (_R * _ROOM_TEMP_K) * 1e-6

    m = openmc.Material(name='process_gas')
    m.set_density('g/cm3', rho_gcc)

    # Atom counts per mol of gas mixture (unnormalised; OpenMC normalises internally)
    # H2 contribution: x_H2 mol × 2 H/mol
    # MTS contribution: x_MTS mol × (1C + 3H + 1Si + 3Cl)/mol
    m.add_element('H',  x_H2 * 2 + x_MTS * 3, percent_type='ao')
    m.add_element('C',  x_MTS * 1,             percent_type='ao')
    m.add_element('Si', x_MTS * 1,             percent_type='ao')
    m.add_element('Cl', x_MTS * 3,             percent_type='ao')

    m.temperature = _ROOM_TEMP_K
    return m


def water(density: float) -> openmc.Material:
    """Light water factory for flooding sweep cases.

    Args:
        density: mass density in g/cm3. Pass 1.0 for ambient liquid water;
                 lower values for partial flooding or steam.
    """
    m = openmc.Material(name=f'water_{density:.4f}gcc')
    m.set_density('g/cm3', density)
    m.add_element('H', 2.0, percent_type='ao')
    m.add_element('O', 1.0, percent_type='ao')
    # c_H_in_H2O S(α,β) omitted: endfb80_hdf5 c_H_in_H2O.h5 has empty kTs (no temperature data).
    # Free-gas treatment is non-conservative (less thermal absorption) but water is only in the
    # injector cooling annulus — not in the neutron path — so the k-eff impact is negligible for
    # Stage 0 NCS.
    m.temperature = _ROOM_TEMP_K
    return m


def air() -> openmc.Material:
    """Standard dry air for the vented, unflooded condition.

    Composition: NIST standard dry air (mol%): N2 78.084, O2 20.946, Ar 0.934.
    CO2 (0.036 mol%) is neglected. Density from ideal gas at 293.6 K, 101325 Pa.

    M_air = 0.78084*28.014 + 0.20946*31.998 + 0.00934*39.948 = 28.949 g/mol
    rho   = 101325 * 28.949e-3 / (8.314 * 293.6) * 1e-3 = 1.202e-3 g/cm3
    """
    M_N2 = 28.014
    M_O2 = 31.998
    M_Ar = 39.948
    x_N2, x_O2, x_Ar = 0.78084, 0.20946, 0.00934
    M_avg = x_N2 * M_N2 + x_O2 * M_O2 + x_Ar * M_Ar

    rho_gcc = 101325.0 * M_avg / (_R * _ROOM_TEMP_K) * 1e-6

    m = openmc.Material(name='air')
    m.set_density('g/cm3', rho_gcc)
    # atom fractions: each N2 contributes 2 N, each O2 contributes 2 O, Ar is monatomic
    m.add_element('N',  x_N2 * 2, percent_type='ao')
    m.add_element('O',  x_O2 * 2, percent_type='ao')
    m.add_element('Ar', x_Ar,     percent_type='ao')
    m.temperature = _ROOM_TEMP_K
    return m


# ---------------------------------------------------------------------------
# Diagnostic check
# ---------------------------------------------------------------------------

def _sab_names(mat: openmc.Material) -> str:
    if not mat._sab:
        return 'none'
    return ', '.join(entry[0] for entry in mat._sab)


def _print_material_table(materials: list[openmc.Material]) -> None:
    header = f"{'Name':<28} {'Density (g/cm3)':>16} {'Temp (K)':>9} {'S(a,b)':>22}"
    print('\n' + header)
    print('-' * len(header))
    for mat in materials:
        sab = _sab_names(mat)
        print(f"{mat.name:<28} {mat.density:>16.4e} {mat.temperature:>9.1f} {sab:>22}")
        dens = mat.get_nuclide_atom_densities()
        for nuc, nd in sorted(dens.items()):
            print(f"    {nuc:<12} {nd:.4e} atoms/b-cm")


def _verify_u235_atom_density(params: types.MappingProxyType) -> None:
    """Print analytic vs OpenMC U-235 atom density for the UCO kernel."""
    import math

    N_A = 6.02214076e23
    M_U235 = 235.044
    M_U238 = 238.051
    M_C = 12.011
    M_O = 15.999

    w = params['materials']['enrichment_wt_pct'] / 100.0   # weight fraction U-235
    rho = params['materials']['kernel_density_gcc']         # g/cm3

    # Atom fraction of U-235 from weight percent enrichment
    x235 = (w / M_U235) / (w / M_U235 + (1.0 - w) / M_U238)
    M_U = x235 * M_U235 + (1.0 - x235) * M_U238

    # Molar mass of one formula unit U(C0.5)(O0.4)
    M_formula = M_U + 0.5 * M_C + 0.4 * M_O

    # Number density of formula units → atom density of U-235
    n_formula = rho * N_A / M_formula          # formula units/cm3
    N_U235_analytic = n_formula * x235         # atoms/cm3

    # OpenMC value from the material object (atoms/b-cm × 1e24 = atoms/cm3)
    kernel = uco_kernel(params)
    nuc_dens = kernel.get_nuclide_atom_densities()
    N_U235_openmc = nuc_dens.get('U235', 0.0) * 1e24

    print('\nU-235 atom density verification (UCO kernel):')
    print(f'  Enrichment (wt%)        : {w*100:.2f}')
    print(f'  Atom fraction U-235     : {x235:.5f}')
    print(f'  M(formula unit) g/mol   : {M_formula:.4f}')
    print(f'  Analytic  N(U-235)      : {N_U235_analytic:.4e}  atoms/cm3')
    print(f'  OpenMC    N(U-235)      : {N_U235_openmc:.4e}  atoms/cm3')
    print(f'  Ratio (OpenMC/analytic) : {N_U235_openmc / N_U235_analytic:.6f}')


if __name__ == '__main__':
    from furnace.params import check_env, load_params

    check_env()
    params = load_params()

    mats = [
        uco_kernel(params),
        buffer_pyc(params),
        ipyc(params),
        sic(params),
        opyc(params),
        graphite_structural(params),
        process_gas(params),
        water(1.0),
        air(),
    ]

    _print_material_table(mats)
    _verify_u235_atom_density(params)
