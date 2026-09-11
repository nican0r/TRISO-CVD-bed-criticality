#!/usr/bin/env python3
"""Add 294 K room-temperature data to endfb80_hdf5/neutron/ for project nuclides.

The project sets all materials to 293.6 K (room-temperature bounding case).
The endfb80_hdf5 library only has data at 900 K / 1200 K / 2500 K, so OpenMC's
'nearest' method currently snaps to 900 K — non-conservative for criticality safety
(900 K data has more Doppler broadening → lower resonance self-shielding → lower k-eff
relative to the true room-temperature cross sections).

This script copies 294 K groups from nndc_hdf5/ into endfb80_hdf5/neutron/ for all
project-relevant nuclides present in both libraries.  After running it, OpenMC will
snap 293.6 K → 294 K, which is the correct room-temperature evaluation point.

Both libraries are ENDF/B-VIII.0, format version [3, 0].  Merging is safe because:
  - atomic_weight_ratio is identical in both for every nuclide checked.
  - Each temperature point (294 K, 900 K, …) uses its own independent energy grid
    stored under energy/<temp>/, so grids from different temperatures never conflict.

Known remaining limitations after running this script:
  C12, C13  — nndc_hdf5 uses C0 (natural element), not individual isotopes.
              Carbon cross sections will still snap to 900 K.
  c_Graphite — thermal S(α,β) in endfb80_hdf5/thermal/ only has 1200 K data.
              Thermal scattering for graphite materials will still use 1200 K.
  c_H_in_H2O — empty kTs in endfb80_hdf5; water S(α,β) already omitted in model.

Usage (only requires h5py, not OpenMC):
    python scripts/add_room_temp.py           # runs with per-file .bak backups
    python scripts/add_room_temp.py --dry-run # lists what would change
    python scripts/add_room_temp.py --verify  # confirms 294 K is present after merge
    python scripts/add_room_temp.py --no-backup  # skip backups (faster)
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

import h5py

NUCLEAR_DATA = Path.home() / "Documents/nuclear-data"
ENDF_NEUTRON = NUCLEAR_DATA / "endfb80_hdf5" / "neutron"
NNDC_DIR     = NUCLEAR_DATA / "nndc_hdf5"
TEMP_KEY     = "294K"

# Elements used in furnace/materials.py via add_element().
# OpenMC expands each element to all isotopes present in the cross_sections.xml.
PROJECT_ELEMENTS = {"U", "C", "O", "Si", "H", "N", "Cl", "Ar"}


def element_of(nuclide: str) -> str:
    return re.match(r"([A-Za-z]+)", nuclide).group(1)


def copy_temp_into(src: h5py.File, dst: h5py.File, nuclide: str) -> dict[str, str]:
    """Copy all TEMP_KEY subgroups from nndc src into endfb80 dst for one nuclide."""
    result: dict[str, str] = {}
    src_nuc = src[nuclide]
    dst_nuc = dst[nuclide]

    # kTs scalar (kT value in eV for this temperature)
    if "kTs" in src_nuc and TEMP_KEY in src_nuc["kTs"]:
        if "kTs" not in dst_nuc:
            dst_nuc.create_group("kTs")
        if TEMP_KEY not in dst_nuc["kTs"]:
            src.copy(f"{nuclide}/kTs/{TEMP_KEY}", dst_nuc["kTs"], name=TEMP_KEY)
            result["kTs"] = "copied"
        else:
            result["kTs"] = "already-exists"

    # Energy grid for this temperature (each temperature has its own grid)
    if "energy" in src_nuc and TEMP_KEY in src_nuc["energy"]:
        if "energy" not in dst_nuc:
            dst_nuc.create_group("energy")
        if TEMP_KEY not in dst_nuc["energy"]:
            src.copy(f"{nuclide}/energy/{TEMP_KEY}", dst_nuc["energy"], name=TEMP_KEY)
            result["energy"] = "copied"
        else:
            result["energy"] = "already-exists"

    # Reaction cross sections (MT numbers stored as reaction_NNN groups)
    n_copied = n_absent = n_exists = 0
    if "reactions" in src_nuc:
        dst_rxn_keys = set(dst_nuc.get("reactions", {}).keys())
        for rxn in src_nuc["reactions"]:
            if rxn not in dst_rxn_keys:
                # Reaction exists in nndc but not endfb80 — skip (different processing)
                n_absent += 1
                continue
            src_rxn = src_nuc[f"reactions/{rxn}"]
            dst_rxn = dst_nuc[f"reactions/{rxn}"]
            if TEMP_KEY in src_rxn:
                if TEMP_KEY not in dst_rxn:
                    src.copy(f"{nuclide}/reactions/{rxn}/{TEMP_KEY}", dst_rxn, name=TEMP_KEY)
                    n_copied += 1
                else:
                    n_exists += 1
    result["reactions"] = f"{n_copied} copied, {n_absent} absent-in-endf, {n_exists} already-exist"

    # Unresolved resonance probability tables
    if "urr" in src_nuc and TEMP_KEY in src_nuc["urr"]:
        if "urr" not in dst_nuc:
            dst_nuc.create_group("urr")
        if TEMP_KEY not in dst_nuc["urr"]:
            src.copy(f"{nuclide}/urr/{TEMP_KEY}", dst_nuc["urr"], name=TEMP_KEY)
            result["urr"] = "copied"
        else:
            result["urr"] = "already-exists"

    return result


def discover_candidates() -> tuple[list[str], list[str]]:
    """Return (to_process, missing_from_nndc) for project-element nuclides."""
    nndc_nucs = {f.stem for f in NNDC_DIR.glob("*.h5")}
    endf_nucs = {f.stem for f in ENDF_NEUTRON.glob("*.h5")}
    project_endf = {n for n in endf_nucs if element_of(n) in PROJECT_ELEMENTS}

    to_process      = sorted(n for n in project_endf if n in nndc_nucs)
    missing_in_nndc = sorted(n for n in project_endf if n not in nndc_nucs)
    return to_process, missing_in_nndc


def verify() -> bool:
    """Check that 294 K was successfully added to all processable nuclides."""
    to_process, missing = discover_candidates()
    ok = True
    for nuc in to_process:
        dst = ENDF_NEUTRON / f"{nuc}.h5"
        with h5py.File(dst, "r") as f:
            temps = list(f[nuc].get("kTs", {}).keys())
            if TEMP_KEY in temps:
                print(f"  {nuc}: OK  (kTs = {temps})")
            else:
                print(f"  {nuc}: MISSING 294 K  (kTs = {temps})", file=sys.stderr)
                ok = False
    if missing:
        print(f"\nNuclides without 294 K (expected, no nndc source): {', '.join(missing)}")
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--dry-run",   action="store_true", help="Show plan without modifying files")
    parser.add_argument("--verify",    action="store_true", help="Check 294 K is present; exit 1 if any missing")
    parser.add_argument("--no-backup", action="store_true", help="Skip .bak copies of modified files")
    args = parser.parse_args()

    if args.verify:
        sys.exit(0 if verify() else 1)

    to_process, missing_in_nndc = discover_candidates()

    print(f"Nuclides to update ({len(to_process)}): {', '.join(to_process)}")
    if missing_in_nndc:
        print(f"\nNuclides NOT in nndc_hdf5 (will remain at 900 K):")
        print(f"  {', '.join(missing_in_nndc)}")
        print(f"  → C12/C13: nndc uses C0 (natural element), no isotopic match")

    if args.dry_run:
        print("\nDry run complete — no files modified.")
        return

    print()
    modified = []
    for nuc in to_process:
        src_path = NNDC_DIR / f"{nuc}.h5"
        dst_path = ENDF_NEUTRON / f"{nuc}.h5"

        # Skip if 294 K already present in target
        with h5py.File(dst_path, "r") as f:
            if TEMP_KEY in f[nuc].get("kTs", {}):
                print(f"  {nuc}: 294 K already present — skipped")
                continue

        if not args.no_backup:
            bak = dst_path.with_suffix(".h5.bak")
            if not bak.exists():
                shutil.copy2(dst_path, bak)

        with h5py.File(src_path, "r") as src, h5py.File(dst_path, "r+") as dst:
            res = copy_temp_into(src, dst, nuc)

        modified.append(nuc)
        print(f"  {nuc}: {res}")

    print(f"\n{'='*60}")
    print(f"Modified {len(modified)} / {len(to_process)} files.")
    print(f"\nRemaining 900 K nuclides: {', '.join(missing_in_nndc) or 'none'}")
    print(f"Thermal scattering (c_Graphite etc.): still at 1200 K — no nndc source")
    print(f"\nRun --verify to confirm, then re-upload to AWS:")
    print(f"  python scripts/add_room_temp.py --verify")
    print(f"  ./aws/hydrate-nuclear-data.sh")


if __name__ == "__main__":
    main()
