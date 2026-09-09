"""Parameter loading, unit conversion, and environment validation.

params.yaml stores dimensions in user-facing units (mm for retort/cone/bed
geometry, µm for TRISO layers).  load_params() converts both to cm and
inserts *_cm keys so downstream OpenMC code never needs the conversion factors.

Chosen return type: recursive types.MappingProxyType (frozen dict).
Rationale: the schema grows across steps; nested dataclasses would require
updating field definitions every time a new parameter is added.
MappingProxyType applied recursively gives the same mutation protection at
every nesting level with zero per-field boilerplate.
"""

from __future__ import annotations

import os
import types
from pathlib import Path
from typing import Any

import yaml

_PARAMS_PATH = Path(__file__).parent.parent / "params.yaml"
_MM_TO_CM = 0.1
_UM_TO_CM = 1e-4

# Only keys whose names contain one of these substrings get a _cm companion.
# Prevents volume (_cc), angle (_deg), pressure (_pa), and temperature (_k)
# keys from receiving a nonsensical cm conversion.
_LENGTH_SUFFIXES = (
    "diameter", "depth", "height", "length", "drop", "thickness", "id", "od", "throat",
)

_REQUIRED_TOP_KEYS = frozenset({"dimensions", "triso", "materials", "gas", "model"})
_REQUIRED_MATERIAL_KEYS = frozenset({
    "enrichment_wt_pct",
    "kernel_density_gcc",
    "buffer_density_gcc",
    "ipyc_density_gcc",
    "sic_density_gcc",
    "opyc_density_gcc",
    "graphite_structural_density_gcc",
})


def _is_length_key(key: str) -> bool:
    k = key.lower()
    return any(s in k for s in _LENGTH_SUFFIXES)


def _add_cm_keys(d: dict, factor: float) -> dict:
    """Return a shallow copy of d with <key>_cm added for every length value."""
    out: dict[str, Any] = {}
    for k, v in d.items():
        out[k] = v
        if isinstance(v, (int, float)) and _is_length_key(k):
            out[f"{k}_cm"] = round(v * factor, 10)
    return out


def _convert_dimensions(dims: dict) -> dict:
    result: dict[str, Any] = {"units": "mm"}
    for name, section in dims.items():
        if name == "units":
            continue
        result[name] = _add_cm_keys(section, _MM_TO_CM) if isinstance(section, dict) else section
    return result


def _convert_triso(triso: dict) -> dict:
    result: dict[str, Any] = {"units": "um"}
    for k, v in triso.items():
        if k == "units":
            continue
        result[k] = v
        if isinstance(v, (int, float)) and _is_length_key(k):
            result[f"{k}_cm"] = round(v * _UM_TO_CM, 10)
    return result


def _freeze(obj: Any) -> Any:
    """Recursively wrap dicts in MappingProxyType and lists in tuples."""
    if isinstance(obj, dict):
        return types.MappingProxyType({k: _freeze(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return tuple(_freeze(item) for item in obj)
    return obj


def _validate(raw: dict) -> None:
    missing_top = _REQUIRED_TOP_KEYS - raw.keys()
    if missing_top:
        raise KeyError(f"params.yaml missing top-level keys: {sorted(missing_top)}")
    missing_mat = _REQUIRED_MATERIAL_KEYS - raw["materials"].keys()
    if missing_mat:
        raise KeyError(f"params.yaml missing materials keys: {sorted(missing_mat)}")
    if raw.get("model", {}).get("seed") is None:
        raise KeyError("params.yaml: model.seed must be set for reproducibility")


def load_params(path: Path | str | None = None) -> types.MappingProxyType:
    """Read params.yaml, validate, add _cm keys, and return a frozen mapping.

    Dimensional conversions applied:
        dimensions.*.<key>_cm  = <key> × 0.1     (mm → cm)
        triso.<key>_cm         = <key> × 1e-4    (µm → cm)

    The returned object is read-only at every nesting level.
    """
    p = Path(path) if path is not None else _PARAMS_PATH
    with p.open() as fh:
        raw = yaml.safe_load(fh)

    _validate(raw)

    converted = dict(raw)
    converted["dimensions"] = _convert_dimensions(raw["dimensions"])
    converted["triso"] = _convert_triso(raw["triso"])
    return _freeze(converted)


def check_env() -> None:
    """Verify OPENMC_CROSS_SECTIONS is set, the file exists, and print library info.

    Raises EnvironmentError with a diagnostic message rather than letting the
    problem surface later as a cryptic OpenMC error inside a long run.

    openmc is imported here (not at module level) so furnace.params is
    importable in environments without OpenMC installed, e.g. for unit tests
    or parameter inspection.
    """
    xs = os.environ.get("OPENMC_CROSS_SECTIONS")
    if not xs:
        raise EnvironmentError(
            "OPENMC_CROSS_SECTIONS is not set.\n"
            "Expected: ~/Documents/nuclear-data/endfb80_hdf5/cross_sections.xml\n"
            "Add it to ~/.zshrc and restart your shell."
        )
    xs_path = Path(xs).expanduser().resolve()
    if not xs_path.exists():
        raise EnvironmentError(
            f"Cross-sections file not found: {xs_path}\n"
            "Check that OPENMC_CROSS_SECTIONS points to an existing file."
        )

    import openmc  # noqa: PLC0415 — deferred import, see docstring

    lib = openmc.data.DataLibrary.from_xml(str(xs_path))
    print(f"OpenMC version  : {openmc.__version__}")
    print(f"Library path    : {xs_path}")
    print(f"Nuclides/tables : {len(lib)}")
