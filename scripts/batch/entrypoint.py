#!/usr/bin/env python
"""In-container dispatcher for AWS Batch array jobs.

Reads a per-array-element parameter row out of a manifest on S3, invokes
`furnace.model.build_model(params, **row)`, runs OpenMC, and uploads the
statepoint plus a small results.json summary back to S3.

Required env (set by the Batch job definition / submitter):
    SWEEP_NAME              — human tag, e.g. step7a_mass
    AWS_BATCH_JOB_ARRAY_INDEX — array child index (Batch injects)
    PARAM_MANIFEST_S3       — s3://bucket/prefix/manifest.json
    RESULTS_PREFIX_S3       — s3://bucket/runs/<sweep>/<submission_id>/

Optional env:
    OPENMC_CROSS_SECTIONS   — path to cross_sections.xml (defaults to /xs/...)
    SMOKE                   — if set, force a tiny 200-particle run and skip S3.

The manifest is a JSON list of dicts. Every dict is forwarded to
build_model as **kwargs; unknown keys will raise, which is the point —
keeping the container sweep-agnostic pushes schema mistakes to submit time.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path
from urllib.parse import urlparse

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import openmc  # noqa: E402

from furnace.params import load_params  # noqa: E402
from furnace.model import build_model, export_and_run  # noqa: E402


def _s3_parse(uri: str) -> tuple[str, str]:
    p = urlparse(uri)
    if p.scheme != "s3":
        raise ValueError(f"not an s3 uri: {uri}")
    return p.netloc, p.path.lstrip("/")


def _s3_download(uri: str, dest: Path) -> None:
    subprocess.run(["aws", "s3", "cp", uri, str(dest)], check=True)


def _s3_upload_dir(local_dir: Path, uri_prefix: str) -> None:
    subprocess.run(
        ["aws", "s3", "cp", "--recursive", str(local_dir), uri_prefix],
        check=True,
    )


def _smoke_row() -> dict:
    # charge_mass_g=5.0 gives ~8500 TRISO particles (vs 162K at default 95g),
    # making local geometry+XS setup tractable while still banking fission sites.
    return {
        "state": "fluidized",
        "stage": "bare_kernel",
        "charge_mass_g": 5.0,
        "n_particles": 500,
        "n_inactive": 5,
        "n_active": 15,
        "seed": 42,
    }


def _load_row() -> tuple[dict, int]:
    if os.environ.get("SMOKE"):
        return _smoke_row(), 0

    manifest_uri = os.environ["PARAM_MANIFEST_S3"]
    # AWS_BATCH_JOB_ARRAY_INDEX is only set for array jobs; a 1-row
    # single-job submission (e.g. the smoke test) implies index 0.
    idx = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX", "0"))

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        _s3_download(manifest_uri, Path(f.name))
        manifest = json.loads(Path(f.name).read_text())

    if not isinstance(manifest, list):
        raise ValueError(f"manifest must be a JSON list, got {type(manifest)}")
    if idx >= len(manifest):
        raise IndexError(f"array index {idx} out of range for manifest of size {len(manifest)}")

    return manifest[idx], idx


def _summarize_statepoint(sp_path: Path, row: dict, stats, wall_seconds: float) -> dict:
    with openmc.StatePoint(str(sp_path)) as sp:
        k = sp.keff
        return {
            "row": row,
            "u235_mass_g": stats.u235_mass_g,
            "pf_achieved": stats.pf_achieved,
            "bed_height_cm": stats.bed_height_cm,
            "V_bulk_cm3": stats.V_bulk_cm3,
            "n_particles_in_bed": stats.n_particles,
            "keff": float(k.nominal_value),
            "keff_std": float(k.std_dev),
            "keff_plus_2sigma": float(k.nominal_value + 2 * k.std_dev),
            "wall_seconds": wall_seconds,
            "sweep_name": os.environ.get("SWEEP_NAME", "unnamed"),
            "git_sha": os.environ.get("GIT_SHA", ""),
        }


def main() -> int:
    row, idx = _load_row()
    print(f"[entrypoint] array index {idx} row: {json.dumps(row)}", flush=True)

    params = load_params()

    t0 = time.time()
    model, stats = build_model(params, **row)

    scratch = Path(tempfile.mkdtemp(prefix=f"run_{idx}_"))
    try:
        sp_path = export_and_run(model, scratch)
        wall = time.time() - t0

        summary = _summarize_statepoint(sp_path, row, stats, wall)
        (scratch / "results.json").write_text(json.dumps(summary, indent=2))
        print(f"[entrypoint] k-eff = {summary['keff']:.5f} ± {summary['keff_std']:.5f} "
              f"({wall:.1f}s wall)", flush=True)

        # Keep only the small artifacts by default; statepoint tends to be a
        # few MB but geometry.xml can be huge for TRISO lattices, so skip it.
        for junk in ("geometry.xml",):
            p = scratch / junk
            if p.exists():
                p.unlink()

        if os.environ.get("SMOKE"):
            print(f"[entrypoint] SMOKE mode, skipping upload; artifacts in {scratch}",
                  flush=True)
            return 0

        prefix = os.environ["RESULTS_PREFIX_S3"].rstrip("/") + f"/{idx:04d}/"
        _s3_upload_dir(scratch, prefix)
        print(f"[entrypoint] uploaded → {prefix}", flush=True)
        return 0
    except Exception:
        traceback.print_exc()
        return 1
    finally:
        if not os.environ.get("SMOKE"):
            shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
