#!/usr/bin/env python
"""Warm the packing cache for the two most-reused configurations.

Covers:
  - fluidized 95 g  : step5 fluidized case (1 job)
  - collapsed 95 g  : step5 collapsed + ALL step8 cases (25 jobs)

step7a large-mass cases (142.5 g – 1900 g) are NOT pre-baked here because
local generation would take 3–40 h per case.  Those Batch jobs generate their
own cache during the run; the 24-hour job timeout (revision 3) provides
sufficient headroom.

Only the packing geometry matters for cache keys, so flood_extent / background
/ particle counts are irrelevant — we pick the cheapest valid combination.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from furnace.params import load_params
from furnace.model import build_model

CASES = [
    # (state, charge_mass_g, seed)
    # None → use params default (95 g)
    ("fluidized", None,  42),   # step5 nominal fluidized
    ("collapsed", None,  42),   # step5 nominal collapsed + ALL step8 cases + step7a 1×
    # Large masses omitted — Batch jobs pack these themselves under 24-hour timeout.
]


def main() -> None:
    params = load_params()

    cache_dir = _REPO_ROOT / "cases" / ".triso_cache"
    files_before = set(cache_dir.glob("*.npz")) if cache_dir.exists() else set()

    for state, mass, seed in CASES:
        label = f"state={state}, mass={mass or 'default'}, seed={seed}"
        print(f"\n{'='*60}")
        print(f"Packing: {label}")
        t0 = time.time()
        build_model(
            params,
            state=state,
            stage="bare_kernel",
            background="gas",
            charge_mass_g=mass,
            n_particles=100,   # irrelevant — packing happens before simulation
            n_inactive=5,
            n_active=10,
            seed=seed,
        )
        elapsed = time.time() - t0
        files_after = set(cache_dir.glob("*.npz")) if cache_dir.exists() else set()
        new_files = files_after - files_before
        print(f"Done in {elapsed:.1f}s — {len(new_files)} new cache file(s) written")
        files_before = files_after

    total = len(set(cache_dir.glob("*.npz"))) if cache_dir.exists() else 0
    print(f"\nCache warmup complete. Total cache files: {total}")


if __name__ == "__main__":
    main()
