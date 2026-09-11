#!/usr/bin/env python
"""Parallel packing-cache warmup using multiprocessing.

Packs each slab of the two critical configurations in parallel:
  - fluidized 95 g  (pf ≈ 0.333)
  - collapsed 95 g  (pf = 0.50)

Each slab pack is independent, so we spawn one worker per slab and let
the OS distribute across all available CPUs.

Usage:
    python scripts/warm_packing_cache_parallel.py [-j N]

    -j N  : max worker processes (default: all logical CPUs minus 1)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from furnace.params import load_params
from furnace.geometry import staircase_bed
from furnace.triso import _stage_outer_radius


def _cache_path(cache_dir: Path, region_ll, region_ur, pf: float,
                r: float, seed: int) -> Path:
    ll = [round(float(v), 8) for v in region_ll]
    ur = [round(float(v), 8) for v in region_ur]
    payload = json.dumps(
        {'ll': ll, 'ur': ur, 'pf': round(pf, 8), 'r': round(r, 8), 'seed': seed},
        sort_keys=True,
    )
    key = hashlib.sha256(payload.encode()).hexdigest()[:16]
    return cache_dir / f"{key}.npz"


def pack_slab_worker(args: tuple) -> tuple[str, int, float, str | None]:
    """Worker: pack one slab and save to cache.

    Returns (label, n_particles, elapsed_seconds, error_or_None).
    """
    (z_bot, z_top, r_slab, pf, r_particle, seed,
     cache_dir_str, slab_idx, state_label) = args

    cache_dir = Path(cache_dir_str)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Build bounding box for cache key (must match _make_cache_key in triso.py)
    region_ll = (-r_slab, -r_slab, z_bot)
    region_ur = (r_slab,  r_slab,  z_top)

    cp = _cache_path(cache_dir, region_ll, region_ur, pf, r_particle, seed)

    label = f"{state_label}/slab_{slab_idx:02d}"
    if cp.exists():
        centers = np.load(cp)['centers']
        return label, len(centers), 0.0, None

    # Estimate particle count; guard against empty slabs.
    V_sphere = (4.0 / 3.0) * math.pi * r_particle ** 3
    bb_vol = (2 * r_slab) ** 2 * (z_top - z_bot)   # bounding box, not slab
    n_est = int(pf * bb_vol / V_sphere)
    if n_est == 0:
        np.savez(cp, centers=np.empty((0, 3), dtype=float))
        return label, 0, 0.0, None

    import openmc  # import inside worker to avoid fork issues
    import openmc.model

    # Thin-slab / narrow-cylinder guard: mirrors the logic in triso.py pack_bed.
    _TRIGGER_DIAMETERS = 2.5
    _NARROW_R_RADII = 15.0
    _MIN_H_DIAMETERS = 5.0

    z_extent = z_top - z_bot
    trigger_h = _TRIGGER_DIAMETERS * 2.0 * r_particle
    narrow_cylinder = r_slab / r_particle < _NARROW_R_RADII
    min_h = _MIN_H_DIAMETERS * 2.0 * r_particle

    t0 = time.time()
    try:
        if (z_extent < trigger_h or narrow_cylinder) and min_h > z_extent:
            z_top_ext = z_bot + min_h
            ext_region = (
                -openmc.ZCylinder(r=r_slab)
                & +openmc.ZPlane(z0=z_bot)
                & -openmc.ZPlane(z0=z_top_ext)
            )
            pf_ext = pf * z_extent / min_h
            all_centers = openmc.model.pack_spheres(
                radius=r_particle,
                region=ext_region,
                pf=pf_ext,
                seed=seed,
            )
            mask = (all_centers[:, 2] - r_particle >= z_bot) & \
                   (all_centers[:, 2] + r_particle <= z_top)
            centers = all_centers[mask]
        else:
            cyl = openmc.ZCylinder(r=r_slab)
            zp_bot = openmc.ZPlane(z0=z_bot)
            zp_top = openmc.ZPlane(z0=z_top)
            region = -cyl & +zp_bot & -zp_top
            centers = openmc.model.pack_spheres(
                radius=r_particle,
                region=region,
                pf=pf,
                seed=seed,
            )
    except Exception as e:
        return label, 0, time.time() - t0, str(e)

    np.savez(cp, centers=centers)
    return label, len(centers), time.time() - t0, None


def build_slab_tasks(params, state: str, seed: int,
                     charge_mass_g: float | None = None) -> list[tuple]:
    """Generate the list of worker args for every non-trivially-empty slab."""
    dim = params['dimensions']
    mdl = params['model']
    stage = 'bare_kernel'

    r_particle = _stage_outer_radius(stage, params)
    V_sphere = (4.0 / 3.0) * math.pi * r_particle ** 3

    pf_static = float(mdl['packing_fraction_static'])
    bed_expansion_ratio = float(mdl['bed_expansion_ratio'])
    pf = pf_static if state == 'collapsed' else pf_static / bed_expansion_ratio

    r_throat = dim['nozzle']['throat_diameter_cm'] / 2.0
    z_cone_top = dim['cone']['vertical_drop_cm']
    r_retort = dim['retort']['id_cm'] / 2.0
    fluidized_height_cm = dim['bed']['fluidized_height_cm']

    if charge_mass_g is None:
        charge_mass_g = float(dim['bed']['charge_mass_g'])

    from furnace.geometry import _particle_effective_density
    rho_eff = _particle_effective_density(stage, params)
    V_solid = charge_mass_g / rho_eff
    V_bulk = V_solid / pf

    n_slabs = int(mdl['n_slabs'])
    slabs, _, _ = staircase_bed(params, n_slabs)

    cache_dir = _REPO_ROOT / "cases" / ".triso_cache"
    state_label = f"{state}_{charge_mass_g:.1f}g"

    tasks = []
    remaining_bulk = V_bulk
    for i, slab in enumerate(slabs):
        if remaining_bulk <= 1e-12:
            break

        r_slab = slab['radius']
        z_bot = slab['z_bot']
        z_top = slab['z_top']
        slab_vol = slab['volume']

        if remaining_bulk >= slab_vol:
            used_vol = slab_vol
            z_fill = z_top
        else:
            h_partial = remaining_bulk / (math.pi * r_slab ** 2)
            z_fill = z_bot + h_partial
            used_vol = remaining_bulk

        remaining_bulk -= used_vol
        seed_i = seed + i

        tasks.append((z_bot, z_fill, r_slab, pf, r_particle, seed_i,
                      str(cache_dir), i, state_label))

    # Overflow into cylinder
    if remaining_bulk > 1e-12:
        h_overflow = remaining_bulk / (math.pi * r_retort ** 2)
        if state == 'fluidized' and h_overflow > fluidized_height_cm:
            raise ValueError("Overflow exceeds fluidized_height_cm")
        z_ov_bot = z_cone_top
        z_ov_top = z_cone_top + h_overflow
        seed_ov = seed + len(slabs)
        tasks.append((z_ov_bot, z_ov_top, r_retort, pf, r_particle, seed_ov,
                      str(cache_dir), len(slabs), f"{state_label}_overflow"))

    return tasks


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('-j', '--jobs', type=int, default=None,
                    help='max parallel workers (default: logical CPUs - 1)')
    args = ap.parse_args()

    params = load_params()
    n_workers = args.jobs or max(1, os.cpu_count() - 1)
    print(f"Using {n_workers} workers")

    seed = int(params['model']['seed'])

    configs = [
        ('fluidized', None, seed, "fluidized 95g"),
        ('collapsed', None, seed, "collapsed 95g"),
    ]

    for state, mass, s, label in configs:
        tasks = build_slab_tasks(params, state, s, charge_mass_g=mass)
        n_total = len(tasks)
        print(f"\n[{label}] {n_total} slabs to process")
        t0 = time.time()
        n_cached = 0
        n_packed = 0
        n_error = 0

        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futs = {pool.submit(pack_slab_worker, t): t for t in tasks}
            for fut in as_completed(futs):
                lbl, n_p, elapsed, err = fut.result()
                if err:
                    print(f"  ERROR {lbl}: {err}")
                    n_error += 1
                elif elapsed == 0.0:
                    n_cached += 1
                else:
                    n_packed += 1
                    print(f"  packed {lbl}: {n_p} particles in {elapsed:.1f}s")

        total = time.time() - t0
        print(f"[{label}] done in {total:.1f}s — "
              f"{n_cached} cached hits, {n_packed} newly packed, {n_error} errors")

    cache_dir = _REPO_ROOT / "cases" / ".triso_cache"
    total_files = len(list(cache_dir.glob("*.npz")))
    print(f"\nTotal cache files: {total_files}")


if __name__ == "__main__":
    main()
