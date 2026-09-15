"""Cone-slab discretisation convergence study — packed-particle approach.

Replaces the former homogenised-surrogate study with a real packed-particle
design that produces physically meaningful bias estimates.

Usage
-----
Generate manifest and print submit command:
    python scripts/run_convergence.py

Generate manifest and immediately submit to AWS Batch:
    python scripts/run_convergence.py --submit

Dry-run geometry check (build one model locally, no OpenMC run):
    python scripts/run_convergence.py --local-check

Study design
------------
EXACT-CONE REFERENCE
    A reference geometry that packs the true frustum by rejection sampling:
    pack the bounding cylinder (r=r_retort, z in [0, z_cone_top]), then discard
    every kernel whose sphere intersects the cone surface (perpendicular
    clearance < r_kernel) or the floor/ceiling planes.  Iterate pf_trial until
    the surviving count matches the charge mass to within 0.5 %.

SEED REPLICATES
    Five independent packing seeds per configuration.  The across-seed standard
    deviation is the error bar for all comparisons, not the per-run MC sigma.

RUN MATRIX
    Fixed: collapsed bed state, bare_kernel stage, nominal charge mass,
    full-flood background (z_flood=999 cm, water at 1.0 g/cm³).
    Full flood is the bounding case for criticality; the convergence study
    measures bias at worst-case moderation rather than at nominal gas conditions.
    Configurations: exact_cone (reference), n=32 (current production), n=16, n=8.
    5 seeds × 4 configurations = 20 jobs.

OUTPUT (after pull_results.py)
    Table: config | mean k-eff | across-seed σ | mean MC σ |
    achieved pf | requested pf | pf deviation | bias vs exact-cone.

    The achieved vs requested packing fraction per configuration is the primary
    diagnostic: interface depletion (fewer particles placed near slab boundaries
    than the bulk packing fraction implies) shows up as achieved pf < requested pf,
    and the deficit grows with decreasing n_slabs.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from furnace.params import check_env, load_params

# ── Constants ────────────────────────────────────────────────────────────────

# Staircase slab counts included in the study.
# n=32 is the current production value; n=16 and n=8 are coarser alternatives.
_STAIRCASE_NS = [32, 16, 8]

# Seeds for the five replicates (consecutive integers from the default seed).
_SEEDS = [42, 43, 44, 45, 46]

# MC settings for each job (match production).
_N_PARTICLES = 20_000
_N_INACTIVE  = 50
_N_ACTIVE    = 250

# z_flood value that submerges the entire retort (all geometry < 72 cm).
# Full flood is the bounding case for criticality; the convergence study runs
# flood-only rather than gas + flood so that bias estimates are made at the
# worst-case moderation condition.
_FLOOD_Z = 999.0

SWEEP_NAME = 'step3_convergence'


# ── Run-matrix builder ────────────────────────────────────────────────────────

def _run_matrix() -> list[dict]:
    """Return a flat list of manifest rows: 4 configs × 5 seeds = 20 jobs (flood only)."""
    n_jobs = (1 + len(_STAIRCASE_NS)) * len(_SEEDS)
    print(f'Run matrix: exact_cone + staircase {_STAIRCASE_NS}, '
          f'{len(_SEEDS)} seeds, full-flood background = {n_jobs} jobs.')

    rows = []
    base = dict(
        state='collapsed',
        stage='bare_kernel',
        n_particles=_N_PARTICLES,
        n_inactive=_N_INACTIVE,
        n_active=_N_ACTIVE,
    )

    for seed in _SEEDS:
        row = dict(**base, seed=seed, use_exact_cone=True, z_flood=_FLOOD_Z)
        row['_config'] = 'exact_cone'
        rows.append(row)

        for n in _STAIRCASE_NS:
            row = dict(**base, seed=seed, n_slabs=n, z_flood=_FLOOD_Z)
            row['_config'] = f'n{n}'
            rows.append(row)

    return rows


# ── Manifest writer ───────────────────────────────────────────────────────────

def _strip_meta(rows: list[dict]) -> list[dict]:
    return [{k: v for k, v in row.items() if not k.startswith('_')} for row in rows]


def write_manifest(rows: list[dict]) -> Path:
    manifest_path = _REPO_ROOT / 'manifests' / f'{SWEEP_NAME}.json'
    manifest_path.parent.mkdir(exist_ok=True)
    manifest_path.write_text(json.dumps(_strip_meta(rows), indent=2))
    print(f'Manifest written: {manifest_path}  ({len(rows)} rows)')
    return manifest_path


# ── Local geometry check ──────────────────────────────────────────────────────

def local_check(params, rows: list[dict]) -> None:
    """Build and export one exact_cone and one staircase model without running OpenMC.

    Uses charge_mass_g=5.0 (matching the smoke-test convention) so that
    pack_spheres finishes in seconds rather than the ~15 min required for the
    full 95 g / ~225k-kernel production pack.  Geometry, cell topology, and
    material wiring are identical to the full run; only particle count differs.
    """
    from furnace.model import build_model

    # Pick one exact_cone row and one staircase row (n=32, most particles per slab)
    check_rows = []
    seen_exact = seen_staircase = False
    for r in _strip_meta(rows):
        if not seen_exact and r.get('use_exact_cone'):
            check_rows.append(('exact_cone', r))
            seen_exact = True
        elif not seen_staircase and 'n_slabs' in r:
            check_rows.append((f"n{r['n_slabs']}", r))
            seen_staircase = True
        if seen_exact and seen_staircase:
            break

    for label, row in check_rows:
        row = dict(row, charge_mass_g=5.0)   # small charge — geometry check only
        print(f'\nLocal check [{label}]: {json.dumps(row)}')
        model, stats = build_model(params, **row)
        with tempfile.TemporaryDirectory(prefix=f'conv_check_{label}_') as tmpdir:
            model.geometry.export_to_xml(str(Path(tmpdir) / 'geometry.xml'))
            model.materials.export_to_xml(str(Path(tmpdir) / 'materials.xml'))
            model.settings.export_to_xml(str(Path(tmpdir) / 'settings.xml'))
            print(f'  OK — n_particles={stats.n_particles}, '
                  f'pf_achieved={stats.pf_achieved:.4f}, bed_top={stats.bed_height_cm:.3f} cm')
    print('\nLocal geometry check passed.')


# ── Summary table ─────────────────────────────────────────────────────────────

def print_matrix_summary(rows: list[dict]) -> None:
    print(f'\n{"#":>4}  {"config":<12}  {"seed":>6}  {"n_slabs":>8}  {"z_flood":>8}')
    print('-' * 46)
    for i, row in enumerate(rows):
        cfg = row.get('_config', '?')
        ns  = row.get('n_slabs', '—')
        zf  = row.get('z_flood', '—')
        s   = row.get('seed', '?')
        print(f'{i:>4}  {cfg:<12}  {s:>6}  {str(ns):>8}  {str(zf):>8}')


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--local-check', action='store_true',
                    help='Build one model locally (no OpenMC run) to verify geometry.')
    ap.add_argument('--submit', action='store_true',
                    help='Submit the manifest to AWS Batch after writing it.')
    ap.add_argument('--dry-run', action='store_true',
                    help='Pass --dry-run to submit_sweep.py (print only, no Batch call).')
    ap.add_argument('--region', default=None,
                    help='AWS region for Batch submission (e.g. sa-east-1).')
    args = ap.parse_args()

    check_env()
    params = load_params()

    rows = _run_matrix()
    print_matrix_summary(rows)
    manifest_path = write_manifest(rows)

    if args.local_check:
        local_check(params, rows)

    submit_cmd = (
        f'python scripts/batch/submit_sweep.py '
        f'--sweep {SWEEP_NAME} '
        f'--manifest {manifest_path}'
    )
    if args.region:
        submit_cmd += f' --region {args.region}'

    print(f'\nSubmit command:\n  {submit_cmd}')
    print('\nPull results after completion:')
    print(f'  python scripts/batch/pull_results.py --sweep {SWEEP_NAME} --submission <id>')

    if args.submit or args.dry_run:
        cmd = submit_cmd.split()
        if args.dry_run:
            cmd.append('--dry-run')
        result = subprocess.run(cmd, cwd=str(_REPO_ROOT))
        return result.returncode

    return 0


if __name__ == '__main__':
    sys.exit(main())
