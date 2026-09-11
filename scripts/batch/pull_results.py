#!/usr/bin/env python
"""Pull sweep results from S3 into results/<sweep>/<submission>/ and aggregate.

Usage:
    python scripts/batch/pull_results.py --sweep step7a_mass                 # latest submission
    python scripts/batch/pull_results.py --sweep step7a_mass --submission ID # specific submission

Writes a CSV at results/<sweep>/<submission>/summary.csv joining every
results.json row with its k-eff, k-eff+2σ, and metadata for plotting.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import boto3


DEFAULT_STACK = "triso-ncs"
_REPO_ROOT = Path(__file__).resolve().parents[2]


def _stack_outputs(stack_name: str, region: str | None) -> dict:
    cfn = boto3.client("cloudformation", region_name=region)
    outs = cfn.describe_stacks(StackName=stack_name)["Stacks"][0].get("Outputs", [])
    return {o["OutputKey"]: o["OutputValue"] for o in outs}


def _latest_submission(bucket: str, sweep: str, region: str | None) -> str:
    s3 = boto3.client("s3", region_name=region)
    prefix = f"runs/{sweep}/"
    paginator = s3.get_paginator("list_objects_v2")
    submissions = set()
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix, Delimiter="/"):
        for cp in page.get("CommonPrefixes", []) or []:
            submissions.add(cp["Prefix"].split("/")[-2])
    if not submissions:
        raise SystemExit(f"no submissions found under s3://{bucket}/{prefix}")
    return sorted(submissions)[-1]


def _sync(bucket: str, sweep: str, submission: str, out_dir: Path) -> None:
    src = f"s3://{bucket}/runs/{sweep}/{submission}/"
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["aws", "s3", "sync", src, str(out_dir)], check=True)


def _aggregate(out_dir: Path) -> Path:
    rows = []
    for rj in sorted(out_dir.glob("*/results.json")):
        try:
            data = json.loads(rj.read_text())
        except json.JSONDecodeError:
            print(f"[pull] skipping malformed {rj}", file=sys.stderr)
            continue
        flat = {
            "array_index": rj.parent.name,
            "keff":            data.get("keff"),
            "keff_std":        data.get("keff_std"),
            "keff_plus_2sigma": data.get("keff_plus_2sigma"),
            "u235_mass_g":     data.get("u235_mass_g"),
            "pf_achieved":     data.get("pf_achieved"),
            "bed_height_cm":   data.get("bed_height_cm"),
            "wall_seconds":    data.get("wall_seconds"),
            "sweep_name":      data.get("sweep_name"),
            "git_sha":         data.get("git_sha"),
        }
        for k, v in (data.get("row") or {}).items():
            flat[f"row.{k}"] = v
        rows.append(flat)

    if not rows:
        raise SystemExit(f"no results.json files found in {out_dir}")

    fieldnames = sorted({k for r in rows for k in r.keys()})
    csv_path = out_dir / "summary.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return csv_path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--submission", default=None,
                    help="Submission ID; defaults to the newest one for this sweep.")
    ap.add_argument("--stack", default=DEFAULT_STACK)
    ap.add_argument("--region", default=None)
    args = ap.parse_args()

    outs = _stack_outputs(args.stack, args.region)
    bucket = outs["ResultsBucket"]

    submission = args.submission or _latest_submission(bucket, args.sweep, args.region)
    print(f"[pull] sweep={args.sweep}  submission={submission}")

    out_dir = _REPO_ROOT / "results" / args.sweep / submission
    _sync(bucket, args.sweep, submission, out_dir)
    csv_path = _aggregate(out_dir)
    print(f"[pull] wrote {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
