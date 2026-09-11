#!/usr/bin/env python
"""Submit a sweep to AWS Batch.

Usage:
    python scripts/batch/submit_sweep.py \
        --sweep step7a_mass \
        --manifest manifests/step7a_mass.json \
        [--dry-run]

The manifest is a JSON list of dicts, each dict a **kwargs row for
`furnace.model.build_model`. Example:

    [
      {"charge_mass_g": 23.75, "state": "fluidized"},
      {"charge_mass_g": 47.5,  "state": "fluidized"},
      {"charge_mass_g": 95.0,  "state": "fluidized"},
      ...
    ]

The script:
    1. Reads CloudFormation stack outputs (bucket name, queue arn, job def arn).
    2. Uploads the manifest to s3://<bucket>/manifests/<sweep>/<sha>.json.
    3. Calls batch:SubmitJob with arrayProperties.size = len(manifest),
       passing the manifest S3 URI and the results prefix as env overrides.
    4. Prints the submitted job ID and where results will land.

Requires: boto3, awscli configured, and the CloudFormation stack deployed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import boto3


DEFAULT_STACK = "triso-ncs"


def _stack_outputs(stack_name: str, region: str | None) -> dict:
    cfn = boto3.client("cloudformation", region_name=region)
    resp = cfn.describe_stacks(StackName=stack_name)
    outs = resp["Stacks"][0].get("Outputs", [])
    return {o["OutputKey"]: o["OutputValue"] for o in outs}


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=Path(__file__).resolve().parents[2],
            text=True,
        ).strip()
    except Exception:
        return "nogit"


def _upload_manifest(bucket: str, sweep: str, manifest_path: Path,
                     region: str | None) -> tuple[str, str]:
    body = manifest_path.read_bytes()
    digest = hashlib.sha256(body).hexdigest()[:12]
    key = f"manifests/{sweep}/{digest}.json"
    s3 = boto3.client("s3", region_name=region)
    s3.put_object(Bucket=bucket, Key=key, Body=body,
                  ContentType="application/json")
    return f"s3://{bucket}/{key}", digest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", required=True,
                    help="Sweep tag, e.g. step7a_mass")
    ap.add_argument("--manifest", required=True, type=Path,
                    help="Path to a JSON manifest (list of build_model kwargs).")
    ap.add_argument("--stack", default=DEFAULT_STACK,
                    help=f"CloudFormation stack name (default: {DEFAULT_STACK})")
    ap.add_argument("--region", default=None,
                    help="AWS region (defaults to boto3 profile)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would be submitted, do nothing.")
    args = ap.parse_args()

    manifest = json.loads(args.manifest.read_text())
    if not isinstance(manifest, list) or not manifest:
        print(f"[submit] manifest must be a non-empty JSON list", file=sys.stderr)
        return 2
    n = len(manifest)

    print(f"[submit] sweep={args.sweep}  n_jobs={n}")
    print(f"[submit] first row: {json.dumps(manifest[0])}")
    print(f"[submit] last  row: {json.dumps(manifest[-1])}")

    if args.dry_run:
        print("[submit] dry-run: not uploading, not submitting")
        return 0

    outs = _stack_outputs(args.stack, args.region)
    bucket = outs["ResultsBucket"]
    queue  = outs["JobQueue"]
    # Strip the revision suffix so Batch always picks the latest active revision.
    jobdef = outs["JobDefinition"].rsplit(":", 1)[0]

    manifest_uri, digest = _upload_manifest(bucket, args.sweep, args.manifest,
                                            args.region)
    submission_id = f"{int(time.time())}_{digest}"
    results_prefix = f"s3://{bucket}/runs/{args.sweep}/{submission_id}/"

    print(f"[submit] manifest → {manifest_uri}")
    print(f"[submit] results  → {results_prefix}")

    submit_kwargs = dict(
        jobName=f"{args.sweep}-{digest}",
        jobQueue=queue,
        jobDefinition=jobdef,
        containerOverrides={
            "environment": [
                {"name": "SWEEP_NAME",        "value": args.sweep},
                {"name": "PARAM_MANIFEST_S3", "value": manifest_uri},
                {"name": "RESULTS_PREFIX_S3", "value": results_prefix},
                {"name": "GIT_SHA",           "value": _git_sha()},
            ],
        },
        tags={
            "Project": "triso-ncs",
            "Sweep": args.sweep,
        },
    )
    if n > 1:
        submit_kwargs["arrayProperties"] = {"size": n}

    batch = boto3.client("batch", region_name=args.region)
    resp = batch.submit_job(**submit_kwargs)

    print(f"[submit] jobId = {resp['jobId']}")
    print(f"[submit] jobName = {resp['jobName']}")
    print()
    print("Monitor:")
    print(f"  aws batch describe-jobs --jobs {resp['jobId']} "
          f"--query 'jobs[0].{{status:status,attempts:attempts}}'")
    print("Pull results when SUCCEEDED:")
    print(f"  python scripts/batch/pull_results.py "
          f"--sweep {args.sweep} --submission {submission_id}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
