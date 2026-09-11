# Running sweeps on AWS Batch

Steps 7 (mass / packing sweeps) and 8 (bottom-up flood accident cases)
launch hundreds of independent OpenMC eigenvalue runs. This document explains
how to run them on AWS Batch instead of the laptop, and what it will cost.

The physics code (`furnace/*`) is unchanged. Only the driver moves — a small
container calls `furnace.model.build_model(**row)` for one row of a JSON
manifest, then uploads the statepoint and a summary JSON to S3.

---

## Contents

- [Architecture at a glance](#architecture-at-a-glance)
- [One-time setup](#one-time-setup)
- [Running a sweep](#running-a-sweep)
- [Manifest format](#manifest-format)
- [Monitoring](#monitoring)
- [Cost](#cost)
- [Teardown](#teardown)
- [Troubleshooting](#troubleshooting)

---

## Architecture at a glance

```
local (Mac)                                 AWS
────────────                                ───
scripts/batch/submit_sweep.py  ── upload ─▶ S3://…/manifests/<sweep>/<sha>.json
                               ── submit ─▶ AWS Batch (EC2 Spot, c7i.2xlarge)
                                              │
                                              ├─ pulls container from ECR
                                              ├─ mounts EFS /xs (nuclear data)
                                              └─ runs entrypoint.py
                                                    └─ writes → S3://…/runs/<sweep>/<sub>/<idx>/
scripts/batch/pull_results.py  ◀── sync ── S3://…/runs/<sweep>/<sub>/
```

All AWS resources are created by **`aws/cloudformation.yaml`**. Everything
lives in a single stack so you can tear it down with one command when the
study is finished.

- **S3 bucket** `triso-cvd-ncs-<account>-<region>` — manifests, results,
  and the source copy of nuclear data. Lifecycle: results transition to
  Infrequent Access after 30 days and expire after 180 days.
- **EFS filesystem** mounted read-only at `/xs` in every job. Hydrated once
  from S3 by a Fargate task.
- **ECR repo** `triso-cvd-ncs` — the job container image.
- **Batch compute environment** on EC2 **Spot** (`c7i.2xlarge` /
  `c6i.2xlarge`, 8 vCPU / 16 GB), min 0 → max 256 vCPU (32 concurrent
  8-vCPU jobs by default).
- **Job queue** `triso-ncs-queue` and **job definition** `triso-ncs-job`.

---

## One-time setup

### Prerequisites

- `aws` CLI configured (`aws configure`) with a profile that has admin
  or equivalent stack-creation permissions in the target account.
- `docker` installed and running.
- `boto3` installed locally (`pip install boto3`).
- Local nuclear-data library at `~/Documents/nuclear-data/endfb80_hdf5/`
  (already set on this machine — see project `CLAUDE.md`).

### 1. Pick a VPC and subnets

The stack needs an existing VPC and at least one **public** subnet.
Using public subnets avoids the ~$30/month standing charge of a NAT
gateway; Spot instances get public IPs and reach ECR/S3 directly.

Discover the default VPC and its public subnets:

```bash
aws ec2 describe-vpcs \
    --filters Name=is-default,Values=true \
    --query 'Vpcs[0].VpcId' --output text
# vpc-0123456789abcdef0

aws ec2 describe-subnets \
    --filters Name=vpc-id,Values=<vpc-id> Name=map-public-ip-on-launch,Values=true \
    --query 'Subnets[].SubnetId' --output text
# subnet-aaa subnet-bbb subnet-ccc
```

### 2. Deploy the stack

```bash
aws cloudformation deploy \
    --template-file aws/cloudformation.yaml \
    --stack-name triso-ncs \
    --parameter-overrides \
        VpcId=vpc-0123456789abcdef0 \
        SubnetIds=subnet-aaa,subnet-bbb,subnet-ccc \
    --capabilities CAPABILITY_IAM
```

First deploy takes ~4–6 minutes.

Check outputs (bucket name, ECR URI, queue ARN, etc.):

```bash
aws cloudformation describe-stacks --stack-name triso-ncs \
    --query 'Stacks[0].Outputs'
```

### 3. Upload nuclear data (5–10 min, one-time)

```bash
./aws/hydrate-nuclear-data.sh
```

This uploads `~/Documents/nuclear-data/` (~5 GB) to
`s3://triso-cvd-ncs-…/nuclear-data/`, then launches a small Fargate
task that syncs S3 → EFS at `/xs/`.

### 4. Build and push the container

```bash
./aws/build-and-push.sh
```

Builds `aws/Dockerfile` as `linux/amd64` (so it runs on x86 Batch hosts
even from an Apple Silicon Mac) and pushes to ECR. The image is
~1 GB; subsequent pushes only send changed layers.

### 5. Sanity test

Two cheap checks before spending on the real sweep.

**a. Container-only smoke test.** Runs a 200-particle nominal case
locally against your bind-mounted nuclear-data dir. Zero AWS cost.

```bash
docker run --rm --platform=linux/amd64 \
    -v "$HOME/Documents/nuclear-data:/xs:ro" \
    -e SMOKE=1 \
    -e AWS_BATCH_JOB_ARRAY_INDEX=0 \
    -e SWEEP_NAME=smoke \
    triso-cvd-ncs:latest
```

You should see `k-eff = …` printed at the end.

**b. Single Batch job.** Submit a one-row manifest to prove IAM, EFS,
S3, and CloudWatch plumbing all work end-to-end.

```bash
cat > /tmp/smoke.json <<'JSON'
[{"state": "fluidized", "n_particles": 200, "n_inactive": 5, "n_active": 15, "seed": 1}]
JSON

python scripts/batch/submit_sweep.py --sweep smoke --manifest /tmp/smoke.json
```

When it finishes (`aws batch describe-jobs --jobs <id>` shows
`SUCCEEDED`), pull results:

```bash
python scripts/batch/pull_results.py --sweep smoke
cat results/smoke/*/summary.csv
```

Cost: well under $0.01.

---

## Running a sweep

Sweeps are driven by JSON manifests — one file per sweep, one row per
Batch array-child job. The submitter does not care what the sweep is
about; it just fans out the manifest.

```bash
python scripts/batch/submit_sweep.py \
    --sweep step7a_mass \
    --manifest manifests/step7a_mass.json
```

Output:

```
[submit] sweep=step7a_mass  n_jobs=10
[submit] manifest → s3://triso-cvd-ncs-…/manifests/step7a_mass/9f3c….json
[submit] results  → s3://triso-cvd-ncs-…/runs/step7a_mass/1734…_9f3c/
[submit] jobId = 6c8f…
```

Pull results (blocks / can be re-run):

```bash
python scripts/batch/pull_results.py --sweep step7a_mass
# → writes results/step7a_mass/<submission>/summary.csv
```

`summary.csv` has one row per array child: `keff`, `keff_std`,
`keff_plus_2sigma`, `u235_mass_g`, plus every key from the manifest row
prefixed with `row.`. Feed it directly into the plotting scripts for
steps 7/8.

Use `--dry-run` on `submit_sweep.py` to see the row count and first/last
rows without spending anything.

---

## Manifest format

Each row is passed straight into
`furnace.model.build_model(params, **row)` inside the container. Valid
keys today (from `furnace/model.py`):

| Key             | Type      | Notes                                       |
|-----------------|-----------|---------------------------------------------|
| `state`         | str       | `"fluidized"` or `"collapsed"`              |
| `stage`         | str       | TRISO deposition stage, e.g. `"bare_kernel"`|
| `background`    | str       | `"gas"` or `"water"`                        |
| `charge_mass_g` | float     | Overrides `dimensions.bed.charge_mass_g`    |
| `n_particles`   | int       | Particles per batch                         |
| `n_inactive`    | int       | Inactive batches                            |
| `n_active`      | int       | Active batches                              |
| `seed`          | int       | OpenMC + packing seed                       |

Steps 7b / 8 (packing fraction, water density, partial flood level)
will need `build_model` to accept new kwargs. Add them to the signature
in `furnace/model.py`; the container picks them up automatically as
long as the manifest key matches.

Example step 7a manifest:

```json
[
  {"state": "fluidized", "charge_mass_g": 23.75},
  {"state": "fluidized", "charge_mass_g": 47.5},
  {"state": "fluidized", "charge_mass_g": 95.0},
  {"state": "fluidized", "charge_mass_g": 190.0},
  {"state": "fluidized", "charge_mass_g": 475.0},
  {"state": "fluidized", "charge_mass_g": 950.0}
]
```

---

## Monitoring

```bash
# All jobs for a sweep:
aws batch list-jobs --job-queue triso-ncs-queue \
    --filters name=JOB_NAME,values=step7a_mass-*

# Detail on one submission:
aws batch describe-jobs --jobs <jobId> \
    --query 'jobs[0].{status:status,attempts:length(attempts)}'

# Live logs from array child N:
aws logs tail /aws/batch/triso-ncs --follow \
    --log-stream-name-prefix job/<jobId>:<N>
```

The console (`Batch → Jobs`) is easier for eyeballing an array in flight.

---

## Cost

Estimates for the full step 7 + 8 study. **All figures are in USD**,
current us-east-1 prices as of 2026 Q3.

### Per-job runtime and price

- Nominal case at production settings (20 k particles × 300 batches,
  6 M histories) is ~10–15 min on a c7i.2xlarge. Denser bed states and
  flooded configs raise collision rate; assume **20 min/job average**.
- c7i.2xlarge Spot: ~$0.06/hr → **$0.02 per job**
- c7i.2xlarge On-Demand: ~$0.11/hr → **$0.037 per job** (ceiling)

### Job counts

| Sweep                           | Jobs |
|---------------------------------|------|
| 7a — mass / batch size          | 10   |
| 7b — packing at fixed mass      | 12   |
| 7b — packing at fixed volume    | 12   |
| 8  — bottom-up flood (z × ρ × state): 2×20×2        | 80 |
| 8  — seed sensitivity, top-10 reactive     | 30 |
| **Total**                       | **~144** |

### Compute cost

- **Spot:** 146 × $0.02 ≈ **$3**
- **On-Demand ceiling:** 146 × $0.037 ≈ **$5.50**

### Fixed / recurring

| Item              | Notes                       | Cost |
|-------------------|-----------------------------|------|
| EFS storage       | 5 GB × $0.30/GB-mo          | $1.50/mo |
| S3 storage        | ~11 GB of statepoints       | $0.25/mo |
| ECR storage       | ~1 GB image                 | $0.10/mo |
| Data egress       | Pulling 11 GB back to Mac   | ~$1 one-time |
| CloudWatch Logs   | Negligible                  | <$0.50 |

### Bottom line

**Realistic total for the complete step 7 + 8 study: $5–10.**
Worst case, if you doubled batch counts across all sweeps, you would
still be under **$25**.

Cost controls in place:
- Compute environment `MinvCpus=0`: nothing runs when no job is queued.
- Job definition `Attempts=2`: retries once on Spot reclamation, then
  gives up rather than looping.
- Job timeout 3 hours per attempt: a runaway won't burn a whole day.
- `submit_sweep.py --dry-run`: preview a sweep before spending.
- Stack tag `Project=triso-ncs` on every resource: enable this tag in
  Cost Explorer to see spend broken out.

---

## Teardown

```bash
aws cloudformation delete-stack --stack-name triso-ncs
```

The Lambda-backed bucket emptier drains S3 first so the delete
succeeds. The ECR repo is set to `EmptyOnDelete=true`, so it goes with
the stack. Total teardown: ~3 minutes.

---

## Troubleshooting

**Jobs stuck in `RUNNABLE`.** Almost always an IAM/instance-profile
mismatch, an ENI limit in the subnet, or the compute environment has
been disabled. Check:
```bash
aws batch describe-compute-environments --compute-environments triso-ncs-spot
aws batch describe-job-queues --job-queues triso-ncs-queue
```

**Job fails immediately with `CannotPullContainerError`.** ECR image
tag is wrong or was never pushed. Rerun `./aws/build-and-push.sh`.

**Job fails with `openmc: error: cross_sections.xml not found`.** EFS
hydration never ran or failed. Rerun `./aws/hydrate-nuclear-data.sh`
and check the CloudWatch log group `/aws/ecs/triso-ncs-hydrator`.

**Spot reclaim shows up as `FAILED` with a `Host EC2*` status reason.**
Batch will already have retried once (per the job definition's retry
strategy). If retries also fail, wait for Spot capacity to improve or
temporarily set `SpotBidPercentage=100` and swap to On-Demand by
switching the CE `Type` to `EC2`.

**`aws ecr get-login-password` fails.** Your local AWS profile expired
(SSO). Refresh with `aws sso login` and retry.
