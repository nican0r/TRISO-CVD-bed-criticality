#!/usr/bin/env bash
# One-time: upload the local nuclear-data library to S3, then sync it onto EFS
# via a Fargate hydrator task. Idempotent — safe to rerun.
#
# Usage:
#   ./aws/hydrate-nuclear-data.sh [--stack triso-ncs] [--src ~/Documents/nuclear-data]
set -euo pipefail

STACK="triso-ncs"
SRC="$HOME/Documents/nuclear-data"
REGION="${AWS_REGION:-${AWS_DEFAULT_REGION:-us-east-1}}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --stack)  STACK="$2"; shift 2 ;;
        --src)    SRC="$2";   shift 2 ;;
        --region) REGION="$2"; shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

if [[ ! -d "$SRC/endfb80_hdf5" ]]; then
    echo "expected $SRC/endfb80_hdf5 to exist (ENDF/B-VIII.0 HDF5 library)" >&2
    exit 1
fi

BUCKET=$(aws cloudformation describe-stacks --stack-name "$STACK" --region "$REGION" \
    --query "Stacks[0].Outputs[?OutputKey=='ResultsBucket'].OutputValue" --output text)
HYDRATE_TASK=$(aws cloudformation describe-stacks --stack-name "$STACK" --region "$REGION" \
    --query "Stacks[0].Outputs[?OutputKey=='HydratorTaskDefinition'].OutputValue" --output text)
CLUSTER=$(aws cloudformation describe-stacks --stack-name "$STACK" --region "$REGION" \
    --query "Stacks[0].Outputs[?OutputKey=='HydratorCluster'].OutputValue" --output text)
SUBNET=$(aws cloudformation describe-stacks --stack-name "$STACK" --region "$REGION" \
    --query "Stacks[0].Outputs[?OutputKey=='HydratorSubnet'].OutputValue" --output text)
SG=$(aws cloudformation describe-stacks --stack-name "$STACK" --region "$REGION" \
    --query "Stacks[0].Outputs[?OutputKey=='HydratorSecurityGroup'].OutputValue" --output text)

echo "[hydrate] bucket=$BUCKET"
echo "[hydrate] uploading $SRC → s3://$BUCKET/nuclear-data/ (this can take 5-10 min)"

aws s3 sync "$SRC" "s3://$BUCKET/nuclear-data/" \
    --region "$REGION" \
    --exclude ".*" --exclude "*/.*"

echo "[hydrate] launching Fargate task to copy S3 → EFS"
TASK_ARN=$(aws ecs run-task \
    --region "$REGION" \
    --cluster "$CLUSTER" \
    --task-definition "$HYDRATE_TASK" \
    --launch-type FARGATE \
    --network-configuration "awsvpcConfiguration={subnets=[$SUBNET],securityGroups=[$SG],assignPublicIp=ENABLED}" \
    --query 'tasks[0].taskArn' --output text)

echo "[hydrate] task: $TASK_ARN"
echo "[hydrate] waiting for task to finish..."
aws ecs wait tasks-stopped --region "$REGION" --cluster "$CLUSTER" --tasks "$TASK_ARN"

EXIT_CODE=$(aws ecs describe-tasks --region "$REGION" --cluster "$CLUSTER" --tasks "$TASK_ARN" \
    --query 'tasks[0].containers[0].exitCode' --output text)

if [[ "$EXIT_CODE" != "0" ]]; then
    echo "[hydrate] task failed with exit code $EXIT_CODE" >&2
    echo "check CloudWatch Logs group /aws/ecs/triso-ncs-hydrator" >&2
    exit 1
fi

echo "[hydrate] done. EFS is populated with nuclear-data/."
