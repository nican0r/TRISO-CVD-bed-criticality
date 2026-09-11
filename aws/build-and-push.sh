#!/usr/bin/env bash
# Build the job container and push it to the ECR repo created by the CFN stack.
#
# Usage:
#   ./aws/build-and-push.sh                    # defaults: stack=triso-ncs, tag=latest
#   ./aws/build-and-push.sh --stack triso-ncs --tag v1
#
# Requires: docker, aws CLI configured, and the CFN stack already deployed.
set -euo pipefail

STACK="triso-ncs"
TAG="latest"
REGION="${AWS_REGION:-${AWS_DEFAULT_REGION:-us-east-1}}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --stack)  STACK="$2"; shift 2 ;;
        --tag)    TAG="$2";   shift 2 ;;
        --region) REGION="$2"; shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

REPO_URI=$(aws cloudformation describe-stacks \
    --stack-name "$STACK" --region "$REGION" \
    --query "Stacks[0].Outputs[?OutputKey=='EcrRepositoryUri'].OutputValue" \
    --output text)

if [[ -z "$REPO_URI" || "$REPO_URI" == "None" ]]; then
    echo "Could not read EcrRepositoryUri from stack $STACK in $REGION." >&2
    echo "Is the CloudFormation stack deployed?" >&2
    exit 1
fi

REGISTRY="${REPO_URI%%/*}"

echo "[build] repo:   $REPO_URI"
echo "[build] region: $REGION"
echo "[build] tag:    $TAG"

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

aws ecr get-login-password --region "$REGION" \
    | docker login --username AWS --password-stdin "$REGISTRY"

# --platform=linux/amd64 forces an x86_64 image so it runs on the Batch EC2
# fleet even when this script is executed from an Apple Silicon Mac.
docker build --platform=linux/amd64 \
    -f aws/Dockerfile \
    -t "$REPO_URI:$TAG" \
    .

docker push "$REPO_URI:$TAG"

echo "[build] pushed $REPO_URI:$TAG"
