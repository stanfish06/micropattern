#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

export ECR_REPO="micropattern-multivi"
export COMPUTE_ENV="micropattern-multivi-compute-env"
export JOB_QUEUE="micropattern-multivi-job-queue"
export JOB_DEF="micropattern-multivi-job-def"
export LOG_GROUP="/aws/batch/micropattern-multivi"

# The MultiVI and original integration environments share these roles.
export DELETE_SHARED_IAM_ROLES="0"

exec "${SCRIPT_DIR}/cleanup.sh"
