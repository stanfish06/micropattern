#!/usr/bin/env bash
set -euo pipefail

if [ -z "${MULTIVI_INPUT_S3:-}" ]; then
    echo "ERROR: Set MULTIVI_INPUT_S3 to the notebook-prepared integration-ready H5MU." >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

export ECR_REPO="micropattern-multivi"
export COMPUTE_ENV="micropattern-multivi-compute-env"
export JOB_QUEUE="micropattern-multivi-job-queue"
export JOB_DEF="micropattern-multivi-job-def"
export LOG_GROUP="/aws/batch/micropattern-multivi"
export JOB_NAME_PREFIX="micropattern-multivi"
export DOCKERFILE="${SCRIPT_DIR}/Dockerfile.multivi"

# Cost-first trial: g5.2xlarge = one 24 GB A10G GPU, 8 vCPU, 32 GiB host RAM.
# If the container is OOM-killed (RAM peak comes from the full-ATAC copy in
# prepare_hvg_tier), bump to g5.4xlarge (64 GiB) by restoring
# INSTANCE_TYPES_JSON=["g5.4xlarge"], MAX_VCPUS/JOB_VCPUS=16, JOB_MEMORY=60000.
export INSTANCE_TYPES_JSON='["g5.2xlarge"]'
export MAX_VCPUS="8"
export JOB_VCPUS="8"
export JOB_MEMORY="28000"
export JOB_GPUS="1"
export JOB_TIMEOUT_SECONDS="${JOB_TIMEOUT_SECONDS:-172800}"

MULTIVI_OUTPUT_S3="${MULTIVI_OUTPUT_S3:-}"
MULTIVI_HVG_VALUES="${MULTIVI_HVG_VALUES:-1500,2000}"
MULTIVI_ATAC_THRESHOLD="${MULTIVI_ATAC_THRESHOLD:-0.05}"
MULTIVI_IMPUTE_CHUNK_SIZE="${MULTIVI_IMPUTE_CHUNK_SIZE:-2000}"
MULTIVI_MAX_EPOCHS="${MULTIVI_MAX_EPOCHS:-500}"
MULTIVI_BATCH_SIZE="${MULTIVI_BATCH_SIZE:-256}"

ENVIRONMENT_JSON=$(jq -nc \
    --arg input "${MULTIVI_INPUT_S3}" \
    --arg output "${MULTIVI_OUTPUT_S3}" \
    --arg hvg_values "${MULTIVI_HVG_VALUES}" \
    --arg threshold "${MULTIVI_ATAC_THRESHOLD}" \
    --arg chunk_size "${MULTIVI_IMPUTE_CHUNK_SIZE}" \
    --arg max_epochs "${MULTIVI_MAX_EPOCHS}" \
    --arg batch_size "${MULTIVI_BATCH_SIZE}" \
    '[
        {name: "MULTIVI_INPUT_S3", value: $input},
        {name: "MULTIVI_HVG_VALUES", value: $hvg_values},
        {name: "MULTIVI_ATAC_THRESHOLD", value: $threshold},
        {name: "MULTIVI_IMPUTE_CHUNK_SIZE", value: $chunk_size},
        {name: "MULTIVI_MAX_EPOCHS", value: $max_epochs},
        {name: "MULTIVI_BATCH_SIZE", value: $batch_size}
    ] + (if $output == "" then [] else [{name: "MULTIVI_OUTPUT_S3", value: $output}] end)'
)
export JOB_CONTAINER_OVERRIDES_JSON
JOB_CONTAINER_OVERRIDES_JSON=$(jq -nc --argjson environment "${ENVIRONMENT_JSON}" \
    '{environment: $environment}')

exec "${SCRIPT_DIR}/deploy.sh" "${1:-ondemand}"
