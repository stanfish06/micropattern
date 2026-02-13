#!/usr/bin/env bash
set -euo pipefail

REGION="us-east-1"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
ECR_REPO="micropattern"

COMPUTE_ENV="micropattern-compute-env"
JOB_QUEUE="micropattern-job-queue"
JOB_DEF="micropattern-job-def"
LOG_GROUP="/aws/batch/micropattern"

BATCH_SERVICE_ROLE="micropattern-batch-service-role"
ECS_EXECUTION_ROLE="micropattern-ecs-execution-role"
ECS_TASK_ROLE="micropattern-ecs-task-role"

echo "=== Step 1: Cancel Running Jobs ==="

for status in SUBMITTED PENDING RUNNABLE STARTING RUNNING; do
    JOBS=$(aws batch list-jobs --job-queue "${JOB_QUEUE}" --job-status "${status}" \
        --query "jobSummaryList[*].jobId" --output text 2>/dev/null || true)
    for job_id in ${JOBS}; do
        echo "Cancelling job ${job_id}..."
        aws batch cancel-job --job-id "${job_id}" --reason "Cleanup" 2>/dev/null || \
        aws batch terminate-job --job-id "${job_id}" --reason "Cleanup" 2>/dev/null || true
    done
done

echo "=== Step 2: Disable and Delete Job Queue ==="

if aws batch describe-job-queues --job-queues "${JOB_QUEUE}" \
    --query "jobQueues[?status!='DELETED']|[0].jobQueueName" --output text 2>/dev/null | grep -q "${JOB_QUEUE}"; then
    echo "Disabling job queue..."
    aws batch update-job-queue --job-queue "${JOB_QUEUE}" --state DISABLED
    echo "Waiting for job queue to stabilize..."
    for i in $(seq 1 20); do
        STATUS=$(aws batch describe-job-queues --job-queues "${JOB_QUEUE}" \
            --query "jobQueues[0].status" --output text 2>/dev/null || echo "DELETED")
        [ "${STATUS}" = "VALID" ] && break
        sleep 15
    done
    echo "Deleting job queue..."
    aws batch delete-job-queue --job-queue "${JOB_QUEUE}"
    echo "Waiting for job queue deletion..."
    for i in $(seq 1 20); do
        if ! aws batch describe-job-queues --job-queues "${JOB_QUEUE}" \
            --query "jobQueues[?status!='DELETED']|[0]" --output text 2>/dev/null | grep -q "${JOB_QUEUE}"; then
            break
        fi
        sleep 15
    done
else
    echo "Job queue ${JOB_QUEUE} not found."
fi

echo "=== Step 3: Disable and Delete Compute Environment ==="

if aws batch describe-compute-environments --compute-environments "${COMPUTE_ENV}" \
    --query "computeEnvironments[?status!='DELETED']|[0].computeEnvironmentName" --output text 2>/dev/null | grep -q "${COMPUTE_ENV}"; then
    echo "Disabling compute environment..."
    aws batch update-compute-environment --compute-environment "${COMPUTE_ENV}" --state DISABLED
    echo "Waiting for compute environment to stabilize..."
    for i in $(seq 1 20); do
        STATUS=$(aws batch describe-compute-environments --compute-environments "${COMPUTE_ENV}" \
            --query "computeEnvironments[0].status" --output text 2>/dev/null || echo "DELETED")
        [ "${STATUS}" = "VALID" ] && break
        sleep 15
    done
    echo "Deleting compute environment..."
    aws batch delete-compute-environment --compute-environment "${COMPUTE_ENV}"
else
    echo "Compute environment ${COMPUTE_ENV} not found."
fi

echo "=== Step 4: Deregister Job Definitions ==="

REVISIONS=$(aws batch describe-job-definitions --job-definition-name "${JOB_DEF}" --status ACTIVE \
    --query "jobDefinitions[*].jobDefinitionArn" --output text 2>/dev/null || true)
for rev in ${REVISIONS}; do
    echo "Deregistering ${rev}..."
    aws batch deregister-job-definition --job-definition "${rev}"
done

echo "=== Step 5: Delete IAM Roles ==="

delete_role() {
    local role_name="$1"
    if ! aws iam get-role --role-name "${role_name}" &>/dev/null; then
        echo "Role ${role_name} not found, skipping."
        return
    fi

    POLICIES=$(aws iam list-attached-role-policies --role-name "${role_name}" \
        --query "AttachedPolicies[*].PolicyArn" --output text 2>/dev/null || true)
    for policy in ${POLICIES}; do
        echo "  Detaching ${policy}..."
        aws iam detach-role-policy --role-name "${role_name}" --policy-arn "${policy}"
    done

    INLINE=$(aws iam list-role-policies --role-name "${role_name}" \
        --query "PolicyNames[*]" --output text 2>/dev/null || true)
    for policy in ${INLINE}; do
        echo "  Deleting inline policy ${policy}..."
        aws iam delete-role-policy --role-name "${role_name}" --policy-name "${policy}"
    done

    echo "Deleting role ${role_name}..."
    aws iam delete-role --role-name "${role_name}"
}

delete_role "${BATCH_SERVICE_ROLE}"
delete_role "${ECS_EXECUTION_ROLE}"
delete_role "${ECS_TASK_ROLE}"
echo "Skipping ecsInstanceRole (may be shared with other services)."

echo "=== Step 6: Delete ECR Repository ==="

if aws ecr describe-repositories --repository-names "${ECR_REPO}" --region "${REGION}" &>/dev/null; then
    echo "Deleting ECR repository ${ECR_REPO}..."
    aws ecr delete-repository --repository-name "${ECR_REPO}" --region "${REGION}" --force
else
    echo "ECR repository ${ECR_REPO} not found."
fi

echo "=== Step 7: Delete CloudWatch Log Group ==="

if aws logs describe-log-groups --log-group-name-prefix "${LOG_GROUP}" \
    --query "logGroups[?logGroupName=='${LOG_GROUP}']" --output text 2>/dev/null | grep -q "${LOG_GROUP}"; then
    echo "Deleting log group ${LOG_GROUP}..."
    aws logs delete-log-group --log-group-name "${LOG_GROUP}"
else
    echo "Log group ${LOG_GROUP} not found."
fi

echo ""
echo "Cleanup complete."
