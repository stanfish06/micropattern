#!/usr/bin/env bash
set -euo pipefail

USE_SPOT="${1:-ondemand}"  # pass "spot" as first arg to use Spot instances

export AWS_DEFAULT_REGION="us-east-2"
REGION="${AWS_DEFAULT_REGION}"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
ECR_REPO="micropattern"
IMAGE_URI="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${ECR_REPO}:latest"
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

S3_BUCKET="stan-sequencing-data"
COMPUTE_ENV="micropattern-compute-env"
JOB_QUEUE="micropattern-job-queue"
JOB_DEF="micropattern-job-def"
LOG_GROUP="/aws/batch/micropattern"

BATCH_SERVICE_ROLE="micropattern-batch-service-role"
ECS_EXECUTION_ROLE="micropattern-ecs-execution-role"
ECS_TASK_ROLE="micropattern-ecs-task-role"
EC2_INSTANCE_ROLE="ecsInstanceRole"

echo "=== Phase A: Container Image ==="

if ! aws ecr describe-repositories --repository-names "${ECR_REPO}" --region "${REGION}" &>/dev/null; then
    echo "Creating ECR repository ${ECR_REPO}..."
    aws ecr create-repository --repository-name "${ECR_REPO}" --region "${REGION}"
else
    echo "ECR repository ${ECR_REPO} already exists."
fi

echo "Logging in to ECR..."
aws ecr get-login-password --region "${REGION}" | \
    sudo docker login --username AWS --password-stdin "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"

echo "Building Docker image..."
sudo docker build -f "${PROJECT_ROOT}/scripts/aws/Dockerfile" -t "${ECR_REPO}" "${PROJECT_ROOT}"

echo "Tagging and pushing image..."
sudo docker tag "${ECR_REPO}:latest" "${IMAGE_URI}"
sudo docker push "${IMAGE_URI}"

echo "=== Phase B: IAM Roles ==="

create_role_if_not_exists() {
    local role_name="$1"
    local trust_policy="$2"
    if ! aws iam get-role --role-name "${role_name}" &>/dev/null; then
        echo "Creating IAM role ${role_name}..."
        aws iam create-role \
            --role-name "${role_name}" \
            --assume-role-policy-document "${trust_policy}"
    else
        echo "IAM role ${role_name} already exists."
    fi
}

BATCH_TRUST='{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"batch.amazonaws.com"},"Action":"sts:AssumeRole"}]}'
ECS_TASKS_TRUST='{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"ecs-tasks.amazonaws.com"},"Action":"sts:AssumeRole"}]}'
EC2_TRUST='{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}'

create_role_if_not_exists "${BATCH_SERVICE_ROLE}" "${BATCH_TRUST}"
aws iam attach-role-policy \
    --role-name "${BATCH_SERVICE_ROLE}" \
    --policy-arn "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole" 2>/dev/null || true

create_role_if_not_exists "${ECS_EXECUTION_ROLE}" "${ECS_TASKS_TRUST}"
aws iam attach-role-policy \
    --role-name "${ECS_EXECUTION_ROLE}" \
    --policy-arn "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy" 2>/dev/null || true

create_role_if_not_exists "${ECS_TASK_ROLE}" "${ECS_TASKS_TRUST}"

S3_POLICY=$(cat <<'POLICY'
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "s3:GetObject",
                "s3:PutObject",
                "s3:ListBucket",
                "s3:DeleteObject"
            ],
            "Resource": [
                "arn:aws:s3:::stan-sequencing-data",
                "arn:aws:s3:::stan-sequencing-data/*"
            ]
        }
    ]
}
POLICY
)
aws iam put-role-policy \
    --role-name "${ECS_TASK_ROLE}" \
    --policy-name "s3-access" \
    --policy-document "${S3_POLICY}"

create_role_if_not_exists "${EC2_INSTANCE_ROLE}" "${EC2_TRUST}"
aws iam attach-role-policy \
    --role-name "${EC2_INSTANCE_ROLE}" \
    --policy-arn "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role" 2>/dev/null || true

if ! aws iam get-instance-profile --instance-profile-name "${EC2_INSTANCE_ROLE}" &>/dev/null; then
    echo "Creating instance profile ${EC2_INSTANCE_ROLE}..."
    aws iam create-instance-profile --instance-profile-name "${EC2_INSTANCE_ROLE}"
    aws iam add-role-to-instance-profile \
        --instance-profile-name "${EC2_INSTANCE_ROLE}" \
        --role-name "${EC2_INSTANCE_ROLE}"
    echo "Waiting for instance profile propagation..."
    sleep 10
else
    echo "Instance profile ${EC2_INSTANCE_ROLE} already exists."
fi

echo "=== Phase C: Networking ==="

VPC_ID=$(aws ec2 describe-vpcs --filters "Name=isDefault,Values=true" --query "Vpcs[0].VpcId" --output text)
echo "Default VPC: ${VPC_ID}"

SUBNETS=$(aws ec2 describe-subnets \
    --filters "Name=vpc-id,Values=${VPC_ID}" \
    --query "Subnets[*].SubnetId" --output json)
echo "Subnets: ${SUBNETS}"

SG_ID=$(aws ec2 describe-security-groups \
    --filters "Name=vpc-id,Values=${VPC_ID}" "Name=group-name,Values=default" \
    --query "SecurityGroups[0].GroupId" --output text)
echo "Security Group: ${SG_ID}"

echo "=== Phase D: Compute Environment ==="

CE_STATUS=$(aws batch describe-compute-environments --compute-environments "${COMPUTE_ENV}" \
    --query "computeEnvironments[?status!='DELETED']|[0].status" --output text 2>/dev/null || echo "None")
if [ "${CE_STATUS}" = "DELETING" ]; then
    echo "Compute environment is DELETING, waiting..."
    while [ "$(aws batch describe-compute-environments --compute-environments "${COMPUTE_ENV}" \
        --query "computeEnvironments[?status=='DELETING']|length(@)" --output text 2>/dev/null)" != "0" ]; do
        sleep 15
    done
    CE_STATUS="None"
fi
if [ "${CE_STATUS}" = "None" ] || [ "${CE_STATUS}" = "" ]; then
    echo "Creating compute environment ${COMPUTE_ENV}..."
    if [ "${USE_SPOT}" = "spot" ]; then
        COMPUTE_TYPE="SPOT"
        ALLOC_STRATEGY="SPOT_PRICE_CAPACITY_OPTIMIZED"
        echo "  Using SPOT instances"
    else
        COMPUTE_TYPE="EC2"
        ALLOC_STRATEGY="BEST_FIT_PROGRESSIVE"
        echo "  Using ON-DEMAND instances"
    fi
    aws batch create-compute-environment \
        --compute-environment-name "${COMPUTE_ENV}" \
        --type MANAGED \
        --compute-resources '{
            "type": "'"${COMPUTE_TYPE}"'",
            "allocationStrategy": "'"${ALLOC_STRATEGY}"'",
            "minvCpus": 0,
            "maxvCpus": 8,
            "instanceTypes": ["g4dn.2xlarge"],
            "subnets": '"${SUBNETS}"',
            "securityGroupIds": ["'"${SG_ID}"'"],
            "instanceRole": "arn:aws:iam::'"${ACCOUNT_ID}"':instance-profile/'"${EC2_INSTANCE_ROLE}"'"
        }' \
        --service-role "arn:aws:iam::${ACCOUNT_ID}:role/${BATCH_SERVICE_ROLE}"
else
    echo "Compute environment ${COMPUTE_ENV} already exists."
fi

echo "Waiting for compute environment to become VALID..."
for i in $(seq 1 20); do
    STATUS=$(aws batch describe-compute-environments \
        --compute-environments "${COMPUTE_ENV}" \
        --query "computeEnvironments[0].status" --output text)
    if [ "${STATUS}" = "VALID" ]; then
        echo "Compute environment is VALID."
        break
    fi
    if [ "$i" -eq 20 ]; then
        echo "ERROR: Compute environment did not become VALID within 10 minutes."
        exit 1
    fi
    echo "  Status: ${STATUS} (attempt ${i}/20, waiting 30s...)"
    sleep 30
done

echo "=== Phase E: Job Queue ==="

if ! aws batch describe-job-queues --job-queues "${JOB_QUEUE}" \
    --query "jobQueues[?status!='DELETED']|[0].jobQueueName" --output text 2>/dev/null | grep -q "${JOB_QUEUE}"; then
    echo "Creating job queue ${JOB_QUEUE}..."
    aws batch create-job-queue \
        --job-queue-name "${JOB_QUEUE}" \
        --priority 1 \
        --compute-environment-order "order=1,computeEnvironment=${COMPUTE_ENV}"
else
    echo "Job queue ${JOB_QUEUE} already exists."
fi

echo "Waiting for job queue to become VALID..."
for i in $(seq 1 20); do
    STATUS=$(aws batch describe-job-queues \
        --job-queues "${JOB_QUEUE}" \
        --query "jobQueues[0].status" --output text)
    if [ "${STATUS}" = "VALID" ]; then
        echo "Job queue is VALID."
        break
    fi
    if [ "$i" -eq 20 ]; then
        echo "ERROR: Job queue did not become VALID within 10 minutes."
        exit 1
    fi
    echo "  Status: ${STATUS} (attempt ${i}/20, waiting 30s...)"
    sleep 30
done

echo "=== Phase F: Job Definition ==="

aws logs create-log-group --log-group-name "${LOG_GROUP}" 2>/dev/null || true

aws batch register-job-definition \
    --job-definition-name "${JOB_DEF}" \
    --type container \
    --container-properties '{
        "image": "'"${IMAGE_URI}"'",
        "resourceRequirements": [
            {"type": "VCPU", "value": "4"},
            {"type": "MEMORY", "value": "28000"},
            {"type": "GPU", "value": "1"}
        ],
        "jobRoleArn": "arn:aws:iam::'"${ACCOUNT_ID}"':role/'"${ECS_TASK_ROLE}"'",
        "executionRoleArn": "arn:aws:iam::'"${ACCOUNT_ID}"':role/'"${ECS_EXECUTION_ROLE}"'",
        "logConfiguration": {
            "logDriver": "awslogs",
            "options": {
                "awslogs-group": "'"${LOG_GROUP}"'",
                "awslogs-region": "'"${REGION}"'",
                "awslogs-stream-prefix": "batch"
            }
        }
    }' \
    --timeout '{"attemptDurationSeconds": 43200}'

echo "=== Phase G: Submit Job ==="

JOB_NAME="micropattern-integration-$(date +%Y%m%d-%H%M%S)"
JOB_ID=$(aws batch submit-job \
    --job-name "${JOB_NAME}" \
    --job-queue "${JOB_QUEUE}" \
    --job-definition "${JOB_DEF}" \
    --query "jobId" --output text)

echo ""
echo "========================================="
echo "Job submitted successfully!"
echo "========================================="
echo "Job Name: ${JOB_NAME}"
echo "Job ID:   ${JOB_ID}"
echo ""
echo "Monitor with:"
echo "  aws batch describe-jobs --jobs ${JOB_ID} --query 'jobs[0].{status:status,reason:statusReason}'"
echo "  aws logs tail ${LOG_GROUP} --follow"
echo ""
echo "Cancel with:"
echo "  aws batch cancel-job --job-id ${JOB_ID} --reason 'Manual cancellation'"
