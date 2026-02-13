FROM ghcr.io/astral-sh/uv:latest AS uv

FROM nvidia/cuda:12.8.1-runtime-ubuntu24.04

COPY --from=uv /uv /uvx /bin/

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.12 python3.12-venv python3.12-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml uv.lock README.md ./
COPY src/ src/
COPY scripts/run_integration_multi_aws.py scripts/

RUN uv sync --frozen --no-dev

CMD ["uv", "run", "python", "scripts/run_integration_multi_aws.py"]
