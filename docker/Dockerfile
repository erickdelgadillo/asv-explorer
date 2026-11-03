FROM python:3.12-slim

LABEL org.opencontainers.image.title="ASV Explorer"
LABEL org.opencontainers.image.description="ASV validation runtime"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

WORKDIR /app

COPY pyproject.toml README.md ./
COPY src/ ./src/
COPY scripts/ ./scripts/

RUN python -m pip install --no-cache-dir .

WORKDIR /workspace