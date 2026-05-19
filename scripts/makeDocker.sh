#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

IMAGE_BASE="${1:-mass_hd}"
IMAGE_NAME="${IMAGE_BASE}:latest"
PLATFORM="${DOCKER_PLATFORM:-linux/amd64}"

if [[ -z "${IMAGE_BASE}" ]]; then
  echo "Error: image name is empty" >&2
  exit 1
fi

echo "==> Building Docker image: ${IMAGE_NAME}"
docker build --platform "${PLATFORM}" -t "${IMAGE_NAME}" "${PROJECT_ROOT}"
