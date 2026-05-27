#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME_INPUT="${1:?Usage: $0 <image_name>}"

CONTAINER_DIR="/home/helu/containers"
NAMESPACE="docker.io/stratustar"
TAG="latest"

BASE_NAME="${IMAGE_NAME_INPUT##*/}"
BASE_NAME="${BASE_NAME%.sif}"

if [[ "${BASE_NAME}" == *:* ]]; then
  IMAGE_TAG="${BASE_NAME#*:}"
  BASE_NAME="${BASE_NAME%%:*}"
  if [[ -n "${IMAGE_TAG}" ]]; then
    TAG="${IMAGE_TAG}"
  fi
fi

if [[ -z "${BASE_NAME}" ]]; then
  echo "Error: invalid image name: ${IMAGE_NAME_INPUT}" >&2
  exit 1
fi

SIF_NAME="${BASE_NAME}.sif"
IMAGE_SRC="docker://${NAMESPACE}/${BASE_NAME}:${TAG}"
SIF_PATH="${CONTAINER_DIR}/${SIF_NAME}"

mkdir -p "${CONTAINER_DIR}"

echo "Pulling ${IMAGE_SRC} -> ${SIF_PATH}"
apptainer pull --force "${SIF_PATH}" "${IMAGE_SRC}"

echo "Done: ${SIF_PATH}"
