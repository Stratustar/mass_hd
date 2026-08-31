#!/usr/bin/env bash
set -euo pipefail

IMAGE_INPUT="${1:?Usage: $0 <image_name> [tag]}"
TAG_INPUT="${2:-}"

DOCKER_NAMESPACE="stratustar"

IMAGE_NAME="${IMAGE_INPUT}"
TAG="${TAG_INPUT}"

if [[ "${IMAGE_NAME}" == *:* ]]; then
  IMAGE_TAG_FROM_NAME="${IMAGE_NAME##*:}"
  IMAGE_NAME="${IMAGE_NAME%:*}"
  if [[ -z "${TAG}" ]]; then
    TAG="${IMAGE_TAG_FROM_NAME}"
  fi
fi

if [[ -z "${TAG}" ]]; then
  TAG="latest"
fi

LOCAL_IMAGE="${IMAGE_NAME}:${TAG}"
REMOTE_IMAGE="docker.io/${DOCKER_NAMESPACE}/${IMAGE_NAME}:${TAG}"

# "Was it built?" -- and NOT via `docker image inspect` alone.
#
# We build with --platform linux/amd64 on an arm64 Mac, so the local image is an
# amd64-only OCI index. Docker 29's containerd image store resolves `image inspect
# <name>` for the HOST platform, finds no arm64 manifest, and answers "No such image"
# for an image `image ls` lists and `push` uploads without complaint. The guard was
# therefore refusing every cross-platform build -- which is every build this project
# makes. Fall back to the tag listing, which is platform-blind.
if ! docker image inspect "${LOCAL_IMAGE}" > /dev/null 2>&1 &&
   ! docker image ls --format '{{.Repository}}:{{.Tag}}' | grep -qx "${LOCAL_IMAGE}"; then
  echo "Error: local image not found: ${LOCAL_IMAGE}" >&2
  echo "Build it first with: ./makeDocker.sh ${IMAGE_NAME}" >&2
  exit 1
fi

docker tag "${LOCAL_IMAGE}" "${REMOTE_IMAGE}"
docker push "${REMOTE_IMAGE}"

echo "Pushed ${LOCAL_IMAGE} -> ${REMOTE_IMAGE}"
