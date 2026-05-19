#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

INPUT_FILE="${1:?Usage: $0 <input_dat> <output_dir> [extra args...] }"
OUTPUT_DIR="${2:?Usage: $0 <input_dat> <output_dir> [extra args...] }"
shift || true
shift || true
EXTRA_ARGS=("$@")

IMAGE_NAME="${IMAGE_NAME:-mass:latest}"
PLATFORM="${DOCKER_PLATFORM:-linux/amd64}"
THREADS="${THREADS:-50}"

if [[ "${INPUT_FILE}" == /* ]]; then
  INPUT_PATH="/work${INPUT_FILE}"
else
  INPUT_PATH="/work/${INPUT_FILE}"
fi

if [[ "${OUTPUT_DIR}" == /* ]]; then
  OUTPUT_PATH="/work${OUTPUT_DIR}"
else
  OUTPUT_PATH="/work/${OUTPUT_DIR}"
fi

HAS_THREADS=0
for ((i = 0; i < ${#EXTRA_ARGS[@]}; i++)); do
  arg="${EXTRA_ARGS[$i]}"
  case "${arg}" in
    -t|--threads|-t*)
      HAS_THREADS=1
      break
      ;;
    --threads=*)
      HAS_THREADS=1
      break
      ;;
  esac
done

echo "==> Running simulation..."
RUN_CMD=(docker run --rm --platform "${PLATFORM}" \
  -v "${REPO_ROOT}:/work" \
  "${IMAGE_NAME}" \
  "${INPUT_PATH}" \
  -o "${OUTPUT_PATH}")

if [[ "${HAS_THREADS}" -eq 0 ]]; then
  RUN_CMD+=(-t"${THREADS}")
fi

if (( ${#EXTRA_ARGS[@]} > 0 )); then
  RUN_CMD+=("${EXTRA_ARGS[@]}")
fi

"${RUN_CMD[@]}"
