#!/usr/bin/env bash
set -euo pipefail

SIF_NAME="${1:?Usage: $0 <sif_name> <cases_dir>}"
CASES_DIR="${2:?Usage: $0 <sif_name> <cases_dir>}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="/home/helu/mass_hd"

if [[ "${CASES_DIR}" == /* ]]; then
  SEARCH_ROOT="${CASES_DIR}"
else
  SEARCH_ROOT="${REPO_ROOT}/${CASES_DIR}"
fi

if [[ ! -d "${SEARCH_ROOT}" ]]; then
  echo "Error: cases directory not found: ${SEARCH_ROOT}" >&2
  exit 1
fi

count=0
while IFS= read -r -d '' dat_file; do
  echo "Submitting: ${dat_file}"
  sbatch "${SCRIPT_DIR}/submit_case.sh" "${SIF_NAME}" "${dat_file}"
  ((count += 1))
done < <(find "${SEARCH_ROOT}" -type f -name '*.dat' -print0 | sort -z)

if [[ ${count} -eq 0 ]]; then
  echo "No .dat files found in ${SEARCH_ROOT}"
  exit 0
fi

echo "Submitted ${count} jobs."
