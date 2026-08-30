#!/usr/bin/env bash
# Submit every run.dat under a cases dir as ONE Slurm job array.
#
# Why an array and not submit_many.sh's loop: a 248-case scan has to hand its results to an
# analysis stage without anyone sitting at a terminal, and dependencies are what make that
# work. With individual jobs the analysis has to name 248 job ids, which Slurm caps
# (MaxJobDependency, commonly 128) -- with an array it names one, and `aftercorr` additionally
# lets task i of a follow-on array start the moment task i of this one succeeds, instead of
# waiting for the slowest case in the batch.
#
# Each task just calls submit_case.sh as a plain script, so the output-path derivation, the
# stale-file cleanup and the in-job plotting all stay in one place.
#
# Usage: submit_array.sh <sif_name> <cases_dir>
# Env:   CONCURRENCY (default 20), ARRAY_TIME (03:00:00), ARRAY_CPUS (32),
#        plus anything submit_case.sh reads (SKIP_PLOTS, PLOT_SCRIPT, ...)
set -euo pipefail

SIF_NAME="${1:?Usage: $0 <sif_name> <cases_dir>}"
CASES_DIR="${2:?Usage: $0 <sif_name> <cases_dir>}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="/home/helu/mass_hd"
LIST_ROOT="/scratch/helu/mass_hd/joblists"
LOG_ROOT="/scratch/helu/mass_hd/analysis_logs"
CONCURRENCY="${CONCURRENCY:-20}"
ARRAY_TIME="${ARRAY_TIME:-03:00:00}"
ARRAY_CPUS="${ARRAY_CPUS:-32}"

if [[ "${CASES_DIR}" == /* ]]; then SEARCH_ROOT="${CASES_DIR}"; else SEARCH_ROOT="${REPO_ROOT}/${CASES_DIR}"; fi
[[ -d "${SEARCH_ROOT}" ]] || { echo "Error: cases dir not found: ${SEARCH_ROOT}" >&2; exit 1; }

mkdir -p "${LIST_ROOT}" "${LOG_ROOT}"
TAG="$(echo "${CASES_DIR}" | tr '/' '_')"
LIST="${LIST_ROOT}/${TAG}.txt"
# The list must live on a SHARED filesystem: /tmp on the login node is not the /tmp a compute
# node sees, and a job that reads its task list from there silently gets an empty path.
find "${SEARCH_ROOT}" -type f -name '*.dat' | sort > "${LIST}"
N="$(wc -l < "${LIST}")"
[[ "${N}" -gt 0 ]] || { echo "No .dat files under ${SEARCH_ROOT}" >&2; exit 1; }

echo "Array: ${N} cases from ${SEARCH_ROOT}"
echo "  list        ${LIST}"
echo "  concurrency ${CONCURRENCY}, time ${ARRAY_TIME}, cpus ${ARRAY_CPUS}"

JOBID=$(sbatch --parsable --export=ALL \
  --job-name="arr_${TAG}" \
  --array="0-$((N-1))%${CONCURRENCY}" \
  --partition=standard --qos=serial --nodes=1 --ntasks=1 \
  --cpus-per-task="${ARRAY_CPUS}" --time="${ARRAY_TIME}" \
  --output="${LOG_ROOT}/${TAG}-%A_%a.out" \
  --error="${LOG_ROOT}/${TAG}-%A_%a.err" \
  --wrap="DAT=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID+1))p\" ${LIST}); \
          echo \"task \$SLURM_ARRAY_TASK_ID -> \$DAT\"; \
          bash ${SCRIPT_DIR}/submit_case.sh ${SIF_NAME} \"\$DAT\"")

echo "${JOBID}"
