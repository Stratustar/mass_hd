#!/usr/bin/env bash
# Run an analysis script over every case of a scan, as one Slurm array, optionally chained
# behind the simulation array that produces the data.
#
# The point is unattended operation: with DEPEND=aftercorr:<simarray> task i starts the moment
# simulation task i succeeds, so per-case analysis overlaps the rest of the batch instead of
# waiting for the slowest case. The task list is the SAME file submit_array.sh wrote, so the
# two arrays agree on which index means which case by construction rather than by both
# re-globbing the directory and hoping the order matches.
#
# The analysis script is called as:  <script> <case_input_dir> <case_output_dir> [extra args]
#
# Usage: submit_analysis_array.sh <script_relpath> <cases_dir> <out_subdir> [extra args...]
# Env:   DEPEND (e.g. aftercorr:12345), CONCURRENCY (20), ARRAY_TIME (02:00:00),
#        ARRAY_CPUS (4), NAME (defaults to the script's basename)
set -euo pipefail

SCRIPT_REL="${1:?Usage: $0 <script_relpath> <cases_dir> <out_subdir> [extra args]}"
CASES_DIR="${2:?Usage: $0 <script_relpath> <cases_dir> <out_subdir> [extra args]}"
OUT_SUBDIR="${3:?Usage: $0 <script_relpath> <cases_dir> <out_subdir> [extra args]}"
shift 3
EXTRA="$*"

REPO_ROOT="/home/helu/mass_hd"
SCRATCH="/scratch/helu/mass_hd"
LOG_ROOT="${SCRATCH}/analysis_logs"
CONDA="/home/helu/miniconda3/bin/conda"
CONCURRENCY="${CONCURRENCY:-20}"
ARRAY_TIME="${ARRAY_TIME:-02:00:00}"
ARRAY_CPUS="${ARRAY_CPUS:-4}"
NAME="${NAME:-$(basename "${SCRIPT_REL}" .py)}"

TAG="$(echo "${CASES_DIR}" | tr '/' '_')"
LIST="${SCRATCH}/joblists/${TAG}.txt"
[[ -f "${LIST}" ]] || { echo "Error: no task list at ${LIST}; run submit_array.sh first" >&2; exit 1; }
N="$(wc -l < "${LIST}")"
mkdir -p "${LOG_ROOT}"

DEP_ARG=()
[[ -n "${DEPEND:-}" ]] && DEP_ARG=(--dependency="${DEPEND}")

echo "Analysis array: ${NAME} over ${N} cases of ${CASES_DIR}"
[[ -n "${DEPEND:-}" ]] && echo "  chained behind ${DEPEND}"

sbatch --parsable --export=ALL "${DEP_ARG[@]}" \
  --job-name="${NAME}" \
  --array="0-$((N-1))%${CONCURRENCY}" \
  --partition=standard --qos=serial --nodes=1 --ntasks=1 \
  --cpus-per-task="${ARRAY_CPUS}" --time="${ARRAY_TIME}" \
  --output="${LOG_ROOT}/${NAME}-%A_%a.out" \
  --error="${LOG_ROOT}/${NAME}-%A_%a.err" \
  --wrap="set -e; \
    DAT=\$(sed -n \"\$((SLURM_ARRAY_TASK_ID+1))p\" ${LIST}); \
    CASE=\$(basename \$(dirname \"\$DAT\")); \
    REL=\$(dirname \"\$DAT\"); REL=\${REL#${REPO_ROOT}/cases/}; \
    IN=${SCRATCH}/cases/\${REL}; \
    OUT=${SCRATCH}/results/cases/\${REL}/${OUT_SUBDIR}; \
    echo \"task \$SLURM_ARRAY_TASK_ID  case \$CASE\"; \
    mkdir -p \"\$OUT\"; \
    ${CONDA} run --no-capture-output -n env1 python ${REPO_ROOT}/${SCRIPT_REL} \"\$IN\" \"\$OUT\" ${EXTRA}"
