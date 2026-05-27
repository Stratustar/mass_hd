#!/usr/bin/env bash
#SBATCH --job-name=mass
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=standard
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

set -euo pipefail

SIF_NAME="${1:?Usage: $0 <sif_name> <input_dat> [output_dir]}"
INPUT_ARG="${2:?Usage: $0 <sif_name> <input_dat> [output_dir]}"
OUTPUT_ARG="${3-}"

# ----------------------------- cluster paths -----------------------------
CONTAINER_DIR="/home/helu/containers"
REPO_ROOT="/home/helu/mass_hd"
SCRATCH_ROOT="/scratch/helu"
CONDA_BIN="/home/helu/miniconda3/bin/conda"
CONDA_ENV="env1"
PLOT_SCRIPT="${REPO_ROOT}/plot/python/plot_hd.py"
PROLIFERATION_SUMMARY_SCRIPT="${REPO_ROOT}/plot/python/proliferation_summary.py"
MASS_BINARY="/opt/mass_hd/mass"
RESULTS_ROOT="${REPO_ROOT}/results"
THREADS="${SLURM_CPUS_PER_TASK}"
SKIP_PLOTS="${SKIP_PLOTS:-0}"
SKIP_SUMMARY="${SKIP_SUMMARY:-0}"
PLOT_HD_ARGS="${PLOT_HD_ARGS:-}"
# ----------------------------------------------------------------------

if [[ "${SIF_NAME}" == *.sif ]]; then
  SIF_FILE="${SIF_NAME}"
else
  SIF_FILE="${SIF_NAME}.sif"
fi
IMAGE_PATH="${CONTAINER_DIR}/${SIF_FILE}"

if [[ "${INPUT_ARG}" == /* ]]; then
  INPUT_PATH="${INPUT_ARG}"
else
  INPUT_PATH="${REPO_ROOT}/${INPUT_ARG}"
fi

INPUT_DIR="$(dirname "${INPUT_PATH}")"
INPUT_BASENAME="$(basename "${INPUT_PATH}")"

if [[ -n "${OUTPUT_ARG}" ]]; then
  if [[ "${OUTPUT_ARG}" == /* ]]; then
    OUTPUT_PATH="${OUTPUT_ARG}"
  else
    OUTPUT_PATH="${SCRATCH_ROOT}/${OUTPUT_ARG}"
  fi
else
  OUTPUT_PATH="${SCRATCH_ROOT}/${INPUT_DIR#/home/helu/}"
fi

if [[ ! -f "${IMAGE_PATH}" ]]; then
  echo "Error: image not found: ${IMAGE_PATH}" >&2
  exit 1
fi

if [[ ! -f "${INPUT_PATH}" ]]; then
  echo "Error: input dat not found: ${INPUT_PATH}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_PATH}"

# Plot folder mirrors the project-relative scratch output path.
if [[ "${OUTPUT_PATH}" == ${SCRATCH_ROOT}/* ]]; then
  REL_OUTPUT_DIR="${OUTPUT_PATH#${SCRATCH_ROOT}/}"
  REL_OUTPUT_DIR="${REL_OUTPUT_DIR#mass_hd/}"
else
  REL_OUTPUT_DIR="$(basename "${OUTPUT_PATH}")"
fi
PLOT_DIR="${RESULTS_ROOT}/${REL_OUTPUT_DIR}"
mkdir -p "${PLOT_DIR}"

export OMP_NUM_THREADS="${THREADS}"
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

find "${OUTPUT_PATH}" -maxdepth 1 -type f -name '*.json' -delete
find "${PLOT_DIR}" -maxdepth 1 -type f \( -name '*.png' -o -name '*.gif' -o -name '*.csv' \) -delete

echo "Job started on $(date)"
echo "Image:       ${IMAGE_PATH}"
echo "Input:       ${INPUT_PATH}"
echo "Sim output:  ${OUTPUT_PATH}"
echo "Plot output: ${PLOT_DIR}"
echo "Threads:     ${THREADS}"
echo "Skip plots:  ${SKIP_PLOTS}"
echo "Skip summary:${SKIP_SUMMARY}"
echo "Plot args:   ${PLOT_HD_ARGS:-<none>}"

apptainer exec \
  --bind "${INPUT_DIR}":/input \
  --bind "${OUTPUT_PATH}":/output \
  "${IMAGE_PATH}" \
  "${MASS_BINARY}" \
  /input/${INPUT_BASENAME} -o /output -t "${THREADS}"

SIM_STATUS=$?
if [[ ${SIM_STATUS} -ne 0 ]]; then
  echo "Simulation failed with status ${SIM_STATUS} on $(date)"
  exit ${SIM_STATUS}
fi

echo "Simulation finished on $(date)"

case "${SKIP_PLOTS}" in
  1|true|TRUE|yes|YES)
    echo "Skipping plotting on request."
    exit 0
    ;;
esac

echo "Running plotting..."

PLOT_HD_ARGV=()
if [[ -n "${PLOT_HD_ARGS}" ]]; then
  read -r -a PLOT_HD_ARGV <<< "${PLOT_HD_ARGS}"
fi

"${CONDA_BIN}" run --no-capture-output -n "${CONDA_ENV}" \
  python "${PLOT_SCRIPT}" "${OUTPUT_PATH}" "${PLOT_DIR}" "${PLOT_HD_ARGV[@]}"

PLOT_STATUS=$?
if [[ ${PLOT_STATUS} -ne 0 ]]; then
  echo "Plotting failed with status ${PLOT_STATUS} on $(date)"
  exit ${PLOT_STATUS}
fi

case "${SKIP_SUMMARY}" in
  1|true|TRUE|yes|YES)
    echo "Skipping proliferation summary on request."
    echo "Plotting finished on $(date)"
    exit 0
    ;;
esac

echo "Running proliferation summary..."

"${CONDA_BIN}" run --no-capture-output -n "${CONDA_ENV}" \
  python "${PROLIFERATION_SUMMARY_SCRIPT}" "${OUTPUT_PATH}" "${PLOT_DIR}"

SUMMARY_STATUS=$?
if [[ ${SUMMARY_STATUS} -ne 0 ]]; then
  echo "Proliferation summary failed with status ${SUMMARY_STATUS} on $(date)"
  exit ${SUMMARY_STATUS}
fi

echo "Plotting finished on $(date)"
exit 0
