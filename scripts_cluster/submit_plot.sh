#!/usr/bin/env bash
#SBATCH --job-name=massplot
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=standard
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Plot-only job: regenerate figures/analysis from an EXISTING simulation archive,
# without re-running the simulation. Mirrors the plotting half of submit_case.sh.
# Usage: sbatch submit_plot.sh <input_dat> [output_dir]

set -euo pipefail

INPUT_ARG="${1:?Usage: $0 <input_dat> [output_dir]}"
OUTPUT_ARG="${2-}"

REPO_ROOT="/home/helu/mass_hd"
SCRATCH_ROOT="/scratch/helu"
CONDA_BIN="/home/helu/miniconda3/bin/conda"
CONDA_ENV="env1"
PLOT_SCRIPT="${REPO_ROOT}/plot/python/plot_hd.py"
RESULTS_ROOT="${REPO_ROOT}/results"
THREADS="${SLURM_CPUS_PER_TASK}"
PLOT_HD_ARGS="${PLOT_HD_ARGS:-}"

if [[ "${INPUT_ARG}" == /* ]]; then
  INPUT_PATH="${INPUT_ARG}"
else
  INPUT_PATH="${REPO_ROOT}/${INPUT_ARG}"
fi
INPUT_DIR="$(dirname "${INPUT_PATH}")"

if [[ -n "${OUTPUT_ARG}" ]]; then
  if [[ "${OUTPUT_ARG}" == /* ]]; then
    OUTPUT_PATH="${OUTPUT_ARG}"
  else
    OUTPUT_PATH="${SCRATCH_ROOT}/${OUTPUT_ARG}"
  fi
else
  OUTPUT_PATH="${SCRATCH_ROOT}/${INPUT_DIR#/home/helu/}"
fi

# Plot folder mirrors the project-relative scratch output path.
if [[ "${OUTPUT_PATH}" == ${SCRATCH_ROOT}/* ]]; then
  REL_OUTPUT_DIR="${OUTPUT_PATH#${SCRATCH_ROOT}/}"
  REL_OUTPUT_DIR="${REL_OUTPUT_DIR#mass_hd/}"
else
  REL_OUTPUT_DIR="$(basename "${OUTPUT_PATH}")"
fi
PLOT_DIR="${RESULTS_ROOT}/${REL_OUTPUT_DIR}"
mkdir -p "${PLOT_DIR}"

if ! ls "${OUTPUT_PATH}"/*.json >/dev/null 2>&1; then
  echo "Error: no frames (*.json) in ${OUTPUT_PATH}" >&2
  exit 1
fi

export OMP_NUM_THREADS="${THREADS}"
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

# Remove stale figures so we regenerate cleanly.
find "${PLOT_DIR}" -maxdepth 1 -type f \( -name '*.png' -o -name '*.gif' -o -name '*.csv' \) -delete

echo "Plot-only job started on $(date)"
echo "Sim output:  ${OUTPUT_PATH}"
echo "Plot output: ${PLOT_DIR}"
echo "Threads:     ${THREADS}"

PLOT_HD_ARGV=()
if [[ -n "${PLOT_HD_ARGS}" ]]; then
  read -r -a PLOT_HD_ARGV <<< "${PLOT_HD_ARGS}"
fi

"${CONDA_BIN}" run --no-capture-output -n "${CONDA_ENV}" \
  python "${PLOT_SCRIPT}" "${OUTPUT_PATH}" "${PLOT_DIR}" "${PLOT_HD_ARGV[@]}"

echo "Plotting finished on $(date)"
exit 0
