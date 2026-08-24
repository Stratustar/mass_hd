#!/usr/bin/env bash
# Submit a confluent-memory analysis script as a Slurm job.
#
# Analyses must NOT be run as `ssh jed '... python ...'`: that attaches the process to the
# ssh session, so closing the laptop sends SIGHUP and kills it mid-run (this is what
# produced the "exit code 255 / Broken pipe" failures on 2026-08-22 and 2026-08-24), and it
# puts the work on the login node. Everything goes through sbatch instead.
#
# Usage:
#   bash scripts_cluster/submit_analysis.sh <script.py> [args...]
# e.g.
#   bash scripts_cluster/submit_analysis.sh cm_phase.py 20260822/cm_regime 20260823/cm_regime_sneg \
#        --out /scratch/helu/mass_hd/results/cases/20260824/cm_phase.json
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ANALYSIS_DIR="${REPO_ROOT}/plot/python/confluent_memory"
LOG_DIR="${MASS_SCRATCH:-/scratch/helu/mass_hd}/analysis_logs"

SCRIPT="${1:?Usage: $0 <script.py> [args...]}"
shift || true
[[ "${SCRIPT}" != /* ]] && SCRIPT="${ANALYSIS_DIR}/${SCRIPT}"
[[ -f "${SCRIPT}" ]] || { echo "no such analysis script: ${SCRIPT}" >&2; exit 1; }

mkdir -p "${LOG_DIR}"
NAME="$(basename "${SCRIPT}" .py)"

sbatch --job-name="${NAME}" \
  --partition="${ANALYSIS_PARTITION:-standard}" --qos="${ANALYSIS_QOS:-serial}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${ANALYSIS_CPUS:-8}" \
  --time="${ANALYSIS_TIME:-04:00:00}" \
  --output="${LOG_DIR}/${NAME}-%j.out" --error="${LOG_DIR}/${NAME}-%j.err" \
  --wrap="cd ${ANALYSIS_DIR} && conda run --no-capture-output -n env1 python ${SCRIPT} $*"
