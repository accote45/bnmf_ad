#!/bin/bash
# submit_b2_1_pipeline.sh
# LSF submission for B2.1: Cox models + cumulative hazard plots
# for META cluster PRS at 3 thresholds.
#
# Usage:
#   bash scripts/b2_analysis/submit_b2_1_pipeline.sh
#   bash scripts/b2_analysis/submit_b2_1_pipeline.sh --hard-assignment

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CONFIG="${PROJECT_ROOT}/config/b2_1_config.yaml"

HARD_ASSIGNMENT=""
SUFFIX=""
if [[ "${1:-}" == "--hard-assignment" ]]; then
  HARD_ASSIGNMENT="--hard-assignment"
  SUFFIX="_hard_assignment"
  echo "  ** HARD ASSIGNMENT MODE **"
fi

LOG_DIR="${PROJECT_ROOT}/logs/b2_1${SUFFIX}"

mkdir -p "${LOG_DIR}"

echo "=== B2.1 Cluster PRS Survival Analysis ==="
echo "  Project root: ${PROJECT_ROOT}"
echo "  Config: ${CONFIG}"
echo ""

echo "--- Submitting b2_1_analysis${SUFFIX} ---"
JOB_OUT=$(bsub -P acc_paul_oreilly -J "b2_1_analysis${SUFFIX}" \
  -q express -n 1 -R "rusage[mem=16000]" -W 1:00 \
  -o "${LOG_DIR}/b2_1_analysis.out" -e "${LOG_DIR}/b2_1_analysis.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   Rscript scripts/b2_analysis/b2_1_analysis.R --config ${CONFIG} ${HARD_ASSIGNMENT}")
JOB_ID=$(echo "${JOB_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB_ID}"

echo ""
echo "=== Job submitted ==="
echo "  Monitor with: bjobs ${JOB_ID}"
echo "  Logs in: ${LOG_DIR}/"
echo "  Results in: ${PROJECT_ROOT}/results/b2_1_analysis${SUFFIX}/"
