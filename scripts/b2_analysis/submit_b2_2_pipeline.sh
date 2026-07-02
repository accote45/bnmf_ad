#!/bin/bash
# submit_b2_2_pipeline.sh
# LSF submission for B2.2 PheWAS: logistic regression of META cluster PRS
# vs ICD10 binary phenotypes + Manhattan plots.
#
# Usage:
#   bash scripts/b2_analysis/submit_b2_2_pipeline.sh
#   bash scripts/b2_analysis/submit_b2_2_pipeline.sh --hard-assignment

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CONFIG="${PROJECT_ROOT}/config/b2_2_config.yaml"

HARD_ASSIGNMENT=""
SUFFIX=""
if [[ "${1:-}" == "--hard-assignment" ]]; then
  HARD_ASSIGNMENT="--hard-assignment"
  SUFFIX="_hard_assignment"
  echo "  ** HARD ASSIGNMENT MODE **"
fi

LOG_DIR="${PROJECT_ROOT}/logs/b2_2${SUFFIX}"

mkdir -p "${LOG_DIR}"

echo "=== B2.2 PheWAS Submission ==="
echo "  Project root: ${PROJECT_ROOT}"
echo "  Config: ${CONFIG}"
echo ""

echo "--- Submitting b2_2_phewas${SUFFIX} ---"
JOB_OUT=$(bsub -P acc_paul_oreilly -J "b2_2_phewas${SUFFIX}" \
  -q premium -n 4 -R "rusage[mem=32000] span[hosts=1]" -W 24:00 \
  -o "${LOG_DIR}/b2_2_phewas.out" -e "${LOG_DIR}/b2_2_phewas.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   Rscript scripts/b2_analysis/b2_2_phewas.R --config ${CONFIG} ${HARD_ASSIGNMENT} && \
   Rscript scripts/b2_analysis/b2_2_manhattan_plots.R ${HARD_ASSIGNMENT}")
JOB_ID=$(echo "${JOB_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB_ID}"

echo ""
echo "=== Job submitted ==="
echo "  Monitor with: bjobs ${JOB_ID}"
echo "  Logs in: ${LOG_DIR}/"
echo "  Results in: ${PROJECT_ROOT}/results/b2_2_analysis${SUFFIX:+/prs_hard_assignment}/"
