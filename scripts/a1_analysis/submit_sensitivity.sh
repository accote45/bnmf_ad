#!/bin/bash
# submit_sensitivity.sh
# LSF submission for A1 sensitivity clustering analysis (META).
# Runs the k-means / hierarchical / half-max / DynamicTreeCut comparison
# (03) followed by the bNMF-vs-hierarchical logistic comparison (08), which
# depends on 03's outputs. Both run in one job so 08 sees 03's results.
#
# Usage:
#   bash scripts/a1_analysis/submit_sensitivity.sh [ANCESTRY]   # default META

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ANCESTRY="${1:-META}"
LOG_DIR="${PROJECT_ROOT}/logs/a1_sensitivity"

mkdir -p "${LOG_DIR}"

echo "=== A1 Sensitivity Clustering Submission ==="
echo "  Project root: ${PROJECT_ROOT}"
echo "  Ancestry:     ${ANCESTRY}"
echo ""

JOB_OUT=$(bsub -P acc_paul_oreilly -J "a1_sens_${ANCESTRY}" \
  -q premium -n 1 -R "rusage[mem=24000]" -W 3:00 \
  -o "${LOG_DIR}/sensitivity_${ANCESTRY}.out" \
  -e "${LOG_DIR}/sensitivity_${ANCESTRY}.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0 gcc/14.2.0; \
   Rscript scripts/a1_analysis/03_sensitivity_clustering.R --config config/a1_config.yaml --ancestry ${ANCESTRY} && \
   Rscript scripts/a1_analysis/08_bnmf_vs_hclust_comparison.R --config config/a1_config.yaml --ancestry ${ANCESTRY}")
JOB_ID=$(echo "${JOB_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB_ID}"
echo ""
echo "Monitor with: bjobs ${JOB_ID}"
echo "Output log:   ${LOG_DIR}/sensitivity_${ANCESTRY}.out"
echo "Error log:    ${LOG_DIR}/sensitivity_${ANCESTRY}.err"
