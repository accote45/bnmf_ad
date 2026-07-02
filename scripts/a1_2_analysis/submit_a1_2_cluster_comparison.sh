#!/bin/bash
# submit_a1_2_cluster_comparison.sh
# LSF submission for A1.2: Compare META bNMF clusters against published studies
#
# Usage:
#   bash scripts/a1_2_analysis/submit_a1_2_cluster_comparison.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
LOG_DIR="${PROJECT_ROOT}/logs/a1_2"

mkdir -p "${LOG_DIR}"

echo "=== A1.2 Cluster Comparison Submission ==="
echo "  Project root: ${PROJECT_ROOT}"
echo ""

JOB_OUT=$(bsub -P acc_paul_oreilly -J a1_2_cluster_comp \
  -q premium -n 1 -R "rusage[mem=16000]" -W 2:00 \
  -o "${LOG_DIR}/a1_2_cluster_comparison.out" \
  -e "${LOG_DIR}/a1_2_cluster_comparison.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0 gcc/14.2.0; \
   Rscript scripts/a1_2_analysis/a1_2_cluster_comparison.R")
JOB_ID=$(echo "${JOB_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB_ID}"
echo ""
echo "Monitor with: bjobs ${JOB_ID}"
echo "Output log:   ${LOG_DIR}/a1_2_cluster_comparison.out"
echo "Error log:    ${LOG_DIR}/a1_2_cluster_comparison.err"
