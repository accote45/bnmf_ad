#!/bin/bash
# submit_b1_2_pipeline.sh
# LSF submission script for B1.2 Pathway PRS vs Cluster PRS Comparison.
# Dependency chaining: METAL prep -> METAL -> PRSet prep -> PRSet -> Association -> Compare
#
# Usage:
#   bash scripts/b1_2_analysis/submit_b1_2_pipeline.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CONFIG="${PROJECT_ROOT}/config/b1_2_config.yaml"
LOG_DIR="${PROJECT_ROOT}/logs/b1_2"

mkdir -p "${LOG_DIR}"

echo "=== B1.2 Pathway PRS Pipeline Submission ==="
echo "  Project root: ${PROJECT_ROOT}"
echo "  Config: ${CONFIG}"
echo ""

# =====================================================================
# Step 1: Prepare METAL input
# =====================================================================

echo "--- Step 1: Submitting METAL input preparation ---"
JOB1_OUT=$(bsub -P acc_paul_oreilly -J b1_2_metal_prep \
  -q express -n 1 -R "rusage[mem=40000]" -W 1:00 \
  -o "${LOG_DIR}/01_metal_prep.out" -e "${LOG_DIR}/01_metal_prep.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   Rscript scripts/b1_2_analysis/01_prepare_metal_input.R --config ${CONFIG}")
JOB1_ID=$(echo "${JOB1_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB1_ID}"

# =====================================================================
# Step 2: Run METAL
# =====================================================================

echo "--- Step 2: Submitting METAL ---"
JOB2_OUT=$(bsub -P acc_paul_oreilly -J b1_2_metal \
  -w "done(${JOB1_ID})" \
  -q express -n 1 -R "rusage[mem=32000]" -W 0:30 \
  -o "${LOG_DIR}/02_metal.out" -e "${LOG_DIR}/02_metal.err" \
  "cd ${PROJECT_ROOT}; bash scripts/b1_2_analysis/02_run_metal.sh")
JOB2_ID=$(echo "${JOB2_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB2_ID}"

# =====================================================================
# Step 3: Prepare PRSet input
# =====================================================================

echo "--- Step 3: Submitting PRSet input preparation ---"
JOB3_OUT=$(bsub -P acc_paul_oreilly -J b1_2_prset_prep \
  -w "done(${JOB2_ID})" \
  -q express -n 1 -R "rusage[mem=20000]" -W 0:30 \
  -o "${LOG_DIR}/03_prset_prep.out" -e "${LOG_DIR}/03_prset_prep.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   Rscript scripts/b1_2_analysis/03_prepare_prset_input.R --config ${CONFIG}")
JOB3_ID=$(echo "${JOB3_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB3_ID}"

# =====================================================================
# Step 4: Run PRSet (long, resource-intensive)
# =====================================================================

echo "--- Step 4: Submitting PRSet ---"
JOB4_OUT=$(bsub -P acc_paul_oreilly -J b1_2_prset \
  -w "done(${JOB3_ID})" \
  -q premium -n 4 -R "rusage[mem=48000]" -W 144:00 \
  -o "${LOG_DIR}/04_prset.out" -e "${LOG_DIR}/04_prset.err" \
  "cd ${PROJECT_ROOT}; bash scripts/b1_2_analysis/04_run_prset.sh")
JOB4_ID=$(echo "${JOB4_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB4_ID}"

# =====================================================================
# Step 5: Pathway PRS association testing
# =====================================================================

echo "--- Step 5: Submitting pathway association testing ---"
JOB5_OUT=$(bsub -P acc_paul_oreilly -J b1_2_assoc \
  -w "done(${JOB4_ID})" \
  -q premium -n 16 -R "rusage[mem=2000]" -W 4:00 \
  -o "${LOG_DIR}/05_assoc.out" -e "${LOG_DIR}/05_assoc.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   OPENBLAS_NUM_THREADS=1 Rscript scripts/b1_2_analysis/05_pathway_prs_association.R --config ${CONFIG}")
JOB5_ID=$(echo "${JOB5_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB5_ID}"

# =====================================================================
# Step 6: Comparison visualization
# =====================================================================

echo "--- Step 6: Submitting comparison visualization ---"
JOB6_OUT=$(bsub -P acc_paul_oreilly -J b1_2_compare \
  -w "done(${JOB5_ID})" \
  -q express -n 1 -R "rusage[mem=4000]" -W 0:30 \
  -o "${LOG_DIR}/06_compare.out" -e "${LOG_DIR}/06_compare.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   Rscript scripts/b1_2_analysis/06_compare_pathway_vs_cluster.R --config ${CONFIG}")
JOB6_ID=$(echo "${JOB6_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB6_ID}"

# =====================================================================
# Summary
# =====================================================================

echo ""
echo "=== All jobs submitted ==="
echo ""
echo "  Dependency chain:"
echo "    Step 1: b1_2_metal_prep (${JOB1_ID})"
echo "      -> Step 2: b1_2_metal (${JOB2_ID})"
echo "           -> Step 3: b1_2_prset_prep (${JOB3_ID})"
echo "                -> Step 4: b1_2_prset (${JOB4_ID}) [LONG - premium queue, 16 cores]"
echo "                     -> Step 5: b1_2_assoc (${JOB5_ID})"
echo "                          -> Step 6: b1_2_compare (${JOB6_ID})"
echo ""
echo "  Monitor with: bjobs -w | grep b1_2"
echo "  Logs in: ${LOG_DIR}/"
echo "  Results in: ${PROJECT_ROOT}/results/b1_2_analysis/"
