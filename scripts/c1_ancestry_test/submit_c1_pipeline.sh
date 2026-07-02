#!/bin/bash
# submit_c1_pipeline.sh
# Master LSF submission script for C1.1 Empiric Ancestry Difference Test.
# Orchestrates all steps with dependency chaining.
#
# Usage:
#   bash scripts/c1_ancestry_test/submit_c1_pipeline.sh
#
# Prerequisites:
#   - Fill in PHENO_FILE and COVAR_FILE below (or in c1_config.yaml)
#   - Ensure 1KG reference panels are downloaded (scripts/gwas_processing/download_1kg_reference.sh)
#   - Ensure published harmonized GWAS files exist in sumstats/harmonized/

set -euo pipefail

# --- Configuration ---
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CONFIG="${PROJECT_ROOT}/config/c1_config.yaml"
LOG_DIR="${PROJECT_ROOT}/logs/c1"
RESULTS_DIR="${PROJECT_ROOT}/results/c1_ancestry_test"

# UKB extracted phenotype/covariate files (produced by 00_extract_phenotypes.R)
export PHENO_FILE="${RESULTS_DIR}/phenotypes.txt"
export COVAR_FILE="${RESULTS_DIR}/covariates.txt"
export PROJECT_ROOT

mkdir -p "${LOG_DIR}"
mkdir -p "${RESULTS_DIR}"

echo "=== C1.1 Pipeline Submission ==="
echo "  Project root: ${PROJECT_ROOT}"
echo "  Config: ${CONFIG}"
echo "  Phenotype file: ${PHENO_FILE}"
echo "  Covariate file: ${COVAR_FILE}"
echo ""

# --- Validate configuration ---
if [ ! -f "${PHENO_FILE}" ]; then
  echo "ERROR: Phenotype file not found: ${PHENO_FILE}"
  echo "Run 00_extract_phenotypes.R first to extract phenotypes from UKB database."
  exit 1
fi
if [ ! -f "${COVAR_FILE}" ]; then
  echo "ERROR: Covariate file not found: ${COVAR_FILE}"
  echo "Run 00_extract_phenotypes.R first to extract covariates from UKB database."
  exit 1
fi

# =====================================================================
# Step 0: Subsample EUR and AFR individuals
# =====================================================================

echo "--- Step 0: Submitting subsampling job ---"
JOB0_OUT=$(bsub -P acc_paul_oreilly -J c1_subsample \
  -q express -n 1 -R "rusage[mem=8000]" -W 0:30 \
  -o "${LOG_DIR}/subsample.out" -e "${LOG_DIR}/subsample.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   Rscript scripts/c1_ancestry_test/00_subsample_ukb.R --config ${CONFIG}")
JOB0_ID=$(echo "${JOB0_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB0_ID}"

# =====================================================================
# Step 1: Run GWAS (array job: 33 elements)
# =====================================================================

echo "--- Step 1: Submitting GWAS array job ---"
JOB1_OUT=$(bsub -P acc_paul_oreilly -J "c1_gwas[1-18]" \
  -w "done(${JOB0_ID})" \
  -q express -n 4 -R "rusage[mem=16000]" -W 4:00 \
  -o "${LOG_DIR}/gwas_%I.out" -e "${LOG_DIR}/gwas_%I.err" \
  -env "all, PROJECT_ROOT=${PROJECT_ROOT}, PHENO_FILE=${PHENO_FILE}, COVAR_FILE=${COVAR_FILE}, N_SUBS=5" \
  "cd ${PROJECT_ROOT}; module load plink2 R/4.2.0; \
   bash scripts/c1_ancestry_test/01_run_gwas.sh")
JOB1_ID=$(echo "${JOB1_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB1_ID} (array 1-18)"

# =====================================================================
# Step 2: Run bNMF for each EUR subsample (array job: 10 elements)
# =====================================================================

echo "--- Step 2: Submitting EUR bNMF array job ---"
JOB2_OUT=$(bsub -P acc_paul_oreilly -J "c1_bnmf_eur[1-5]" \
  -w "done(${JOB1_ID})" \
  -q express -n 1 -R "rusage[mem=32000]" -W 12:00 \
  -o "${LOG_DIR}/bnmf_eur_%I.out" -e "${LOG_DIR}/bnmf_eur_%I.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0 plink/1.90b6.21; \
   Rscript scripts/c1_ancestry_test/02_run_bnmf_subsample.R \
     --subsample-id \$LSB_JOBINDEX --config ${CONFIG}")
JOB2_ID=$(echo "${JOB2_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB2_ID} (array 1-5)"

# =====================================================================
# Step 3: Run bNMF for AFR (single job)
# =====================================================================

echo "--- Step 3: Submitting AFR bNMF job ---"
JOB3_OUT=$(bsub -P acc_paul_oreilly -J c1_bnmf_afr \
  -w "done(${JOB1_ID})" \
  -q express -n 1 -R "rusage[mem=32000]" -W 12:00 \
  -o "${LOG_DIR}/bnmf_afr.out" -e "${LOG_DIR}/bnmf_afr.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0 plink/1.90b6.21; \
   Rscript scripts/c1_ancestry_test/03_run_bnmf_afr.R --config ${CONFIG}")
JOB3_ID=$(echo "${JOB3_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB3_ID}"

# =====================================================================
# Step 4: Compute Jaccard distances + Hungarian matching
# =====================================================================

echo "--- Step 4: Submitting Jaccard analysis job ---"
JOB4_OUT=$(bsub -P acc_paul_oreilly -J c1_jaccard \
  -w "done(${JOB2_ID}) && done(${JOB3_ID})" \
  -q express -n 1 -R "rusage[mem=16000]" -W 1:00 \
  -o "${LOG_DIR}/jaccard.out" -e "${LOG_DIR}/jaccard.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   Rscript scripts/c1_ancestry_test/04_compute_jaccard.R --config ${CONFIG}")
JOB4_ID=$(echo "${JOB4_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB4_ID}"

# =====================================================================
# Step 5: Generate visualizations
# =====================================================================

echo "--- Step 5: Submitting visualization job ---"
JOB5_OUT=$(bsub -P acc_paul_oreilly -J c1_visualize \
  -w "done(${JOB4_ID})" \
  -q express -n 1 -R "rusage[mem=8000]" -W 0:30 \
  -o "${LOG_DIR}/visualize.out" -e "${LOG_DIR}/visualize.err" \
  "cd ${PROJECT_ROOT}; module load R/4.2.0; \
   Rscript scripts/c1_ancestry_test/05_visualize_results.R --config ${CONFIG}")
JOB5_ID=$(echo "${JOB5_OUT}" | grep -oP '\d+' | head -1)
echo "  Job ID: ${JOB5_ID}"

# =====================================================================
# Step 6 (optional): LD simulation — validate gene/pathway convergence
# Set RUN_SIMULATION=1 to enable (e.g., RUN_SIMULATION=1 bash submit_c1_pipeline.sh)
# =====================================================================

JOB6_ID=""
if [ "${RUN_SIMULATION:-0}" = "1" ]; then
  echo "--- Step 6: Submitting LD simulation job ---"
  JOB6_OUT=$(bsub -P acc_paul_oreilly -J c1_ld_simulation \
    -w "done(${JOB5_ID})" \
    -q express -n 1 -R "rusage[mem=16000]" -W 2:00 \
    -o "${LOG_DIR}/ld_simulation.out" -e "${LOG_DIR}/ld_simulation.err" \
    "cd ${PROJECT_ROOT}; module load R/4.2.0 plink/1.90b6.21; \
     Rscript scripts/c1_ancestry_test/06_simulate_ld_null.R --config ${CONFIG}")
  JOB6_ID=$(echo "${JOB6_OUT}" | grep -oP '\d+' | head -1)
  echo "  Job ID: ${JOB6_ID}"
else
  echo "--- Step 6: LD simulation SKIPPED (set RUN_SIMULATION=1 to enable) ---"
fi

# =====================================================================
# Summary
# =====================================================================

echo ""
echo "=== All jobs submitted ==="
echo ""
echo "  Dependency chain:"
echo "    Step 0: c1_subsample (${JOB0_ID})"
echo "      -> Step 1: c1_gwas[1-33] (${JOB1_ID})"
echo "           -> Step 2: c1_bnmf_eur[1-10] (${JOB2_ID})"
echo "           -> Step 3: c1_bnmf_afr (${JOB3_ID})"
echo "                -> Step 4: c1_jaccard (${JOB4_ID})"
echo "                     -> Step 5: c1_visualize (${JOB5_ID})"
if [ -n "${JOB6_ID}" ]; then
echo "                          -> Step 6: c1_ld_simulation (${JOB6_ID})"
fi
echo ""
echo "  Monitor with: bjobs -w | grep c1_"
echo "  Logs in: ${LOG_DIR}/"
echo "  Results in: ${RESULTS_DIR}/"
