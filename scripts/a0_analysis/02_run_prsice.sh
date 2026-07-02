#!/bin/bash
# 02_run_prsice.sh
# Run PRSice2 C+T for 6 traits using imputed UKB genotype data.
#
# PRSice parameters:
#   --clump-kb 250, --clump-r2 0.1 (standard C+T)
#   --bar-levels at multiple p-value thresholds
#   --fastscore --no-regress (raw scores only, regression in step 03)
#
# Usage:
#   # Run all traits:
#   bash scripts/a0_analysis/02_run_prsice.sh
#
#   # Run a single trait:
#   bash scripts/a0_analysis/02_run_prsice.sh T2D
#
#   # Submit as LSF array job:
#   bsub -J "prsice[1-6]" -n 4 -R "rusage[mem=8000]" -W 4:00 \
#     -q premium -P acc_paul_oreilly \
#     -o logs/prsice_%J_%I.stdout -e logs/prsice_%J_%I.stderr \
#     bash scripts/a0_analysis/02_run_prsice.sh LSF_ARRAY

set -euo pipefail

# --- Configuration ---
PROJECT_DIR="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
PRSICE_BIN="${PROJECT_DIR}/tools/PRSice_linux"
TARGET_PREFIX="/sc/arion/projects/paul_oreilly/data/Biobanks/UKB/genotyped/ukb18177"
SUMSTATS_DIR="${PROJECT_DIR}/results/a0_analysis/prs_ct/sumstats_prsice"
PHENO_FILE="${PROJECT_DIR}/results/a0_analysis/prs_ct/phenotypes_combined.txt"
COVAR_FILE="${PROJECT_DIR}/results/a0_analysis/prs_ct/covariates.txt"
OUT_DIR="${PROJECT_DIR}/results/a0_analysis/prs_ct/prsice_output"

# Optional sample restriction: PLINK keep file (FID IID). Restricts BOTH the
# clumping LD reference and scoring to this subset. Leave empty for full sample.
# Set to the 478,636 ancestry-assigned cohort to match downstream B1 analyses.
KEEP_FILE="${PROJECT_DIR}/results/a0_analysis/prs_ct/keep_478636.keep"

# PRSice C+T parameters
CLUMP_KB=250
CLUMP_R2=0.1
BAR_LEVELS="5e-8,1e-5,1e-3,0.01,0.05,0.1,0.5,1"
THREADS=4

# Trait array
TRAITS=("T2D" "CAD" "Angina" "MI" "Stroke" "PAD")

mkdir -p "${OUT_DIR}" "${PROJECT_DIR}/logs"

echo "=== PRSice C+T PRS Computation ==="
echo "  PRSice:      ${PRSICE_BIN}"
echo "  Target:      ${TARGET_PREFIX}"
echo "  Sumstats:    ${SUMSTATS_DIR}"
echo "  Pheno:       ${PHENO_FILE}"
echo "  Covar:       ${COVAR_FILE}"
echo "  Output:      ${OUT_DIR}"
echo "  Clump:       kb=${CLUMP_KB}, r2=${CLUMP_R2}"
echo "  Thresholds:  ${BAR_LEVELS}"
echo "  Threads:     ${THREADS}"
echo ""

# --- Determine which trait(s) to run ---
if [[ "${1:-}" == "LSF_ARRAY" ]]; then
    # LSF array job mode: use LSB_JOBINDEX
    IDX=$((LSB_JOBINDEX - 1))
    TRAIT_LIST=("${TRAITS[$IDX]}")
elif [[ -n "${1:-}" && "${1:-}" != "LSF_ARRAY" ]]; then
    # Single trait mode
    TRAIT_LIST=("$1")
else
    # Run all traits
    TRAIT_LIST=("${TRAITS[@]}")
fi

# --- Run PRSice for each trait ---
for TRAIT in "${TRAIT_LIST[@]}"; do
    echo "--- Running PRSice for ${TRAIT} ---"

    BASE_FILE="${SUMSTATS_DIR}/${TRAIT}_prsice.txt"

    if [[ ! -f "${BASE_FILE}" ]]; then
        echo "  ERROR: Base file not found: ${BASE_FILE}"
        continue
    fi

    # Use --extract with .valid file if it exists (dedup from prior run)
    EXTRACT_FLAG=""
    VALID_FILE="${OUT_DIR}/${TRAIT}.valid"
    if [[ -f "${VALID_FILE}" ]]; then
        echo "  Using --extract ${VALID_FILE}"
        EXTRACT_FLAG="--extract ${VALID_FILE}"
    fi

    # Restrict target (clumping + scoring) to a sample subset if KEEP_FILE is set
    KEEP_FLAG=""
    if [[ -n "${KEEP_FILE}" && -f "${KEEP_FILE}" ]]; then
        echo "  Using --keep ${KEEP_FILE}"
        KEEP_FLAG="--keep ${KEEP_FILE}"
    fi

    ${PRSICE_BIN} \
        --base "${BASE_FILE}" \
        --snp SNP --a1 A1 --a2 A2 --pvalue P --stat BETA \
        --target "${TARGET_PREFIX}" \
        --pheno "${PHENO_FILE}" \
        --pheno-col "${TRAIT}" \
        --cov "${COVAR_FILE}" \
        --cov-col age,age2,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Batch \
        --cov-factor Batch \
        --binary-target T \
        --beta \
        --clump-kb "${CLUMP_KB}" \
        --clump-r2 "${CLUMP_R2}" \
        --bar-levels "${BAR_LEVELS}" \
        --fastscore \
        --no-regress \
        --thread "${THREADS}" \
        ${EXTRACT_FLAG} \
        ${KEEP_FLAG} \
        --out "${OUT_DIR}/${TRAIT}" \
        2>&1 | tee "${PROJECT_DIR}/logs/prsice_${TRAIT}.log"

    echo "  Done: ${TRAIT}"
    echo ""
done

echo "=== PRSice complete ==="
