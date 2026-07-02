#!/bin/bash
# 01_run_gwas.sh
# Run PLINK2 GWAS for a single ancestry/subsample/trait combination.
# Supports both single merged genotype files and per-chromosome imputed data.
#
# Interactive usage:
#   bash scripts/c1_ancestry_test/01_run_gwas.sh --ancestry EUR --subsample 1 --trait LDL
#   bash scripts/c1_ancestry_test/01_run_gwas.sh --ancestry AFR --trait BMI
#
# LSF array usage (legacy):
#   Uses $LSB_JOBINDEX when no --ancestry/--trait flags are provided.

set -euo pipefail

module load plink2/2.3 2>/dev/null || module load plink2 2>/dev/null || true

# --- Configuration ---
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
RESULTS_DIR="${PROJECT_ROOT}/results/c1_ancestry_test"
PHENO_FILE="${PHENO_FILE:-${RESULTS_DIR}/phenotypes.txt}"
COVAR_FILE="${COVAR_FILE:-${RESULTS_DIR}/covariates.txt}"
COVAR_NAMES="age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"

# Genotype data: per-chromosome imputed (default) or single merged
GENO_PREFIX="${GENO_PREFIX:-/sc/arion/projects/paul_oreilly/data/Biobanks/UKB/imputed/qced_data/chr}"
GENO_TYPE="${GENO_TYPE:-per_chromosome}"  # "per_chromosome" or "single"

TRAITS=("LDL" "glucose" "BMI")
N_SUBS=${N_SUBS:-5}

# --- Parse arguments ---
ANCESTRY=""
SUB_IDX=""
TRAIT=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --ancestry) ANCESTRY="$2"; shift 2 ;;
    --subsample) SUB_IDX="$2"; shift 2 ;;
    --trait) TRAIT="$2"; shift 2 ;;
    --geno-prefix) GENO_PREFIX="$2"; shift 2 ;;
    --geno-type) GENO_TYPE="$2"; shift 2 ;;
    *) shift ;;
  esac
done

# Fall back to LSF array index if no args provided
if [ -z "${ANCESTRY}" ] && [ -n "${LSB_JOBINDEX:-}" ]; then
  N_EUR_JOBS=$(( N_SUBS * 3 ))
  IDX=${LSB_JOBINDEX}
  if [ ${IDX} -le ${N_EUR_JOBS} ]; then
    ANCESTRY="EUR"
    TRAIT_IDX=$(( (${IDX} - 1) / ${N_SUBS} ))
    SUB_IDX=$(( (${IDX} - 1) % ${N_SUBS} + 1 ))
    TRAIT=${TRAITS[${TRAIT_IDX}]}
  else
    ANCESTRY="AFR"
    TRAIT_IDX=$(( ${IDX} - ${N_EUR_JOBS} - 1 ))
    TRAIT=${TRAITS[${TRAIT_IDX}]}
  fi
fi

if [ -z "${ANCESTRY}" ] || [ -z "${TRAIT}" ]; then
  echo "Usage: bash 01_run_gwas.sh --ancestry EUR|AFR --trait LDL|glucose|BMI [--subsample N]"
  exit 1
fi

# --- Determine output directory and keep file ---
if [ "${ANCESTRY}" = "EUR" ]; then
  if [ -z "${SUB_IDX}" ]; then
    echo "ERROR: --subsample required for EUR ancestry"
    exit 1
  fi
  SUB_LABEL=$(printf "sub_%02d" ${SUB_IDX})
  KEEP_FILE="${RESULTS_DIR}/subsamples/${SUB_LABEL}/subsample_$(printf '%02d' ${SUB_IDX}).keep"
  OUT_DIR="${RESULTS_DIR}/subsamples/${SUB_LABEL}"
else
  KEEP_FILE="${RESULTS_DIR}/afr/afr_all.keep"
  OUT_DIR="${RESULTS_DIR}/afr"
fi

mkdir -p "${OUT_DIR}"

echo "=== C1.1 GWAS ==="
echo "Ancestry: ${ANCESTRY}"
echo "Trait: ${TRAIT}"
echo "Keep file: ${KEEP_FILE}"
echo "Geno type: ${GENO_TYPE}"
echo "Output dir: ${OUT_DIR}"

# --- Validate inputs ---
if [ ! -f "${KEEP_FILE}" ]; then
  echo "ERROR: Keep file not found: ${KEEP_FILE}"
  exit 1
fi

if [ ! -f "${PHENO_FILE}" ]; then
  echo "ERROR: Phenotype file not found: ${PHENO_FILE}"
  exit 1
fi

N_SAMPLES=$(wc -l < "${KEEP_FILE}")
echo "N samples in keep file: ${N_SAMPLES}"

# --- Run PLINK2 GWAS ---
OUT_PREFIX="${OUT_DIR}/${TRAIT}_gwas"
COMBINED_GLM="${OUT_PREFIX}.${TRAIT}.glm.linear"

if [ "${GENO_TYPE}" = "per_chromosome" ]; then
  # Per-chromosome imputed data: run per chr, concatenate
  echo "Running per-chromosome GWAS (chr1-22)..."
  HEADER_WRITTEN=false

  for CHR in $(seq 1 22); do
    BFILE="${GENO_PREFIX}${CHR}"
    if [ ! -f "${BFILE}.bed" ]; then
      echo "  WARNING: Skipping chr${CHR} (${BFILE}.bed not found)"
      continue
    fi

    CHR_PREFIX="${OUT_DIR}/${TRAIT}_gwas_chr${CHR}"
    echo -n "  chr${CHR}..."

    plink2 \
      --bfile "${BFILE}" \
      --keep "${KEEP_FILE}" \
      --pheno "${PHENO_FILE}" \
      --pheno-name "${TRAIT}" \
      --covar "${COVAR_FILE}" \
      --covar-name ${COVAR_NAMES} \
      --covar-variance-standardize \
      --freq \
      --glm hide-covar cols=+a1freq \
      --out "${CHR_PREFIX}" \
      --threads 4 \
      2>/dev/null

    CHR_GLM="${CHR_PREFIX}.${TRAIT}.glm.linear"
    if [ -f "${CHR_GLM}" ]; then
      if [ "${HEADER_WRITTEN}" = false ]; then
        cp "${CHR_GLM}" "${COMBINED_GLM}"
        HEADER_WRITTEN=true
      else
        tail -n +2 "${CHR_GLM}" >> "${COMBINED_GLM}"
      fi
      # Clean up per-chr files
      rm -f "${CHR_PREFIX}".*
      echo " done"
    else
      echo " no output (possibly skipped by plink2)"
    fi
  done

else
  # Single merged file
  echo "Running PLINK2 GWAS (single file)..."
  plink2 \
    --bfile "${GENO_PREFIX}" \
    --keep "${KEEP_FILE}" \
    --pheno "${PHENO_FILE}" \
    --pheno-name "${TRAIT}" \
    --covar "${COVAR_FILE}" \
    --covar-name ${COVAR_NAMES} \
    --covar-variance-standardize \
    --freq \
    --glm hide-covar cols=+a1freq \
    --out "${OUT_PREFIX}" \
    --threads 4
fi

# --- Check output ---
if [ ! -f "${COMBINED_GLM}" ]; then
  COMBINED_GLM=$(ls ${OUT_PREFIX}*.glm.linear 2>/dev/null | head -1)
fi

if [ -z "${COMBINED_GLM}" ] || [ ! -f "${COMBINED_GLM}" ]; then
  echo "ERROR: GWAS output not found. Expected: ${OUT_PREFIX}.*.glm.linear"
  exit 1
fi

N_VARIANTS=$(wc -l < "${COMBINED_GLM}")
echo "GWAS complete: ${N_VARIANTS} lines"

# --- Format output ---
echo "Formatting GWAS output..."
Rscript "${PROJECT_ROOT}/scripts/c1_ancestry_test/01_format_gwas.R" \
  --input "${COMBINED_GLM}" \
  --output "${OUT_DIR}/${TRAIT}_formatted.txt.gz" \
  --trait "${TRAIT}"

# Clean up combined raw file
rm -f "${COMBINED_GLM}"

echo "=== Done: ${ANCESTRY} ${TRAIT} ==="
