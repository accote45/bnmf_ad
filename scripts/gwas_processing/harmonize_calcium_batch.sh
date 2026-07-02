#!/bin/bash
#BSUB -J harmonize_calcium
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_calcium_%J.out
#BSUB -e logs/lsf/harmonize_calcium_%J.err

# Harmonize 8 Calcium GWAS files: UKB (3), Verma/MVP (4), Sakaue/BBJ (1).
# All files are GRCh37 — no liftover needed.
# Sakaue file has no N column and needs --n-override 83980.
# No reference GWAS for allele alignment (first Calcium batch).

set -euo pipefail

PROJECT_ROOT="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
cd "${PROJECT_ROOT}"

SCRIPTS="${PROJECT_ROOT}/scripts/gwas_processing"
SUMSTATS="${PROJECT_ROOT}/sumstats"
HARMONIZED="${SUMSTATS}/harmonized"
BUILD_DIR="${SUMSTATS}/build_checks"
QC_DIR="${SUMSTATS}/qc"
PREFERRED_BUILD="19"
N_CORES=32

export OPENBLAS_NUM_THREADS=1

mkdir -p "${HARMONIZED}" "${BUILD_DIR}" "${QC_DIR}" logs/lsf

# N overrides for files missing sample size column
declare -A N_OVERRIDES=(
    ["Sakaue_NatureGenetics_2021.Calcium.EAS"]=83980
)

# Files to harmonize: raw_filename|ancestry
FILES=(
    # UKB Nature 2025 (GWAS Catalog format, has beta + n)
    "UKB_Nature_2025.Calcium.EUR.GRCh37.tsv.gz|EUR"
    "UKB_Nature_2025.Calcium.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.Calcium.SAS.GRCh37.tsv.gz|SAS"
    # Verma/MVP Science 2024 (GWAS Catalog meta-analysis format, has beta + n)
    "Verma_Science_2024.Calcium.META.GRCh37.tsv.gz|META"
    "Verma_Science_2024.Calcium.EUR.GRCh37.tsv.gz|EUR"
    "Verma_Science_2024.Calcium.AFR.GRCh37.tsv.gz|AFR"
    "Verma_Science_2024.Calcium.AMR.GRCh37.tsv.gz|AMR"
    # Sakaue/BBJ 2021 (has OR+beta, no N column — needs --n-override 83980)
    "Sakaue_NatureGenetics_2021.Calcium.EAS.GRCh37.txt.gz|EAS"
)

TOTAL=${#FILES[@]}
COUNT=0
FAILED=0

for ENTRY in "${FILES[@]}"; do
    IFS='|' read -r RAW_FNAME ANCESTRY <<< "${ENTRY}"
    COUNT=$((COUNT + 1))

    RAW_FILE="${SUMSTATS}/${RAW_FNAME}"

    # Derive base name: strip .txt.gz, .tsv.gz, or .txt extension
    BASE=$(echo "${RAW_FNAME}" | sed -E 's/\.(txt\.gz|tsv\.gz|txt|tsv)$//')
    BUILD_JSON="${BUILD_DIR}/${BASE}.build.json"
    OUT_FILE="${HARMONIZED}/${BASE}.processed.txt.gz"
    QC_JSON="${QC_DIR}/${BASE}.qc.json"

    if [ -f "${OUT_FILE}" ]; then
        echo "[${COUNT}/${TOTAL}] SKIP (exists): ${BASE}"
        continue
    fi

    if [ ! -f "${RAW_FILE}" ]; then
        echo "[${COUNT}/${TOTAL}] SKIP (not found): ${RAW_FNAME}"
        continue
    fi

    echo ""
    echo "[${COUNT}/${TOTAL}] Processing: ${BASE} (ancestry=${ANCESTRY})"

    # Build check
    if [ ! -f "${BUILD_JSON}" ]; then
        echo "  Running build check..."
        python "${SCRIPTS}/check_build.py" \
            --input-file "${RAW_FILE}" \
            --output-file "${BUILD_JSON}"
    fi

    # Fix Unknown builds — override to labeled build
    if grep -q '"actual_build": "Unknown"' "${BUILD_JSON}" 2>/dev/null; then
        # Extract labeled build number from filename (GRCh37 -> 19, GRCh38 -> 38)
        LABELED_BUILD=$(echo "${BASE}" | grep -oP 'GRCh\K[0-9]+')
        if [ "${LABELED_BUILD}" = "37" ]; then LABELED_BUILD="19"; fi
        echo "  Build detection failed — overriding to GRCh${LABELED_BUILD} (labeled build)"
        python -c "
import json
with open('${BUILD_JSON}') as f:
    d = json.load(f)
d['actual_build'] = '${LABELED_BUILD}'
d['override_reason'] = 'Build detection failed; using labeled build from filename'
with open('${BUILD_JSON}', 'w') as f:
    json.dump(d, f, indent=2)
"
    fi

    # Determine N override flag
    N_FLAG=""
    for key in "${!N_OVERRIDES[@]}"; do
        if [[ "${BASE}" == *"${key}"* ]]; then
            N_FLAG="--n-override ${N_OVERRIDES[${key}]}"
            echo "    Using N override: ${N_OVERRIDES[${key}]}"
            break
        fi
    done

    # Build harmonize command (no --ref-gwas: first Calcium batch, no reference available)
    HARM_CMD="python ${SCRIPTS}/harmonize_sumstats.py \
        --input-file ${RAW_FILE} \
        --build-info ${BUILD_JSON} \
        --output-file ${OUT_FILE} \
        --qc-json ${QC_JSON} \
        --preferred-build ${PREFERRED_BUILD} \
        --n-cores ${N_CORES} ${N_FLAG}"

    echo "  Harmonizing..."

    if eval "${HARM_CMD}"; then
        echo "  Done: ${OUT_FILE}"
    else
        echo "  FAILED: ${BASE}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Calcium batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized Calcium files: $(ls ${HARMONIZED}/*Calcium*.processed.txt.gz 2>/dev/null | wc -l)"
