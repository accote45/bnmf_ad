#!/bin/bash
#BSUB -J harmonize_afib
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_afib_%J.out
#BSUB -e logs/lsf/harmonize_afib_%J.err

# Harmonize 12 Atrial Fibrillation GWAS files:
#   Yuan (3: META, EUR, EUR-noUKB)
#   UKB (2: AFR, SAS)
#   Verma/MVP (4: META, EUR, AFR, AMR)
#   Miyazawa/BBJ+UKB (1: META)
#   Sakaue/BBJ (1: EAS)
#   Roselli (1: META)
# All are GWAS Catalog harmonized format, all GRCh37.

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

# Files to harmonize: raw_filename|ancestry
FILES=(
    # Yuan 2025 (BETA, GWAS Catalog format, has n)
    "Yuan_NatureComms_2025.Afib.META.GRCh37.tsv.gz|META"
    "Yuan_NatureComms_2025.Afib.EUR.GRCh37.tsv.gz|EUR"
    "Yuan_NatureComms_2025.Afib.EUR.GRCh37.noukb.tsv.gz|EUR"
    # UKB 2025 (GWAS Catalog format)
    "UKB_Nature_2025.Afib.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.Afib.SAS.GRCh37.tsv.gz|SAS"
    # Verma/MVP 2024 (OR, GWAS Catalog format)
    "Verma_Science_2024.Afib.META.GRCh37.tsv.gz|META"
    "Verma_Science_2024.Afib.EUR.GRCh37.tsv.gz|EUR"
    "Verma_Science_2024.Afib.AFR.GRCh37.tsv.gz|AFR"
    "Verma_Science_2024.Afib.AMR.GRCh37.tsv.gz|AMR"
    # Miyazawa/BBJ+UKB 2023 (BETA, GWAS Catalog format)
    "Miyazawa_NatureGenetics_2023.Afib.META.GRCh37.tsv.gz|META"
    # Sakaue/BBJ 2021 (OR, GWAS Catalog format)
    "Sakaue_NatGenet_2021.Afib.EAS.GRCh37.tsv.gz|EAS"
    # Roselli 2018 (OR, GWAS Catalog format)
    "Roselli_NatGenet_2018.Afib.META.GRCh37.tsv.gz|META"
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

    # Harmonize (no reference GWAS — first Afib batch, no pre-existing harmonized files)
    echo "  Harmonizing..."
    if python "${SCRIPTS}/harmonize_sumstats.py" \
        --input-file "${RAW_FILE}" \
        --build-info "${BUILD_JSON}" \
        --output-file "${OUT_FILE}" \
        --qc-json "${QC_JSON}" \
        --preferred-build "${PREFERRED_BUILD}" \
        --n-cores "${N_CORES}"; then
        echo "  Done: ${OUT_FILE}"
    else
        echo "  FAILED: ${BASE}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Afib batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized Afib files: $(ls ${HARMONIZED}/*Afib*.processed.txt.gz 2>/dev/null | wc -l)"
