#!/bin/bash
#BSUB -J harm_batch1_afr
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 36:00
#BSUB -o logs/lsf/harm_batch1_afr_%J.out
#BSUB -e logs/lsf/harm_batch1_afr_%J.err

# Batch 1: AFR harmonization (17 files, ~19.5 GB)

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

mkdir -p "${HARMONIZED}" "${BUILD_DIR}" "${QC_DIR}"

declare -A T2D_REF=(
    ["AFR"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.AFR.GRCh37.txt.gz"
)

# Files to harmonize: raw_filename|ancestry
FILES=(
    "UKB_Nature_2025.ALP.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.Creatinine.AFR.GRCh37.tsv.gz|AFR"
    "Verma_Science_2024.Creatinine.AFR.GRCh37.txt.gz|AFR"
    "UKB_Nature_2025.CystatinC.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.Eosinophilcount.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.IGF1.AFR.GRCh37.tsv.gz|AFR"
    "Singh_NatureComms_2023.MAP.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.MeanReticVol.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.Phosphate.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.PltDistWidth.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.Protein.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.RandomGlucose.AFR.GRCh37.tsv.gz|AFR"
    "Verma_Science_2024.RandomGlucose.AFR.GRCh37.txt.gz|AFR"
    "Karczewski_NatGenet_2025.ReticulocyteCount.AFR.GRCh37.txt.gz|AFR"
    "UKB_Nature_2025.SHBG.AFR.GRCh37.tsv.gz|AFR"
    "Pagadala_NatureComms_2025.Testosterone.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.Urea.AFR.GRCh37.tsv.gz|AFR"
)

TOTAL=${#FILES[@]}
COUNT=0

for ENTRY in "${FILES[@]}"; do
    IFS='|' read -r RAW_FNAME ANCESTRY <<< "${ENTRY}"
    COUNT=$((COUNT + 1))

    RAW_FILE="${SUMSTATS}/${RAW_FNAME}"
    BASE=$(echo "${RAW_FNAME}" | sed 's/\(GRCh[0-9]*\).*/\1/')
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
    echo "[${COUNT}/${TOTAL}] Processing: ${BASE}"

    # Build check
    if [ ! -f "${BUILD_JSON}" ]; then
        echo "  Running build check..."
        python "${SCRIPTS}/check_build.py" \
            --input-file "${RAW_FILE}" \
            --output-file "${BUILD_JSON}"
    fi

    # Fix Unknown builds
    if grep -q '"actual_build": "Unknown"' "${BUILD_JSON}" 2>/dev/null; then
        echo "  Build detection failed — overriding to GRCh37 (labeled build)"
        python -c "
import json
with open('${BUILD_JSON}') as f:
    d = json.load(f)
d['actual_build'] = '19'
d['override_reason'] = 'Build detection failed; using labeled GRCh37'
with open('${BUILD_JSON}', 'w') as f:
    json.dump(d, f, indent=2)
"
    fi

    # Harmonize with ancestry-matched T2D reference
    REF_GWAS="${T2D_REF[${ANCESTRY}]}"
    echo "  Harmonizing (ref: ${ANCESTRY} T2D)..."
    python "${SCRIPTS}/harmonize_sumstats.py" \
        --input-file "${RAW_FILE}" \
        --build-info "${BUILD_JSON}" \
        --output-file "${OUT_FILE}" \
        --qc-json "${QC_JSON}" \
        --preferred-build "${PREFERRED_BUILD}" \
        --n-cores "${N_CORES}" \
        --ref-gwas "${REF_GWAS}"

    echo "  Done: ${OUT_FILE}"
done

echo ""
echo "=== Batch 1 (AFR) harmonization complete ==="
