#!/bin/bash
#BSUB -J harm_batch5_ldl_sas
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=128000]"
#BSUB -W 24:00
#BSUB -o logs/lsf/harm_batch5_ldl_sas_%J.out
#BSUB -e logs/lsf/harm_batch5_ldl_sas_%J.err

# Batch 5: Graham LDL SAS (1 file)

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
    ["SAS"]="${HARMONIZED}/Mahajan_NatureGenetics_2022.t2d.SAS.GRCh37.processed.txt.gz"
)

FILES=(
    "Graham_Nature_2021.LDL.SAS.GRCh37.txt.gz|SAS"
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

    if [ ! -f "${BUILD_JSON}" ]; then
        echo "  Running build check..."
        python "${SCRIPTS}/check_build.py" \
            --input-file "${RAW_FILE}" \
            --output-file "${BUILD_JSON}"
    fi

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
echo "=== Batch 5 (Graham LDL SAS) harmonization complete ==="
