#!/bin/bash
#BSUB -J harm_batch3_eas_sas
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=32000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harm_batch3_eas_sas_%J.out
#BSUB -e logs/lsf/harm_batch3_eas_sas_%J.err

# Batch 3: EAS + SAS harmonization (11 files, ~4.2 GB)

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
    ["EAS"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.EAS.GRCh37.processed.txt.gz"
    ["SAS"]="${HARMONIZED}/Mahajan_NatureGenetics_2022.t2d.SAS.GRCh37.processed.txt.gz"
)

# Files to harmonize: raw_filename|ancestry
FILES=(
    "Verma_Science_2024.Creatinine.EAS.GRCh37.txt.gz|EAS"
    "Sakaue_NatGenet_2021.Eosinophilcount.EAS.GRCh37.txt.gz|EAS"
    "Sakaue_NatGenet_2021.RandomGlucose.EAS.GRCh37.txt.gz|EAS"
    "Verma_Science_2024.RandomGlucose.EAS.GRCh37.txt.gz|EAS"
    "Karczewski_NatureGenetics_2025.ALT.SAS.GRCh37.txt.gz|SAS"
    "Karczewski_NatureGenetics_2025.ApoA.SAS.GRCh37.txt.gz|SAS"
    "Graham_Nature_2021.HDL.SAS.GRCh37.txt.gz|SAS"
    "Chen_Cell_2020.MonocyteCount.SAS.GRCh37.txt.gz|SAS"
    "UKB_Nature_2025.RandomGlucose.SAS.GRCh37.tsv.gz|SAS"
    "Graham_Nature_2021.TotalCholesterol.SAS.GRCh37.txt.gz|SAS"
    "Graham_Nature_2021.Triglycerides.SAS.GRCh37.txt.gz|SAS"
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
echo "=== Batch 3 (EAS + SAS) harmonization complete ==="
