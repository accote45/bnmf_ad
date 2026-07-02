#!/bin/bash
#BSUB -J harmonize_plt
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 8:00
#BSUB -o logs/lsf/harmonize_plt_%J.out
#BSUB -e logs/lsf/harmonize_plt_%J.err

# Harmonize all Platelet Count (PltCount) GWAS files (16 total):
#   Karczewski_NatGenet_2025 — EUR, META, SAS, AFR
#   Verma_Science_2024       — EUR, EAS, AMR, META
#   Sakaue_NatGenet_2021     — EAS
#   Chen_CellGenomics_2023   — EAS
#   Chen_Cell_2020           — EUR, AFR, EAS, SAS, AMR
#   Gurdasani_Cell_2019      — AFR

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

# T2D references for allele alignment (ancestry-matched)
declare -A T2D_REF=(
    ["EUR"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.EUR.GRCh37.processed.txt.gz"
    ["META"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"
    ["EAS"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.EAS.GRCh37.processed.txt.gz"
    ["SAS"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"
    ["AMR"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"
    ["AFR"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"
)

# Files to harmonize: raw_filename|ancestry
FILES=(
    # Karczewski NatGenet 2025
    "Karczewski_NatGenet_2025.PltCount.EUR.GRCh37.txt.gz|EUR"
    "Karczewski_NatGenet_2025.PltCount.META.GRCh37.txt.gz|META"
    "Karczewski_NatGenet_2025.PltCount.SAS.GRCh37.txt.gz|SAS"
    "Karczewski_NatGenet_2025.PltCount.AFR.GRCh37.txt.gz|AFR"
    # Verma Science 2024
    "Verma_Science_2024.PltCount.EUR.GRCh37.txt.gz|EUR"
    "Verma_Science_2024.PltCount.EAS.GRCh37.txt.gz|EAS"
    "Verma_Science_2024.PltCount.AMR.GRCh37.txt.gz|AMR"
    "Verma_Science_2024.PltCount.META.GRCh37.txt.gz|META"
    # Sakaue NatGenet 2021
    "Sakaue_NatGenet_2021.PltCount.EAS.GRCh37.txt.gz|EAS"
    # Chen CellGenomics 2023
    "Chen_CellGenomics_2023.PltCount.EAS.GRCh37.tsv.gz|EAS"
    # Chen Cell 2020
    "Chen_Cell_2020.PltCount.EUR.GRCh37.txt.gz|EUR"
    "Chen_Cell_2020.PltCount.AFR.GRCh37.txt.gz|AFR"
    "Chen_Cell_2020.PltCount.EAS.GRCh37.txt.gz|EAS"
    "Chen_Cell_2020.PltCount.SAS.GRCh37.txt.gz|SAS"
    "Chen_Cell_2020.PltCount.AMR.GRCh37.txt.gz|AMR"
    # Gurdasani Cell 2019 (custom METAL format)
    "Gurdasani_Cell_2019.PltCount.AFR.GRCh37.txt.gz|AFR"
)

TOTAL=${#FILES[@]}
COUNT=0
FAILED=0

for ENTRY in "${FILES[@]}"; do
    IFS='|' read -r RAW_FNAME ANCESTRY <<< "${ENTRY}"
    COUNT=$((COUNT + 1))

    RAW_FILE="${SUMSTATS}/${RAW_FNAME}"

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

    HARM_CMD="python ${SCRIPTS}/harmonize_sumstats.py \
        --input-file ${RAW_FILE} \
        --build-info ${BUILD_JSON} \
        --output-file ${OUT_FILE} \
        --qc-json ${QC_JSON} \
        --preferred-build ${PREFERRED_BUILD} \
        --n-cores ${N_CORES}"

    if [[ -v "T2D_REF[${ANCESTRY}]" ]] && [ -f "${T2D_REF[${ANCESTRY}]}" ]; then
        HARM_CMD="${HARM_CMD} --ref-gwas ${T2D_REF[${ANCESTRY}]}"
        echo "  Harmonizing (ref: ${ANCESTRY} T2D)..."
    else
        echo "  Harmonizing (no ${ANCESTRY} T2D reference available)..."
    fi

    if eval "${HARM_CMD}"; then
        echo "  Done: ${OUT_FILE}"
    else
        echo "  FAILED: ${BASE}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== PltCount batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized files: $(ls ${HARMONIZED}/*PltCount*.processed.txt.gz 2>/dev/null | wc -l)"
