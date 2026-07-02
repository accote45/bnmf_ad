#!/bin/bash
#BSUB -J harmonize_cad
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_cad_%J.out
#BSUB -e logs/lsf/harmonize_cad_%J.err

# Harmonize 9 CAD GWAS files: Verma (5), Karczewski (1), Nikpay (1), Aragam (2).
# Verma AFR and META are GRCh38 and need liftover to GRCh19.

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

# CAD references for allele alignment (ancestry-matched, already harmonized)
declare -A CAD_REF=(
    ["EUR"]="${HARMONIZED}/Tcheandjieu_NatureMed_2023.CAD.EUR.GRCh37.processed.txt.gz"
    ["AFR"]="${HARMONIZED}/Tcheandjieu_NatureMed_2023.CAD.AFR.GRCh37.processed.txt.gz"
    ["META"]="${HARMONIZED}/Tcheandjieu_NatureMed_2023.CAD.META.GRCh37.processed.txt.gz"
    ["EAS"]="${HARMONIZED}/Sakaue_NatGenet_2021.CAD.EAS.GRCh37.processed.txt.gz"
)
# No SAS or AMR CAD reference available

# Files to harmonize: raw_filename|ancestry
FILES=(
    # Verma 2024 (OR only, GWAS Catalog format, has n)
    "Verma_Science_2024.CAD.EUR.GRCh37.txt.gz|EUR"
    "Verma_Science_2024.CAD.AFR.GRCh38.txt.gz|AFR"
    "Verma_Science_2024.CAD.AMR.GRCh37.txt.gz|AMR"
    "Verma_Science_2024.CAD.EAS.GRCh37.txt.gz|EAS"
    "Verma_Science_2024.CAD.META.GRCh38.txt.gz|META"
    # Karczewski 2025 (BETA, GWAS Catalog format)
    "Karczewski_NatGenet_2025.CAD.SAS.GRCh37.txt.gz|SAS"
    # Nikpay 2015 (has both beta and OR)
    "Nikpay_NatureGenetics_2015.CAD.EUR.GRCh37.txt.gz|EUR"
    # Aragam 2022 (has both beta and OR, has n)
    "Aragam_NatureGenetics_2022.CAD.EUR.GRCh37.tsv.gz|EUR"
    "Aragam_NatureGenetics_2022.CAD.META.GRCh37.tsv.gz|META"
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

    # Build harmonize command
    HARM_CMD="python ${SCRIPTS}/harmonize_sumstats.py \
        --input-file ${RAW_FILE} \
        --build-info ${BUILD_JSON} \
        --output-file ${OUT_FILE} \
        --qc-json ${QC_JSON} \
        --preferred-build ${PREFERRED_BUILD} \
        --n-cores ${N_CORES}"

    # Add ancestry-matched CAD reference if available
    if [[ -v "CAD_REF[${ANCESTRY}]" ]] && [ -f "${CAD_REF[${ANCESTRY}]}" ]; then
        HARM_CMD="${HARM_CMD} --ref-gwas ${CAD_REF[${ANCESTRY}]}"
        echo "  Harmonizing (ref: ${ANCESTRY} CAD)..."
    else
        echo "  Harmonizing (no ${ANCESTRY} CAD reference available)..."
    fi

    if eval "${HARM_CMD}"; then
        echo "  Done: ${OUT_FILE}"
    else
        echo "  FAILED: ${BASE}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== CAD batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized CAD files: $(ls ${HARMONIZED}/*CAD*.processed.txt.gz 2>/dev/null | wc -l)"
