#!/bin/bash
#BSUB -J harmonize_hypertension
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 24:00
#BSUB -o logs/lsf/harmonize_hypertension_%J.out
#BSUB -e logs/lsf/harmonize_hypertension_%J.err

# Harmonize 9 Hypertension GWAS files from 3 studies:
#   Karczewski NatGenet 2025 (3: EUR, SAS, META) — tab-sep, beta, no N
#   Singh NatureComms 2023 (1: AFR) — tab-sep, beta, has N
#   Verma Science 2024 (5: META, EUR, AFR, EAS, AMR) — tab-sep, odds_ratio,
#       SE all NA (derived from ci_upper/ci_lower), has N + case/control counts
# No reference GWAS — first Hypertension batch.

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
    # Karczewski NatGenet 2025 (tab-sep, beta, no N)
    "Karczewski_NatGenet_2025.Hypertension.EUR.GRCh37.txt.gz|EUR"
    "Karczewski_NatGenet_2025.Hypertension.SAS.GRCh37.txt.gz|SAS"
    "Karczewski_NatGenet_2025.Hypertension.META.GRCh37.txt.gz|META"
    # Singh NatureComms 2023 (tab-sep, beta, has N)
    "Singh_NatureComms_2023.Hypertension.AFR.GRCh37.txt.gz|AFR"
    # Verma Science 2024 (tab-sep, odds_ratio, SE from CI bounds, has N)
    "Verma_Science_2024.Hypertension.META.GRCh37.txt.gz|META"
    "Verma_Science_2024.Hypertension.EUR.GRCh37.txt.gz|EUR"
    "Verma_Science_2024.Hypertension.AFR.GRCh37.txt.gz|AFR"
    "Verma_Science_2024.Hypertension.EAS.GRCh37.txt.gz|EAS"
    "Verma_Science_2024.Hypertension.AMR.GRCh37.txt.gz|AMR"
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

    # Harmonize (no reference GWAS — first Hypertension batch)
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
echo "=== Hypertension batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized Hypertension files: $(ls ${HARMONIZED}/*Hypertension*.processed.txt.gz 2>/dev/null | wc -l)"
