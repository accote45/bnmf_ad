#!/bin/bash
#BSUB -J harmonize_basophil
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_basophil_%J.out
#BSUB -e logs/lsf/harmonize_basophil_%J.err

# Harmonize 7 Basophil Count GWAS files:
#   Sakaue/BBJ 2021 (2: META, EAS) — OR + beta, GWAS Catalog format
#   UKB 2025 (3: EUR, AFR, SAS) — GWAS Catalog format, has n
#   Vuckovic 2020 (1: EUR) — GWAS Catalog format, detects as build 38 → liftover
#   Gurdasani 2019 (1: AFR) — Raw space-delimited, compound snpid

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
    # Sakaue/BBJ 2021 (OR + BETA, GWAS Catalog format)
    "Sakaue_NatGenet_2021.Basophilcount.META.GRCh37.txt.gz|META"
    "Sakaue_NatGenet_2021.Basophilcount.EAS.GRCh37.txt.gz|EAS"
    # UKB 2025 (GWAS Catalog format, has n)
    "UKB_Nature_2025.Basophilcount.EUR.GRCh37.tsv.gz|EUR"
    "UKB_Nature_2025.Basophilcount.AFR.GRCh37.tsv.gz|AFR"
    "UKB_Nature_2025.Basophilcount.SAS.GRCh37.tsv.gz|SAS"
    # Vuckovic 2020 (GWAS Catalog format, build 38 detected → liftover needed)
    "Vuckovic_Cell_2020.Basophilcount.EUR.GRCh37.txt.gz|EUR"
    # Gurdasani 2019 (Raw, compound snpid, beta_fe/se_fe/pval_fe)
    "Gurdasani_Cell_2019.Basophilcount.AFR.GRCh37.txt.gz|AFR"
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

    # Harmonize (no reference GWAS — first basophil batch, no pre-existing harmonized files)
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
echo "=== Basophil batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized Basophil files: $(ls ${HARMONIZED}/*Basophilcount*.processed.txt.gz 2>/dev/null | wc -l)"
