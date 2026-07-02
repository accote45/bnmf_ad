#!/bin/bash
#BSUB -J harmonize_hgb
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 24:00
#BSUB -o logs/lsf/harmonize_hgb_%J.out
#BSUB -e logs/lsf/harmonize_hgb_%J.err

# Harmonize 15 HGB GWAS files:
#   Chen Cell 2020 (5: AFR, AMR, EAS, EUR, SAS) — beta + OR, no N
#   Chen CellGenomics 2023 (1: EAS) — beta, has N
#   Gurdasani Cell 2019 (1: AFR) — beta_fe, space-sep, compound snpid, no N
#   Jacobs NatureComms 2024 (1: SAS) — beta, no N
#   Sakaue NatGenet 2021 (1: EAS) — beta + OR, no N
#   Verma Science 2024 (5: META, EUR, AFR, EAS, AMR) — beta, has N
#   Vuckovic Cell 2020 (1: EUR) — beta + OR, no N
# No reference GWAS — first HGB batch.

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
    # Chen Cell 2020 (beta + OR, no N)
    "Chen_Cell_2020.Hgb.AFR.GRCh37.txt.gz|AFR"
    "Chen_Cell_2020.Hgb.AMR.GRCh37.txt.gz|AMR"
    "Chen_Cell_2020.Hgb.EAS.GRCh37.txt.gz|EAS"
    "Chen_Cell_2020.Hgb.EUR.GRCh37.txt.gz|EUR"
    "Chen_Cell_2020.Hgb.SAS.GRCh37.txt.gz|SAS"
    # Chen CellGenomics 2023 (beta, has N)
    "Chen_CellGenomics_2023.Hgb.EAS.GRCh37.tsv.gz|EAS"
    # Gurdasani Cell 2019 (beta_fe, space-sep, compound snpid, no N)
    "Gurdasani_Cell_2019.Hgb.AFR.GRCh37.txt.gz|AFR"
    # Jacobs NatureComms 2024 (beta, no N)
    "Jacobs_NatureComms_2024.Hgb.SAS.GRCh37.txt.gz|SAS"
    # Sakaue NatGenet 2021 (beta + OR, no N)
    "Sakaue_NatGenet_2021.Hgb.EAS.GRCh37.txt.gz|EAS"
    # Verma Science 2024 (beta, has N)
    "Verma_Science_2024.Hgb.META.GRCh37.txt.gz|META"
    "Verma_Science_2024.Hgb.EUR.GRCh37.txt.gz|EUR"
    "Verma_Science_2024.Hgb.AFR.GRCh37.txt.gz|AFR"
    "Verma_Science_2024.Hgb.EAS.GRCh37.txt.gz|EAS"
    "Verma_Science_2024.Hgb.AMR.GRCh37.txt.gz|AMR"
    # Vuckovic Cell 2020 (beta + OR, no N)
    "Vuckovic_Cell_2020.Hgb.EUR.GRCh37.txt.gz|EUR"
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

    # Harmonize (no reference GWAS — first HGB batch)
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
echo "=== HGB batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized HGB files: $(ls ${HARMONIZED}/*Hgb*.processed.txt.gz 2>/dev/null | wc -l)"
