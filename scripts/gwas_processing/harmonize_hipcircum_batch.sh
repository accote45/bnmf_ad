#!/bin/bash
#BSUB -J harmonize_hipcircum
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 24:00
#BSUB -o logs/lsf/harmonize_hipcircum_%J.out
#BSUB -e logs/lsf/harmonize_hipcircum_%J.err

# Harmonize 5 Hip Circumference GWAS files from 3 studies:
#   Gurdasani Cell 2019 (1: AFR) — space-sep, compound snpid, beta_fe/se_fe/pval_fe, no N
#   Chen CellGenomics 2023 (1: EAS) — tab-sep, beta, has N
#   Karczewski NatGenet 2025 (3: EUR, SAS, META) — tab-sep, beta, no N
# No reference GWAS — first Hipcircum batch.
# Note: Gurdasani file may need 128GB if it fails on 64GB (see harmonize_rerun_gurdasani_height.sh).

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
    # Gurdasani Cell 2019 (space-sep, compound snpid, beta_fe/se_fe/pval_fe, no N)
    "Gurdasani_Cell_2019.Hipcircum.AFR.GRCh37.txt.gz|AFR"
    # Chen CellGenomics 2023 (tab-sep, beta, has N)
    "Chen_CellGenomics_2023.Hipcircum.EAS.GRCh37.tsv.gz|EAS"
    # Karczewski NatGenet 2025 (tab-sep, beta, no N)
    "Karczewski_NatGenet_2025.Hipcircum.EUR.GRCh37.txt.gz|EUR"
    "Karczewski_NatGenet_2025.Hipcircum.SAS.GRCh37.txt.gz|SAS"
    "Karczewski_NatGenet_2025.Hipcircum.META.GRCh37.txt.gz|META"
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

    # Harmonize (no reference GWAS — first Hipcircum batch)
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
echo "=== Hipcircum batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized Hipcircum files: $(ls ${HARMONIZED}/*Hipcircum*.processed.txt.gz 2>/dev/null | wc -l)"
