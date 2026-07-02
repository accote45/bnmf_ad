#!/bin/bash
#BSUB -J harmonize_height
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 24:00
#BSUB -o logs/lsf/harmonize_height_%J.out
#BSUB -e logs/lsf/harmonize_height_%J.err

# Harmonize 15 Height GWAS files:
#   Karczewski NatGenet 2025 (4: META, EUR, AFR, SAS) — beta, no N
#   Verma Science 2024 (5: META, EUR, AFR, EAS, AMR) — beta, has N
#   Chen CellGenomics 2023 (1: EAS) — beta, has N
#   Sohail Nature 2023 (1: AMR) — beta, has N
#   FernandezRhodes HGGAdv 2022 (1: AMR) — beta + OR, no N
#   Sakaue NatGenet 2021 (1: EAS) — beta + OR, no N
#   Graff AJHG 2023 (1: AFR) — beta, has N (totalsamplesize)
#   Gurdasani Cell 2019 (1: AFR) — beta_fe, space-sep, compound snpid, no N
# No reference GWAS — first Height batch.

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
    # Karczewski NatGenet 2025 (beta, GWAS Catalog format, no N)
    "Karczewski_NatGenet_2025.Height.META.GRCh37.txt.gz|META"
    "Karczewski_NatGenet_2025.Height.EUR.GRCh37.txt.gz|EUR"
    "Karczewski_NatGenet_2025.Height.AFR.GRCh37.txt.gz|AFR"
    "Karczewski_NatGenet_2025.Height.SAS.GRCh37.txt.gz|SAS"
    # Verma Science 2024 (beta, has n)
    "Verma_Science_2024.Height.META.GRCh37.txt.gz|META"
    "Verma_Science_2024.Height.EUR.GRCh37.txt.gz|EUR"
    "Verma_Science_2024.Height.AFR.GRCh37.txt.gz|AFR"
    "Verma_Science_2024.Height.EAS.GRCh37.txt.gz|EAS"
    "Verma_Science_2024.Height.AMR.GRCh37.txt.gz|AMR"
    # Chen CellGenomics 2023 (beta, has n)
    "Chen_CellGenomics_2023.Height.EAS.GRCh37.tsv.gz|EAS"
    # Sohail Nature 2023 (beta, has n)
    "Sohail_Nature_2023.Height.AMR.GRCh37.txt.gz|AMR"
    # FernandezRhodes HGGAdv 2022 (beta + OR, no N)
    "FernandezRhodes_HGGAdv_2022.Height.AMR.GRCh37.txt.gz|AMR"
    # Sakaue NatGenet 2021 (beta + OR, no N)
    "Sakaue_NatGenet_2021.Height.EAS.GRCh37.txt.gz|EAS"
    # Graff AJHG 2023 (beta, has N as totalsamplesize)
    "Graff_AJHG_2023.Height.AFR.GRCh37.txt.gz|AFR"
    # Gurdasani Cell 2019 (beta_fe, space-sep, compound snpid, no N)
    "Gurdasani_Cell_2019.Height.AFR.GRCh37.txt.gz|AFR"
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

    # Harmonize (no reference GWAS — first Height batch)
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
echo "=== Height batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized Height files: $(ls ${HARMONIZED}/*Height*.processed.txt.gz 2>/dev/null | wc -l)"
