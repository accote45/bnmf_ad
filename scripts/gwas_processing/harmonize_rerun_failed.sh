#!/bin/bash
#BSUB -J harmonize_rerun
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 24:00
#BSUB -o logs/lsf/harmonize_rerun_%J.out
#BSUB -e logs/lsf/harmonize_rerun_%J.err

# Re-run 15 GWAS files that previously produced 0 variants due to:
#   - Beta + OR column conflict (Sakaue, Nikpay): fixed in column_map.py
#   - All-NA SE with valid OR (Verma CAD/atherosclerosis): fixed in harmonize_sumstats.py
#
# Excludes 3 unrecoverable p-value-only files (Ward ALT/AST, Downie FastingGlucose).

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

# N overrides for files missing sample size column
declare -A N_OVERRIDES=(
    ["Sakaue_NatureGenetics_2021.Calcium.EAS"]=83980
)

# Files to re-harmonize: raw_filename|ancestry|trait_group
# trait_group is used to determine reference GWAS (CAD files get CAD ref)
FILES=(
    # Group 1: Beta + OR conflict (column_map.py fix)
    "Sakaue_NatureGenetics_2021.Calcium.EAS.GRCh37.txt.gz|EAS|calcium"
    "Sakaue_NatGenet_2021.Creatinine.EAS.GRCh37.txt.gz|EAS|creatinine"
    "Sakaue_NatGenet_2021.Basophilcount.EAS.GRCh37.txt.gz|EAS|basophil"
    "Sakaue_NatGenet_2021.Basophilcount.META.GRCh37.txt.gz|META|basophil"
    "Nikpay_NatureGenetics_2015.CAD.EUR.GRCh37.txt.gz|EUR|cad"
    # Group 2: All-NA SE with valid OR (harmonize_sumstats.py fix)
    "Verma_Science_2024.CAD.EUR.GRCh37.txt.gz|EUR|cad"
    "Verma_Science_2024.CAD.AFR.GRCh38.txt.gz|AFR|cad"
    "Verma_Science_2024.CAD.AMR.GRCh37.txt.gz|AMR|cad"
    "Verma_Science_2024.CAD.EAS.GRCh37.txt.gz|EAS|cad"
    "Verma_Science_2024.CAD.META.GRCh38.txt.gz|META|cad"
    "Verma_Science_2024.atherosclerosis.META.GRCh37.tsv.gz|META|atherosclerosis"
    "Verma_Science_2024.atherosclerosis.EUR.GRCh37.tsv.gz|EUR|atherosclerosis"
    "Verma_Science_2024.atherosclerosis.AFR.GRCh37.tsv.gz|AFR|atherosclerosis"
    "Verma_Science_2024.atherosclerosis.EAS.GRCh37.tsv.gz|EAS|atherosclerosis"
    "Verma_Science_2024.atherosclerosis.AMR.GRCh37.tsv.gz|AMR|atherosclerosis"
)

TOTAL=${#FILES[@]}
COUNT=0
FAILED=0

for ENTRY in "${FILES[@]}"; do
    IFS='|' read -r RAW_FNAME ANCESTRY TRAIT_GROUP <<< "${ENTRY}"
    COUNT=$((COUNT + 1))

    RAW_FILE="${SUMSTATS}/${RAW_FNAME}"

    # Derive base name: strip .txt.gz, .tsv.gz, or .txt extension
    BASE=$(echo "${RAW_FNAME}" | sed -E 's/\.(txt\.gz|tsv\.gz|txt|tsv)$//')
    BUILD_JSON="${BUILD_DIR}/${BASE}.build.json"
    OUT_FILE="${HARMONIZED}/${BASE}.processed.txt.gz"
    QC_JSON="${QC_DIR}/${BASE}.qc.json"

    if [ ! -f "${RAW_FILE}" ]; then
        echo "[${COUNT}/${TOTAL}] SKIP (not found): ${RAW_FNAME}"
        continue
    fi

    # Remove previous 0-variant output so we can overwrite
    rm -f "${OUT_FILE}" "${QC_JSON}"

    echo ""
    echo "[${COUNT}/${TOTAL}] Processing: ${BASE} (ancestry=${ANCESTRY})"

    # Build check (reuse existing build JSON)
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

    # Determine N override flag
    N_FLAG=""
    for key in "${!N_OVERRIDES[@]}"; do
        if [[ "${BASE}" == *"${key}"* ]]; then
            N_FLAG="--n-override ${N_OVERRIDES[${key}]}"
            echo "  Using N override: ${N_OVERRIDES[${key}]}"
            break
        fi
    done

    # Build harmonize command
    HARM_CMD="python ${SCRIPTS}/harmonize_sumstats.py \
        --input-file ${RAW_FILE} \
        --build-info ${BUILD_JSON} \
        --output-file ${OUT_FILE} \
        --qc-json ${QC_JSON} \
        --preferred-build ${PREFERRED_BUILD} \
        --n-cores ${N_CORES} ${N_FLAG}"

    # Add ancestry-matched CAD reference if this is a CAD file
    if [ "${TRAIT_GROUP}" = "cad" ] && [[ -v "CAD_REF[${ANCESTRY}]" ]] && [ -f "${CAD_REF[${ANCESTRY}]}" ]; then
        HARM_CMD="${HARM_CMD} --ref-gwas ${CAD_REF[${ANCESTRY}]}"
        echo "  Harmonizing (ref: ${ANCESTRY} CAD)..."
    else
        echo "  Harmonizing (no reference GWAS)..."
    fi

    if eval "${HARM_CMD}"; then
        echo "  Done: ${OUT_FILE}"
    else
        echo "  FAILED: ${BASE}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Rerun harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
