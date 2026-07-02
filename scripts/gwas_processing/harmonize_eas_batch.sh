#!/bin/bash
#BSUB -J harmonize_eas
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_eas_%J.out
#BSUB -e logs/lsf/harmonize_eas_%J.err

# Harmonize all EAS GWAS summary statistics files
# Step 1: Build check + harmonize T2D EAS (used as allele alignment reference)
# Step 2: Build check + harmonize all remaining EAS trait files

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

# ============================================================
# Step 1: Harmonize T2D EAS reference (no allele alignment ref)
# ============================================================
echo "=== Step 1: Harmonizing T2D EAS reference ==="

T2D_RAW="${SUMSTATS}/Suzuki_Nature_2024.t2d.EAS.GRCh37.txt"
T2D_BUILD="${BUILD_DIR}/Suzuki_Nature_2024.t2d.EAS.GRCh37.build.json"
T2D_OUT="${HARMONIZED}/Suzuki_Nature_2024.t2d.EAS.GRCh37.processed.txt.gz"
T2D_QC="${QC_DIR}/Suzuki_Nature_2024.t2d.EAS.GRCh37.qc.json"

if [ ! -f "${T2D_OUT}" ]; then
    # Build check
    if [ ! -f "${T2D_BUILD}" ]; then
        echo "  Running build check for T2D EAS..."
        python "${SCRIPTS}/check_build.py" \
            --input-file "${T2D_RAW}" \
            --output-file "${T2D_BUILD}"
    fi

    # Harmonize (no --ref-gwas since this IS the reference)
    echo "  Harmonizing T2D EAS..."
    python "${SCRIPTS}/harmonize_sumstats.py" \
        --input-file "${T2D_RAW}" \
        --build-info "${T2D_BUILD}" \
        --output-file "${T2D_OUT}" \
        --qc-json "${T2D_QC}" \
        --preferred-build "${PREFERRED_BUILD}" \
        --n-cores "${N_CORES}"
else
    echo "  T2D EAS already harmonized: ${T2D_OUT}"
fi

# ============================================================
# Step 2: Harmonize all other EAS trait GWAS files
# ============================================================
echo ""
echo "=== Step 2: Harmonizing EAS trait GWAS files ==="

# N overrides for files missing sample size column
declare -A N_OVERRIDES=(
    ["Mishra_Nature_2022.stroke.EAS"]=254159
)

# List of all EAS raw files (excluding T2D which is already done)
EAS_FILES=(
    "Sakaue_NatGenet_2021.CAD.EAS.GRCh37.txt.gz"
    "Mishra_Nature_2022.stroke.EAS.GRCh37.harmonized.tsv.gz"
    "Enzan_NatureComms_2025.HF.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.LDL.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.HDL.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.TotalCholesterol.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.Triglycerides.EAS.GRCh37.txt.gz"
    "Sakaue_NatureGenetics_2021.BMI.EAS.GRCh37.txt.gz"
    "Sakaue_NatureGenetics_2021.SBP.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.DBP.EAS.GRCh37.txt.gz"
    "Chen_NatGenet_2021.FastingGlucose.EAS.GRCh37.txt.gz"
    "Kanai_NatureGenetics_2021.ALT.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.GGT.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.Bilirubin.EAS.GRCh37.txt.gz"
    "Sakaue_NatureGenetics_2021.ALP.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.CRP.EAS.GRCh37.txt.gz"
    "Sakaue_NatGenet_2021.MonocyteCount.EAS.GRCh37.txt.gz"
)

TOTAL=${#EAS_FILES[@]}
COUNT=0

for RAW_FNAME in "${EAS_FILES[@]}"; do
    COUNT=$((COUNT + 1))
    RAW_FILE="${SUMSTATS}/${RAW_FNAME}"

    if [ ! -f "${RAW_FILE}" ]; then
        echo "  [${COUNT}/${TOTAL}] SKIP (not found): ${RAW_FNAME}"
        continue
    fi

    # Derive base name (everything up to and including GRCh##)
    # e.g. Sakaue_NatGenet_2021.CAD.EAS.GRCh37
    BASE=$(echo "${RAW_FNAME}" | sed 's/\(GRCh[0-9]*\).*/\1/')
    BUILD_JSON="${BUILD_DIR}/${BASE}.build.json"
    OUT_FILE="${HARMONIZED}/${BASE}.processed.txt.gz"
    QC_JSON="${QC_DIR}/${BASE}.qc.json"

    if [ -f "${OUT_FILE}" ]; then
        echo "  [${COUNT}/${TOTAL}] SKIP (exists): ${BASE}"
        continue
    fi

    echo ""
    echo "  [${COUNT}/${TOTAL}] Processing: ${BASE}"

    # Build check
    if [ ! -f "${BUILD_JSON}" ]; then
        echo "    Running build check..."
        python "${SCRIPTS}/check_build.py" \
            --input-file "${RAW_FILE}" \
            --output-file "${BUILD_JSON}"
    fi

    # Fix Unknown builds: all EAS files are labeled GRCh37 (Biobank Japan = GRCh37)
    if grep -q '"actual_build": "Unknown"' "${BUILD_JSON}" 2>/dev/null; then
        echo "    Build detection failed — overriding to GRCh37 (labeled build)"
        python -c "
import json
with open('${BUILD_JSON}') as f:
    d = json.load(f)
d['actual_build'] = '19'
d['override_reason'] = 'Build detection failed; using labeled GRCh37 (Biobank Japan data)'
with open('${BUILD_JSON}', 'w') as f:
    json.dump(d, f, indent=2)
"
    fi

    # Determine N override flag
    N_FLAG=""
    for key in "${!N_OVERRIDES[@]}"; do
        if [[ "${BASE}" == *"${key}"* ]]; then
            N_FLAG="--n-override ${N_OVERRIDES[${key}]}"
            echo "    Using N override: ${N_OVERRIDES[${key}]}"
            break
        fi
    done

    # Harmonize with T2D EAS as allele alignment reference
    echo "    Harmonizing..."
    python "${SCRIPTS}/harmonize_sumstats.py" \
        --input-file "${RAW_FILE}" \
        --build-info "${BUILD_JSON}" \
        --output-file "${OUT_FILE}" \
        --qc-json "${QC_JSON}" \
        --preferred-build "${PREFERRED_BUILD}" \
        --n-cores "${N_CORES}" \
        --ref-gwas "${T2D_OUT}" \
        ${N_FLAG}

    echo "    Done: ${OUT_FILE}"
done

echo ""
echo "=== EAS harmonization complete ==="
echo "Harmonized files in: ${HARMONIZED}"
ls -1 "${HARMONIZED}"/*EAS* 2>/dev/null | wc -l
echo "EAS files harmonized"
