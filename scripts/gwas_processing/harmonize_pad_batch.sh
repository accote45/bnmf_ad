#!/bin/bash
# Parallel launcher: submits 4 separate LSF jobs for PAD harmonization
# Usage: bash scripts/gwas_processing/harmonize_pad_batch.sh

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

mkdir -p "${HARMONIZED}" "${BUILD_DIR}" "${QC_DIR}" logs/lsf

# T2D references for allele alignment (by ancestry)
declare -A T2D_REF=(
    ["EUR"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.EUR.GRCh37.txt.gz"
    ["AFR"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.AFR.GRCh37.txt.gz"
    ["META"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"
    ["EAS"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.EAS.GRCh37.processed.txt.gz"
)

# PAD files to harmonize: raw_filename|ancestry
FILES=(
    "Verma_Science_2024.PAD.EUR.GRCh37.txt.gz|EUR"
    "Verma_Science_2024.PAD.AFR.GRCh37.txt.gz|AFR"
    "Verma_Science_2024.PAD.META.GRCh37.txt.gz|META"
    "Sakaue_NatGenet_2021.PAD.EAS.GRCh37.txt.gz|EAS"
)

for ENTRY in "${FILES[@]}"; do
    IFS='|' read -r RAW_FNAME ANCESTRY <<< "${ENTRY}"

    RAW_FILE="${SUMSTATS}/${RAW_FNAME}"
    BASE=$(echo "${RAW_FNAME}" | sed 's/\(GRCh[0-9]*\).*/\1/')
    BUILD_JSON="${BUILD_DIR}/${BASE}.build.json"
    OUT_FILE="${HARMONIZED}/${BASE}.processed.txt.gz"
    QC_JSON="${QC_DIR}/${BASE}.qc.json"
    REF_GWAS="${T2D_REF[${ANCESTRY}]}"

    # Short name for job
    SHORT=$(echo "${BASE}" | sed 's/.*2024\.\|.*2021\.//' | sed 's/\.GRCh[0-9]*//')

    if [ -f "${OUT_FILE}" ]; then
        echo "SKIP (exists): ${BASE}"
        continue
    fi

    if [ ! -f "${RAW_FILE}" ]; then
        echo "SKIP (not found): ${RAW_FNAME}"
        continue
    fi

    echo "Submitting: ${BASE}"
    bsub -J "harm_${SHORT}" \
         -P acc_paul_oreilly \
         -q premium \
         -n 32 \
         -R "rusage[mem=64000]" \
         -W 4:00 \
         -o "logs/lsf/harm_${SHORT}_%J.out" \
         -e "logs/lsf/harm_${SHORT}_%J.err" \
         bash -c "
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
cd ${PROJECT_ROOT}

# Build check
if [ ! -f '${BUILD_JSON}' ]; then
    echo 'Running build check for ${BASE}...'
    python '${SCRIPTS}/check_build.py' \
        --input-file '${RAW_FILE}' \
        --output-file '${BUILD_JSON}'
fi

# Harmonize
echo 'Harmonizing ${BASE} (ref: ${ANCESTRY} T2D)...'
python '${SCRIPTS}/harmonize_sumstats.py' \
    --input-file '${RAW_FILE}' \
    --build-info '${BUILD_JSON}' \
    --output-file '${OUT_FILE}' \
    --qc-json '${QC_JSON}' \
    --preferred-build '${PREFERRED_BUILD}' \
    --n-cores '${N_CORES}' \
    --ref-gwas '${REF_GWAS}'

echo 'Done: ${OUT_FILE}'
"
done

echo ""
echo "All jobs submitted. Monitor with: bjobs -J 'harm_*'"
