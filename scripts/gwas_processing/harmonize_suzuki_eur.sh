#!/bin/bash
#BSUB -J harmonize_suzuki_eur
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 8:00
#BSUB -o logs/lsf/harmonize_suzuki_eur_%J.out
#BSUB -e logs/lsf/harmonize_suzuki_eur_%J.err

# Harmonize Suzuki_Nature_2024 EUR T2D GWAS so it can serve as the
# ancestry-matched reference for downstream EUR harmonizations (e.g. IGF-1).
# Uses META T2D as the allele-alignment reference.

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

RAW_FILE="${SUMSTATS}/Suzuki_Nature_2024.t2d.EUR.GRCh37.txt"
BASE="Suzuki_Nature_2024.t2d.EUR.GRCh37"
BUILD_JSON="${BUILD_DIR}/${BASE}.build.json"
OUT_FILE="${HARMONIZED}/${BASE}.processed.txt.gz"
QC_JSON="${QC_DIR}/${BASE}.qc.json"
REF_GWAS="${HARMONIZED}/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"

if [ -f "${OUT_FILE}" ]; then
    echo "SKIP (already exists): ${OUT_FILE}"
    exit 0
fi

echo "Processing: ${BASE}"

# Build check
if [ ! -f "${BUILD_JSON}" ]; then
    echo "  Running build check..."
    python "${SCRIPTS}/check_build.py" \
        --input-file "${RAW_FILE}" \
        --output-file "${BUILD_JSON}"
fi

# Fix Unknown builds
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

echo "  Harmonizing (ref: META T2D)..."
python "${SCRIPTS}/harmonize_sumstats.py" \
    --input-file "${RAW_FILE}" \
    --build-info "${BUILD_JSON}" \
    --output-file "${OUT_FILE}" \
    --qc-json "${QC_JSON}" \
    --preferred-build "${PREFERRED_BUILD}" \
    --n-cores "${N_CORES}" \
    --ref-gwas "${REF_GWAS}"

echo "  Done: ${OUT_FILE}"
echo ""
echo "=== Suzuki EUR T2D harmonization complete ==="
