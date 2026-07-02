#!/bin/bash
#BSUB -J harmonize_gurdasani_height
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=128000]"
#BSUB -W 24:00
#BSUB -o logs/lsf/harmonize_gurdasani_height_%J.out
#BSUB -e logs/lsf/harmonize_gurdasani_height_%J.err

# Re-run Gurdasani_Cell_2019.Height.AFR with 128GB memory (was OOM-killed at 64GB).
# This file has ~25M variants in a non-standard space-separated format with compound SNPIDs.

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

RAW_FNAME="Gurdasani_Cell_2019.Height.AFR.GRCh37.txt.gz"
RAW_FILE="${SUMSTATS}/${RAW_FNAME}"
BASE="Gurdasani_Cell_2019.Height.AFR.GRCh37"
BUILD_JSON="${BUILD_DIR}/${BASE}.build.json"
OUT_FILE="${HARMONIZED}/${BASE}.processed.txt.gz"
QC_JSON="${QC_DIR}/${BASE}.qc.json"

# Remove previous failed output if any
rm -f "${OUT_FILE}" "${QC_JSON}"

echo "Processing: ${BASE} (ancestry=AFR)"

# Build check (reuse existing build JSON)
if [ ! -f "${BUILD_JSON}" ]; then
    echo "  Running build check..."
    python "${SCRIPTS}/check_build.py" \
        --input-file "${RAW_FILE}" \
        --output-file "${BUILD_JSON}"
fi

# Fix Unknown builds — override to labeled build
if grep -q '"actual_build": "Unknown"' "${BUILD_JSON}" 2>/dev/null; then
    LABELED_BUILD="19"
    echo "  Build detection failed — overriding to GRCh37 (labeled build)"
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

# Harmonize
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
    exit 1
fi

echo ""
echo "=== Gurdasani Height rerun complete ==="
