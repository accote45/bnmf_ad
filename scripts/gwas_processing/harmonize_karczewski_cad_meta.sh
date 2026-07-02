#!/bin/bash
#BSUB -J harmonize_karczewski_cad_meta
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_karczewski_cad_meta_%J.out
#BSUB -e logs/lsf/harmonize_karczewski_cad_meta_%J.err

# Harmonize Karczewski NatGenet 2025 CAD META (GRCh37, no liftover needed).

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

# --- Input file ---
RAW_FNAME="Karczewski_NatGenet_2025.CAD.META.GRCh37.txt.gz"
RAW_FILE="${SUMSTATS}/${RAW_FNAME}"

if [ ! -f "${RAW_FILE}" ]; then
    echo "ERROR: ${RAW_FILE} not found"
    exit 1
fi

# --- Build check ---
BASE="Karczewski_NatGenet_2025.CAD.META.GRCh37"
BUILD_JSON="${BUILD_DIR}/${BASE}.build.json"
OUT_FILE="${HARMONIZED}/${BASE}.processed.txt.gz"
QC_JSON="${QC_DIR}/${BASE}.qc.json"

if [ ! -f "${BUILD_JSON}" ]; then
    echo "Running build check..."
    python "${SCRIPTS}/check_build.py" \
        --input-file "${RAW_FILE}" \
        --output-file "${BUILD_JSON}"
fi

# Fix Unknown builds — override to labeled build
if grep -q '"actual_build": "Unknown"' "${BUILD_JSON}" 2>/dev/null; then
    LABELED_BUILD="19"
    echo "Build detection failed — overriding to GRCh37 (labeled build)"
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

# --- Harmonize with META CAD reference ---
REF_GWAS="${HARMONIZED}/Tcheandjieu_NatureMed_2023.CAD.META.GRCh37.processed.txt.gz"

HARM_CMD="python ${SCRIPTS}/harmonize_sumstats.py \
    --input-file ${RAW_FILE} \
    --build-info ${BUILD_JSON} \
    --output-file ${OUT_FILE} \
    --qc-json ${QC_JSON} \
    --preferred-build ${PREFERRED_BUILD} \
    --n-cores ${N_CORES}"

if [ -f "${REF_GWAS}" ]; then
    HARM_CMD="${HARM_CMD} --ref-gwas ${REF_GWAS}"
    echo "Harmonizing with META CAD reference..."
else
    echo "WARNING: META CAD reference not found, harmonizing without allele alignment..."
fi

if eval "${HARM_CMD}"; then
    echo ""
    echo "=== Done ==="
    echo "  Output: ${OUT_FILE}"
    echo "  QC:     ${QC_JSON}"
else
    echo "FAILED: harmonization of ${BASE}"
    exit 1
fi
