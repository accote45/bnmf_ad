#!/bin/bash
#BSUB -J harmonize_finngen_2018
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_finngen_2018_%J.out
#BSUB -e logs/lsf/harmonize_finngen_2018_%J.err

# Harmonize FinnGen 2018 EUR sumstats: T2D and CAD (GRCh37, no liftover needed).

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

# --- Define datasets: filename and corresponding META reference ---
declare -A REF_MAP
REF_MAP["FinnGen_2018.T2D.EUR.GRCh37"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"
REF_MAP["FinnGen_2018.CAD.EUR.GRCh37"]="${HARMONIZED}/Tcheandjieu_NatureMed_2023.CAD.META.GRCh37.processed.txt.gz"

FAILED=0

for BASE in "${!REF_MAP[@]}"; do
    echo ""
    echo "========================================"
    echo "Processing: ${BASE}"
    echo "========================================"

    RAW_FILE="${SUMSTATS}/${BASE}.txt.gz"
    BUILD_JSON="${BUILD_DIR}/${BASE}.build.json"
    OUT_FILE="${HARMONIZED}/${BASE}.processed.txt.gz"
    QC_JSON="${QC_DIR}/${BASE}.qc.json"
    REF_GWAS="${REF_MAP[${BASE}]}"

    if [ ! -f "${RAW_FILE}" ]; then
        echo "ERROR: ${RAW_FILE} not found"
        FAILED=$((FAILED + 1))
        continue
    fi

    # --- Build check ---
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

    # --- Harmonize ---
    HARM_CMD="python ${SCRIPTS}/harmonize_sumstats.py \
        --input-file ${RAW_FILE} \
        --build-info ${BUILD_JSON} \
        --output-file ${OUT_FILE} \
        --qc-json ${QC_JSON} \
        --preferred-build ${PREFERRED_BUILD} \
        --n-cores ${N_CORES}"

    if [ -f "${REF_GWAS}" ]; then
        HARM_CMD="${HARM_CMD} --ref-gwas ${REF_GWAS}"
        echo "Harmonizing with reference: $(basename ${REF_GWAS})"
    else
        echo "WARNING: Reference not found ($(basename ${REF_GWAS})), harmonizing without allele alignment..."
    fi

    if eval "${HARM_CMD}"; then
        echo ""
        echo "  Done: ${OUT_FILE}"
        echo "  QC:   ${QC_JSON}"
    else
        echo "FAILED: harmonization of ${BASE}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "========================================"
if [ ${FAILED} -eq 0 ]; then
    echo "All datasets harmonized successfully."
else
    echo "WARNING: ${FAILED} dataset(s) failed."
    exit 1
fi
