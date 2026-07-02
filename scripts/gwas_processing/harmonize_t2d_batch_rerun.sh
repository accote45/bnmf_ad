#!/bin/bash
#BSUB -J harmonize_t2d_rerun
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_t2d_rerun_%J.out
#BSUB -e logs/lsf/harmonize_t2d_rerun_%J.err

# Re-run harmonization for 8 T2D GWAS files that failed due to wall-time limit
# in job 235564468. Increased wall time from 4h to 12h.
# The 4 Vujkovic files require GRCh38->19 liftover on ~37M variants each.

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

# T2D references for allele alignment (ancestry-matched)
declare -A T2D_REF=(
    ["EUR"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.EUR.GRCh37.txt.gz"
    ["AFR"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.AFR.GRCh37.txt.gz"
    ["META"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz"
    ["EAS"]="${HARMONIZED}/Suzuki_Nature_2024.t2d.EAS.GRCh37.processed.txt.gz"
)
# No SAS or AMR T2D reference available

# Only the 8 files that were not completed in the previous run
FILES=(
    # Vujkovic 2020 (OR only, has num_samples) — require GRCh38->19 liftover
    "Vujkovic_NatureGenetics_2020.t2d.AFR.GRCh37.txt.gz|AFR"
    "Vujkovic_NatureGenetics_2020.t2d.EAS.GRCh37.txt.gz|EAS"
    "Vujkovic_NatureGenetics_2020.t2d.AMR.GRCh37.txt.gz|AMR"
    "Vujkovic_NatureGenetics_2020.t2d.META.GRCh37.txt.gz|META"
    # GWAS Catalog harmonized format (no N)
    "HuertaChagoya_Diabetologia_2023.t2d.AMR.GRCh37.tsv.gz|AMR"
    "Loh_CommunBio_2022.t2d.SAS.GRCh37.tsv.gz|SAS"
    # Chen Diabetologia — METAL format, compound IDs
    "Chen_Diabetologia_2019.t2d.AFR.GRCh37.txt|AFR"
    # Morris 2012 — OR only, N_CASES+N_CONTROLS, no EAF
    "Morris_NatureGenetics_2012.t2d.EUR.GRCh37.txt|EUR"
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

    # Fix Unknown builds — override to GRCh37 (all files are labeled GRCh37)
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

    # Build harmonize command
    HARM_CMD="python ${SCRIPTS}/harmonize_sumstats.py \
        --input-file ${RAW_FILE} \
        --build-info ${BUILD_JSON} \
        --output-file ${OUT_FILE} \
        --qc-json ${QC_JSON} \
        --preferred-build ${PREFERRED_BUILD} \
        --n-cores ${N_CORES}"

    # Add ancestry-matched T2D reference if available
    if [[ -v "T2D_REF[${ANCESTRY}]" ]] && [ -f "${T2D_REF[${ANCESTRY}]}" ]; then
        HARM_CMD="${HARM_CMD} --ref-gwas ${T2D_REF[${ANCESTRY}]}"
        echo "  Harmonizing (ref: ${ANCESTRY} T2D)..."
    else
        echo "  Harmonizing (no ${ANCESTRY} T2D reference available)..."
    fi

    if eval "${HARM_CMD}"; then
        echo "  Done: ${OUT_FILE}"
    else
        echo "  FAILED: ${BASE}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== T2D batch harmonization rerun complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized files: $(ls ${HARMONIZED}/*t2d*.processed.txt.gz 2>/dev/null | wc -l)"
