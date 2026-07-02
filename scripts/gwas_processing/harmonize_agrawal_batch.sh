#!/bin/bash
#BSUB -J harmonize_agrawal
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 32
#BSUB -R "rusage[mem=64000]"
#BSUB -W 12:00
#BSUB -o logs/lsf/harmonize_agrawal_%J.out
#BSUB -e logs/lsf/harmonize_agrawal_%J.err

# Harmonize 27 Agrawal NatureComms 2022 fat depot GWAS files.
# All BOLT-LMM output, META ancestry, GRCh37 — no liftover needed.
# Traits: ASAT, GFAT, VAT, ASATGFAT, VATASAT, VATGFAT
# Variants: sex-stratified (male/female) and BMI-adjusted

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

# No allele reference — no existing fat depot reference GWAS

# Files to harmonize: raw_filename|ancestry
FILES=(
    # ASAT
    "Agrawal_NatureComms_2022.ASAT.META.GRCh37.txt.gz|META"
    "Agrawal_NatureComms_2022.ASAT.META.GRCh37.female.txt.gz|META"
    "Agrawal_NatureComms_2022.ASAT.META.GRCh37.male.txt.gz|META"
    "Agrawal_NatureComms_2022.ASAT.META.GRCh37.adjbmi.txt.gz|META"
    "Agrawal_NatureComms_2022.ASAT.META.GRCh37.female.adjbmi.txt.gz|META"
    "Agrawal_NatureComms_2022.ASAT.META.GRCh37.male.adjbmi.txt.gz|META"
    # GFAT
    "Agrawal_NatureComms_2022.GFAT.META.GRCh37.txt.gz|META"
    "Agrawal_NatureComms_2022.GFAT.META.GRCh37.female.txt.gz|META"
    "Agrawal_NatureComms_2022.GFAT.META.GRCh37.male.txt.gz|META"
    "Agrawal_NatureComms_2022.GFAT.META.GRCh37.adjbmi.txt.gz|META"
    "Agrawal_NatureComms_2022.GFAT.META.GRCh37.adjbmi.female.txt.gz|META"
    "Agrawal_NatureComms_2022.GFAT.META.GRCh37.adjbmi.male.txt.gz|META"
    # VAT
    "Agrawal_NatureComms_2022.VAT.META.GRCh37.txt.gz|META"
    "Agrawal_NatureComms_2022.VAT.META.GRCh37.female.txt.gz|META"
    "Agrawal_NatureComms_2022.VAT.META.GRCh37.male.txt.gz|META"
    "Agrawal_NatureComms_2022.VAT.META.GRCh37.adjbmi.txt.gz|META"
    "Agrawal_NatureComms_2022.VAT.META.GRCh37.adjbmi.female.txt.gz|META"
    "Agrawal_NatureComms_2022.VAT.META.GRCh37.adjbmi.male.txt.gz|META"
    # ASATGFAT
    "Agrawal_NatureComms_2022.ASATGFAT.META.GRCh37.txt.gz|META"
    "Agrawal_NatureComms_2022.ASATGFAT.META.GRCh37.female.txt.gz|META"
    "Agrawal_NatureComms_2022.ASATGFAT.META.GRCh37.male.txt.gz|META"
    # VATASAT
    "Agrawal_NatureComms_2022.VATASAT.META.GRCh37.txt.gz|META"
    "Agrawal_NatureComms_2022.VATASAT.META.GRCh37.female.txt.gz|META"
    "Agrawal_NatureComms_2022.VATASAT.META.GRCh37.male.txt.gz|META"
    # VATGFAT
    "Agrawal_NatureComms_2022.VATGFAT.META.GRCh37.txt.gz|META"
    "Agrawal_NatureComms_2022.VATGFAT.META.GRCh37.female.txt.gz|META"
    "Agrawal_NatureComms_2022.VATGFAT.META.GRCh37.male.txt.gz|META"
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
        # Extract labeled build number from filename (GRCh37 -> 19, GRCh38 -> 38)
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

    # Build harmonize command (no ref GWAS for fat depot traits)
    HARM_CMD="python ${SCRIPTS}/harmonize_sumstats.py \
        --input-file ${RAW_FILE} \
        --build-info ${BUILD_JSON} \
        --output-file ${OUT_FILE} \
        --qc-json ${QC_JSON} \
        --preferred-build ${PREFERRED_BUILD} \
        --n-cores ${N_CORES}"

    echo "  Harmonizing (no reference GWAS)..."

    if eval "${HARM_CMD}"; then
        echo "  Done: ${OUT_FILE}"
    else
        echo "  FAILED: ${BASE}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Agrawal fat depot batch harmonization complete ==="
echo "  Total: ${TOTAL}, Processed: ${COUNT}, Failed: ${FAILED}"
echo "  Harmonized Agrawal files: $(ls ${HARMONIZED}/*Agrawal*.processed.txt.gz 2>/dev/null | wc -l)"
