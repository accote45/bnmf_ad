#!/bin/bash
# 05a_clump_gwas_loci.sh
# Prepare clump input files and run plink --clump for 6 traits to get
# genome-wide significant (P < 5e-8) independent lead SNPs.
#
# Usage:
#   # Step 1 — prepare input files (run interactively):
#   bash scripts/a0_analysis/05a_clump_gwas_loci.sh prep
#
#   # Step 2 — clump (LSF array):
#   bsub -J "clump[1-6]" -n 1 -R "rusage[mem=8000]" -W 1:00 \
#     -q express -P acc_paul_oreilly \
#     -o logs/clump_%J_%I.stdout -e logs/clump_%J_%I.stderr \
#     bash scripts/a0_analysis/05a_clump_gwas_loci.sh clump LSF_ARRAY
#
#   # Or clump a single trait:
#   bash scripts/a0_analysis/05a_clump_gwas_loci.sh clump t2d

set -euo pipefail

# --- Configuration ---
PROJECT_DIR="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
REF_PREFIX="${PROJECT_DIR}/reference/1kg_eur/1000G.EUR.QC"
SUMSTATS_DIR="${PROJECT_DIR}/sumstats/harmonized"
OUT_DIR="${PROJECT_DIR}/results/a0_analysis/loci_overlap"
CLUMP_INPUT_DIR="${OUT_DIR}/clump_input"
CLUMPED_DIR="${OUT_DIR}/clumped"

# Trait names (lowercase keys) and their sumstats files.
# T2D = Suzuki 2024 and CAD = Aragam 2022 (matched to the a0 heritability/rg
# analysis). Both lack usable rs-IDs in the RSID column (Suzuki RSID = ".",
# Aragam rs-IDs are a newer build), so they are mapped to 1KG rs-IDs by
# CHR:POS parsed from VAR_ID — see POS_MAPPED_TRAITS below.
declare -A TRAIT_FILES
TRAIT_FILES[t2d]="Suzuki_Nature_2024.t2d.EUR.GRCh37.processed.txt.gz"
TRAIT_FILES[cad]="Aragam_NatureGenetics_2022.CAD.EUR.GRCh37.processed.txt.gz"
TRAIT_FILES[mi]="Verma_Science_2024.MI.EUR.GRCh37.processed.txt.gz"
TRAIT_FILES[angina]="Verma_Science_2024.Angina.EUR.GRCh37.processed.txt.gz"
TRAIT_FILES[stroke]="Mishra_Nature_2022.stroke.EUR.GRCh37.processed.txt.gz"
TRAIT_FILES[pad]="Verma_Science_2024.PAD.EUR.GRCh37.processed.txt.gz"

TRAIT_ORDER=(t2d cad mi angina stroke pad)

# Traits whose clump input is built by mapping VAR_ID (CHR_POS_A1_A2) -> CHR:POS
# -> rs-ID via the 1KG lookup, rather than reading the RSID column directly.
POS_MAPPED_TRAITS=" t2d cad "

# Clumping parameters (override CLUMP_R2 via env var if needed)
CLUMP_P1="5e-8"
CLUMP_P2="1"
CLUMP_R2="${CLUMP_R2:-0.05}"   # LD-based clumping: absorb variants with r2 > 0.05
CLUMP_KB="500"

# MHC region to exclude (GRCh37)
MHC_CHR=6
MHC_START=25000000
MHC_END=35000000

# ============================================================
# MODE: prep
# ============================================================
do_prep() {
    echo "=== Preparing clump input files ==="
    mkdir -p "${CLUMP_INPUT_DIR}" "${CLUMPED_DIR}"

    # Step 1: Build CHR:POS -> rsID lookup from 1KG bim files
    local LOOKUP="${CLUMP_INPUT_DIR}/1kg_pos_to_rsid.txt"
    if [[ -f "${LOOKUP}" ]]; then
        echo "  Lookup file already exists: ${LOOKUP}"
    else
        echo "  Building CHR:POS -> rsID lookup from 1KG bim files..."
        for CHR in $(seq 1 22); do
            awk -v OFS="\t" '{print $1":"$4, $2}' "${REF_PREFIX}.${CHR}.bim"
        done | sort -k1,1 > "${LOOKUP}"
        echo "  Done. $(wc -l < "${LOOKUP}") variants in lookup."
    fi

    # Step 2: Create SNP/P files for each trait
    for TRAIT in "${TRAIT_ORDER[@]}"; do
        local INFILE="${SUMSTATS_DIR}/${TRAIT_FILES[$TRAIT]}"
        local OUTFILE="${CLUMP_INPUT_DIR}/${TRAIT}_for_clump.txt"

        if [[ -f "${OUTFILE}" ]]; then
            echo "  ${TRAIT}: clump input already exists, skipping."
            continue
        fi

        echo "  ${TRAIT}: creating clump input..."

        if [[ "${POS_MAPPED_TRAITS}" == *" ${TRAIT} "* ]]; then
            # RSID column is missing/unreliable — map to 1KG rs-IDs by position.
            # Columns: VAR_ID(1) RSID(2) Effect_Allele(3) P_VALUE(4) ...
            # VAR_ID format: CHR_POS_A1_A2 (e.g., 1_79033_A_G) -> key CHR:POS
            zcat "${INFILE}" | awk -F'\t' 'NR>1 {split($1,a,"_"); print a[1]":"a[2]"\t"$4}' | \
                sort -k1,1 | \
                join -t$'\t' -1 1 -2 1 - "${LOOKUP}" | \
                awk -F'\t' 'BEGIN{print "SNP\tP"} {print $3"\t"$2}' > "${OUTFILE}"

            local N_GWS_MAPPED=$(awk -F'\t' 'NR>1 && $2+0 < 5e-8' "${OUTFILE}" | wc -l)
            echo "    Mapped to rs-IDs by position. ${N_GWS_MAPPED} GWS SNPs (P < 5e-8) in clump input."
        else
            # CVD traits have rs-IDs in RSID column — extract directly
            zcat "${INFILE}" | awk -F'\t' 'BEGIN{print "SNP\tP"} NR>1 {print $2"\t"$4}' > "${OUTFILE}"

            local N_GWS=$(awk -F'\t' 'NR>1 && $2+0 < 5e-8' "${OUTFILE}" | wc -l)
            echo "    ${N_GWS} GWS SNPs (P < 5e-8) in clump input."
        fi
    done

    echo "=== Prep complete ==="
}

# ============================================================
# MODE: clump
# ============================================================
do_clump() {
    local TRAIT="$1"
    echo "=== Clumping ${TRAIT} ==="

    module load plink/1.90b6.21 2>/dev/null || true

    local INPUT_FILE="${CLUMP_INPUT_DIR}/${TRAIT}_for_clump.txt"
    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "ERROR: Clump input not found: ${INPUT_FILE}"
        echo "Run '05a_clump_gwas_loci.sh prep' first."
        exit 1
    fi

    # Clump per chromosome
    for CHR in $(seq 1 22); do
        local OUT_PREFIX="${CLUMPED_DIR}/${TRAIT}_chr${CHR}"

        # Skip if already done
        if [[ -f "${OUT_PREFIX}.clumped" ]]; then
            echo "  chr${CHR}: already clumped, skipping."
            continue
        fi

        echo "  chr${CHR}: clumping..."
        plink \
            --bfile "${REF_PREFIX}.${CHR}" \
            --clump "${INPUT_FILE}" \
            --clump-snp-field SNP \
            --clump-field P \
            --clump-p1 "${CLUMP_P1}" \
            --clump-p2 "${CLUMP_P2}" \
            --clump-r2 "${CLUMP_R2}" \
            --clump-kb "${CLUMP_KB}" \
            --out "${OUT_PREFIX}" \
            2>&1 || true  # plink exits non-zero if no significant SNPs on this chr
    done

    # Concatenate results, excluding MHC
    local CONCAT="${CLUMPED_DIR}/${TRAIT}_clumped.txt"
    echo "  Concatenating clumped results (excluding MHC chr6:${MHC_START}-${MHC_END})..."

    # Header: CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2
    local FIRST=true
    for CHR in $(seq 1 22); do
        local CLUMPED_FILE="${CLUMPED_DIR}/${TRAIT}_chr${CHR}.clumped"
        if [[ ! -f "${CLUMPED_FILE}" ]]; then
            continue
        fi

        if ${FIRST}; then
            # Write header from first file
            head -1 "${CLUMPED_FILE}" > "${CONCAT}"
            FIRST=false
        fi

        # Append data rows, excluding MHC
        tail -n +2 "${CLUMPED_FILE}" | \
            awk -v chr="${MHC_CHR}" -v start="${MHC_START}" -v end="${MHC_END}" \
            '{ if (!($1 == chr && $4 >= start && $4 <= end)) print }' >> "${CONCAT}"
    done

    if ${FIRST}; then
        echo "  WARNING: No clumped results found for ${TRAIT}."
    else
        local N_LEAD=$(tail -n +2 "${CONCAT}" | grep -cve '^\s*$' || true)
        echo "  ${TRAIT}: ${N_LEAD} independent lead SNPs (P < ${CLUMP_P1}, excl. MHC)"
    fi

    echo "=== ${TRAIT} clumping complete ==="
}

# ============================================================
# Main
# ============================================================
MODE="${1:-}"

if [[ "${MODE}" == "prep" ]]; then
    do_prep

elif [[ "${MODE}" == "clump" ]]; then
    TRAIT_ARG="${2:-}"

    if [[ "${TRAIT_ARG}" == "LSF_ARRAY" ]]; then
        IDX="${LSB_JOBINDEX}"
        TRAIT="${TRAIT_ORDER[$((IDX-1))]}"
        do_clump "${TRAIT}"
    elif [[ -n "${TRAIT_ARG}" ]]; then
        do_clump "${TRAIT_ARG}"
    else
        echo "Usage: $0 clump <trait|LSF_ARRAY>"
        exit 1
    fi

else
    echo "Usage: $0 <prep|clump> [trait|LSF_ARRAY]"
    exit 1
fi
