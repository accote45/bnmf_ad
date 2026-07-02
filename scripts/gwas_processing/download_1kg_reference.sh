#!/usr/bin/env bash
# download_1kg_reference.sh
# Download 1000 Genomes Phase 3 data and create per-chromosome
# EUR and AFR plink1 reference panels for LD clumping (GRCh37).
#
# Usage:
#   bash scripts/gwas_processing/download_1kg_reference.sh
#
# Requirements: plink2, wget/curl
# Output: reference/1kg_eur/ and reference/1kg_afr/ with per-chr bed/bim/fam

set -euo pipefail

# --- Configuration ---
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
REF_DIR="${PROJECT_ROOT}/reference"
DOWNLOAD_DIR="${REF_DIR}/download"
EUR_DIR="${REF_DIR}/1kg_eur"
AFR_DIR="${REF_DIR}/1kg_afr"

# 1000G Phase 3 plink2 files from cog-genomics.org
PGEN_URL="https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"
PVAR_URL="https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1"
PSAM_URL="https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1"

# Load plink2 module
module load plink2/v2.00a5.14 2>/dev/null || module load plink2/2.3 2>/dev/null || {
    echo "ERROR: Could not load plink2 module"
    exit 1
}

echo "=== 1000G Phase 3 Reference Panel Setup ==="
echo "Project root: ${PROJECT_ROOT}"
echo "Output dirs:  ${EUR_DIR}, ${AFR_DIR}"
echo ""

# --- Step 1: Download raw files ---
mkdir -p "${DOWNLOAD_DIR}" "${EUR_DIR}" "${AFR_DIR}"

if [ ! -f "${DOWNLOAD_DIR}/all_phase3.pgen.zst" ]; then
    echo "Downloading pgen.zst (~4.5 GB)..."
    wget -O "${DOWNLOAD_DIR}/all_phase3.pgen.zst" "${PGEN_URL}"
else
    echo "pgen.zst already downloaded, skipping."
fi

if [ ! -f "${DOWNLOAD_DIR}/all_phase3.pvar.zst" ]; then
    echo "Downloading pvar.zst..."
    wget -O "${DOWNLOAD_DIR}/all_phase3.pvar.zst" "${PVAR_URL}"
else
    echo "pvar.zst already downloaded, skipping."
fi

if [ ! -f "${DOWNLOAD_DIR}/phase3_corrected.psam" ]; then
    echo "Downloading psam..."
    wget -O "${DOWNLOAD_DIR}/phase3_corrected.psam" "${PSAM_URL}"
else
    echo "psam already downloaded, skipping."
fi

# plink2 --pfile expects all_phase3.psam alongside the pgen/pvar
if [ ! -f "${DOWNLOAD_DIR}/all_phase3.psam" ]; then
    ln -sf "${DOWNLOAD_DIR}/phase3_corrected.psam" "${DOWNLOAD_DIR}/all_phase3.psam"
fi

# --- Step 2: Decompress pgen.zst ---
if [ ! -f "${DOWNLOAD_DIR}/all_phase3.pgen" ]; then
    echo "Decompressing pgen.zst..."
    plink2 --zst-decompress "${DOWNLOAD_DIR}/all_phase3.pgen.zst" > "${DOWNLOAD_DIR}/all_phase3.pgen"
else
    echo "pgen already decompressed, skipping."
fi

# --- Step 3: Create sample lists for EUR and AFR ---
echo "Creating ancestry-specific sample lists..."

# Extract EUR sample IDs (SuperPop == EUR) — IID-only format with header for plink2
echo "#IID" > "${DOWNLOAD_DIR}/eur_samples.txt"
awk -F'\t' 'NR > 1 && $5 == "EUR" {print $1}' "${DOWNLOAD_DIR}/phase3_corrected.psam" \
    >> "${DOWNLOAD_DIR}/eur_samples.txt"
EUR_N=$(awk 'NR > 1' "${DOWNLOAD_DIR}/eur_samples.txt" | wc -l)
echo "  EUR samples: ${EUR_N}"

# Extract AFR sample IDs (SuperPop == AFR)
echo "#IID" > "${DOWNLOAD_DIR}/afr_samples.txt"
awk -F'\t' 'NR > 1 && $5 == "AFR" {print $1}' "${DOWNLOAD_DIR}/phase3_corrected.psam" \
    >> "${DOWNLOAD_DIR}/afr_samples.txt"
AFR_N=$(awk 'NR > 1' "${DOWNLOAD_DIR}/afr_samples.txt" | wc -l)
echo "  AFR samples: ${AFR_N}"

# --- Step 4: Create per-chromosome plink1 files for each ancestry ---
for ANCESTRY in EUR AFR; do
    if [ "${ANCESTRY}" == "EUR" ]; then
        OUT_DIR="${EUR_DIR}"
        KEEP_FILE="${DOWNLOAD_DIR}/eur_samples.txt"
    else
        OUT_DIR="${AFR_DIR}"
        KEEP_FILE="${DOWNLOAD_DIR}/afr_samples.txt"
    fi

    echo ""
    echo "--- Processing ${ANCESTRY} ---"

    for CHR in $(seq 1 22); do
        OUTFILE="${OUT_DIR}/1000G.${ANCESTRY}.QC.${CHR}"

        if [ -f "${OUTFILE}.bed" ] && [ -f "${OUTFILE}.bim" ] && [ -f "${OUTFILE}.fam" ]; then
            echo "  Chr ${CHR}: already exists, skipping."
            continue
        fi

        echo "  Chr ${CHR}: extracting ${ANCESTRY} samples, converting to plink1..."
        plink2 \
            --pfile "${DOWNLOAD_DIR}/all_phase3" vzs \
            --chr "${CHR}" \
            --keep "${KEEP_FILE}" \
            --snps-only just-acgt \
            --max-alleles 2 \
            --rm-dup exclude-all \
            --maf 0.005 \
            --geno 0.01 \
            --hwe 1e-6 \
            --make-bed \
            --out "${OUTFILE}" \
            --threads 4 \
            --memory 8000 \
            2>&1 | tail -1
    done

    # Report total SNP count
    TOTAL_SNPS=0
    for CHR in $(seq 1 22); do
        if [ -f "${OUT_DIR}/1000G.${ANCESTRY}.QC.${CHR}.bim" ]; then
            N=$(wc -l < "${OUT_DIR}/1000G.${ANCESTRY}.QC.${CHR}.bim")
            TOTAL_SNPS=$((TOTAL_SNPS + N))
        fi
    done
    echo "  ${ANCESTRY} total SNPs across chr1-22: ${TOTAL_SNPS}"
done

# --- Step 5: Verify ---
echo ""
echo "=== Verification ==="
for ANCESTRY in EUR AFR; do
    if [ "${ANCESTRY}" == "EUR" ]; then
        DIR="${EUR_DIR}"
    else
        DIR="${AFR_DIR}"
    fi

    N_FILES=$(ls "${DIR}"/1000G.${ANCESTRY}.QC.*.bed 2>/dev/null | wc -l)
    echo "${ANCESTRY}: ${N_FILES}/22 chromosome files created in ${DIR}"
done

echo ""
echo "=== Done! ==="
echo "EUR reference: ${EUR_DIR}/1000G.EUR.QC.{1-22}.{bed,bim,fam}"
echo "AFR reference: ${AFR_DIR}/1000G.AFR.QC.{1-22}.{bed,bim,fam}"
echo ""
echo "You can now delete the download directory to save space:"
echo "  rm -rf ${DOWNLOAD_DIR}"
