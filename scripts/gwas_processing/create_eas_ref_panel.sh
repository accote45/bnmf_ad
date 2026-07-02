#!/bin/bash
#BSUB -J create_eas_ref
#BSUB -P acc_paul_oreilly
#BSUB -q express
#BSUB -n 4
#BSUB -R "rusage[mem=16000]"
#BSUB -W 02:00
#BSUB -o logs/lsf/create_eas_ref_%J.out
#BSUB -e logs/lsf/create_eas_ref_%J.err

# Create 1000G EAS reference panel for LD clumping
# Extracts EAS samples from Phase 3 pgen, creates per-chromosome plink1 bed/bim/fam

set -euo pipefail

module load plink2/v2.00a5.14
module load plink/1.90b6.21

PROJECT_ROOT="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
DOWNLOAD_DIR="${PROJECT_ROOT}/reference/download"
OUTPUT_DIR="${PROJECT_ROOT}/reference/1kg_eas"
SAMPLES_FILE="${DOWNLOAD_DIR}/eas_samples.txt"

mkdir -p "${OUTPUT_DIR}"

echo "=== Creating EAS reference panel ==="
echo "Samples: $(wc -l < ${SAMPLES_FILE}) EAS individuals"

# Decompress pvar if needed
if [ ! -f "${DOWNLOAD_DIR}/all_phase3.pvar" ]; then
    echo "Decompressing pvar.zst..."
    plink2 --zst-decompress "${DOWNLOAD_DIR}/all_phase3.pvar.zst" > "${DOWNLOAD_DIR}/all_phase3.pvar"
fi

for CHR in $(seq 1 22); do
    echo "Processing chr${CHR}..."

    # Extract EAS samples for this chromosome, convert to plink1 format
    # Apply same QC as EUR/AFR panels: SNPs only, MAF > 0.01, geno < 0.05
    plink2 \
        --pfile "${DOWNLOAD_DIR}/all_phase3" \
        --keep "${SAMPLES_FILE}" \
        --chr "${CHR}" \
        --snps-only just-acgt \
        --max-alleles 2 \
        --maf 0.01 \
        --geno 0.05 \
        --make-bed \
        --out "${OUTPUT_DIR}/1000G.EAS.QC.${CHR}" \
        --threads 4 \
        --memory 14000

    echo "  chr${CHR}: $(wc -l < ${OUTPUT_DIR}/1000G.EAS.QC.${CHR}.bim) variants"
done

echo "=== EAS reference panel complete ==="
echo "Output: ${OUTPUT_DIR}/1000G.EAS.QC.{1-22}.{bed,bim,fam}"
