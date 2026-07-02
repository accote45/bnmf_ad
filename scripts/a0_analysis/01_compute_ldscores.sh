#!/usr/bin/env bash
# 01_compute_ldscores.sh — Compute EUR LD scores from 1KG EUR PLINK files.
#
# Run once to generate per-chromosome .l2.ldscore.gz files used by LDSC.
module load python/2.7.17
# Output: reference/1kg_eur_ldscores/
#
# Prerequisites:
#   - ldsc.py in PATH (or conda activate ldsc)
#   - PLINK files at reference/1kg_eur/1000G.EUR.QC.{chr}.bed/bim/fam
#
# Usage:
#   bash scripts/a0_analysis/01_compute_ldscores.sh
#   (or submit as an LSF job — see Snakefile_a0 rule compute_ldscores)

set -euo pipefail

PLINK_PREFIX="reference/1kg_eur/1000G.EUR.QC"
HM3_SNP_LIST="reference/hapmap3/hapmap3_snps.tsv"
OUT_DIR="reference/1kg_eur_ldscores"

mkdir -p "${OUT_DIR}"

# Extract rsid list from HapMap3 for --extract (ldsc expects rsids, one per line)
HM3_RSID_LIST="${OUT_DIR}/hm3_rsids.txt"
if [ ! -f "${HM3_RSID_LIST}" ]; then
    echo "Extracting HapMap3 rsids..."
    awk 'NR > 1 && $2 != "." {print $2}' "${HM3_SNP_LIST}" > "${HM3_RSID_LIST}"
    echo "  $(wc -l < ${HM3_RSID_LIST}) rsids written to ${HM3_RSID_LIST}"
fi

for CHR in $(seq 1 22); do
    PLINK_CHR="${PLINK_PREFIX}.${CHR}"
    OUT_CHR="${OUT_DIR}/${CHR}"

    if [ -f "${OUT_CHR}.l2.ldscore.gz" ]; then
        echo "CHR ${CHR}: already computed, skipping."
        continue
    fi

    echo "Computing LD scores for CHR ${CHR}..."
    python /sc/arion/projects/psychgen/projects/prs/sample_overlap/software/ldsc/ldsc.py \
        --bfile "${PLINK_CHR}" \
        --l2 \
        --ld-wind-kb 1000 \
        --extract "${HM3_RSID_LIST}" \
        --out "${OUT_CHR}"

    echo "  Done: ${OUT_CHR}.l2.ldscore.gz"
done

echo ""
echo "All chromosomes complete. LD scores in: ${OUT_DIR}/"
