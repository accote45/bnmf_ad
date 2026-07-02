#!/bin/bash
# 04_run_prset.sh
# Run PRSet (pathway-based PRS via PRSice) on meta-analyzed T2D+CAD GWAS.
# Full EUR sample, LD clumping enabled, Synthetic phenotype for threshold optimization.
#
# Usage:
#   bash scripts/b1_2_analysis/04_run_prset.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PRSET_DIR="${PROJECT_ROOT}/results/b1_2_analysis/prset"

echo "=== B1.2 Step 4: Run PRSet ==="
echo "  Project root: ${PROJECT_ROOT}"
echo "  PRSet dir: ${PRSET_DIR}"

PRSICE_BIN="${PROJECT_ROOT}/tools/PRSice_linux"

echo ""
echo "--- Running PRSet ---"

${PRSICE_BIN} \
    --base "${PRSET_DIR}/meta_t2d_cad_prsice.txt" \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --pvalue P \
    --stat BETA \
    --chr CHR \
    --bp BP \
    --beta \
    --target "/sc/arion/projects/paul_oreilly/data/Biobanks/UKB/imputed/qced_data/chr#" \
    --keep "${PRSET_DIR}/eur_all.keep" \
    --pheno "${PRSET_DIR}/phenotypes_with_synthetic.txt" \
    --pheno-col Synthetic \
    --cov "${PROJECT_ROOT}/results/a0_analysis/prs_ct/covariates.txt" \
    --cov-col age,age2,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Batch \
    --cov-factor Batch \
    --binary-target T \
    --gtf "/sc/arion/projects/paul_oreilly/lab/kestir01/Homo_sapiens.GRCh37.75.gtf" \
    --msigdb "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/pathway_databases/msigdb/qced_data/c2.all.v2023.2.Hs.symbols.gmt_filtered.txt" \
    --wind-3 35000 \
    --wind-5 35000 \
    --print-snp \
    --thread 4 \
    --seed 47 \
    --ultra \
    --out "${PRSET_DIR}/prset_output"

echo ""
echo "--- PRSet output files ---"
ls -lh "${PRSET_DIR}/prset_output"* 2>/dev/null || echo "  No output files found"

echo ""
echo "=== Step 4 complete ==="
