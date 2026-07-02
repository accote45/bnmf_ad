#!/usr/bin/env bash
# ============================================================================
# setup_aou_workspace.sh
# Create multiancestry_polygenic directory structure within an existing
# All of Us Researcher Workbench workspace.
#
# Usage (run from AoU terminal):
#   bash setup_aou_workspace.sh
#
# Creates the tree under ~/workspaces/duplicateofstatinresponse/multiancestry_polygenic/
# ============================================================================
set -euo pipefail

BASE_DIR="${HOME}/workspaces/duplicateofstatinresponse"
PROJECT_DIR="${BASE_DIR}/multiancestry_polygenic"

echo "============================================"
echo "Creating multiancestry_polygenic directory"
echo "  at: ${PROJECT_DIR}"
echo "============================================"
echo ""

dirs=(
  "${PROJECT_DIR}/scripts"
  "${PROJECT_DIR}/config"
  "${PROJECT_DIR}/genotypes/EUR"
  "${PROJECT_DIR}/genotypes/AFR"
  "${PROJECT_DIR}/genotypes/META"
  "${PROJECT_DIR}/phenotypes"
  "${PROJECT_DIR}/prs/score_files"
  "${PROJECT_DIR}/prs/output"
  "${PROJECT_DIR}/prs/tmp"
  "${PROJECT_DIR}/tools/liftover"
)

for d in "${dirs[@]}"; do
  mkdir -p "$d"
  echo "  Created: $d"
done

echo ""
echo "============================================"
echo "Done!"
echo "============================================"
echo ""
echo "Next steps:"
echo "  1. Upload W_matrix_{ANC}.tsv + filtered_variants_{ANC}.tsv into:"
echo "       ${PROJECT_DIR}/genotypes/{EUR,AFR,META}/"
echo ""
echo "  2. Install liftOver (run in terminal):"
echo "       LIFTOVER_DIR=${PROJECT_DIR}/tools/liftover"
echo '       wget -O "$LIFTOVER_DIR/liftOver" https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver'
echo '       chmod +x "$LIFTOVER_DIR/liftOver"'
echo '       wget -O "$LIFTOVER_DIR/hg19ToHg38.over.chain.gz" https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
echo ""
echo "  3. Run PRS computation:"
echo "       Rscript scripts/compute_prs.R"
echo ""
echo "  4. Run phenotype wrangling:"
echo "       Rscript scripts/wrangle_aou_phenotypes.R"
echo ""
echo "Existing workspace resources:"
echo "  PRScsx:     ${BASE_DIR}/PRScsx/"
echo "  LD panels:  ${BASE_DIR}/genome_useful/ldblk_1kg_{afr,eur}"
echo "  Genotypes:  ${BASE_DIR}/genotype/plink/chr{1..22}.bed"
