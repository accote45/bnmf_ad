#!/usr/bin/env bash
# 03_run_ldsc.sh — Run LDSC genetic correlation between EUR T2D and EUR CAD.
#
# Prerequisites:
module load python/2.7.17
#   - ldsc.py in PATH (or conda activate ldsc)
#   - Munged sumstats from 02_munge_sumstats.sh
#   - LD scores computed by 01_compute_ldscores.sh (reference/1kg_eur_ldscores/)
#
# Usage:
#   bash scripts/a0_analysis/03_run_ldsc.sh
#   (or via Snakefile_a0 rule run_ldsc)

set -euo pipefail

RESULTS_DIR="results/a0_analysis"
LD_DIR="reference/1kg_eur_ldscores/"

echo "Running LDSC genetic correlation: EUR T2D × EUR CAD"

python /sc/arion/projects/psychgen/projects/prs/sample_overlap/software/ldsc/ldsc.py \
    --rg "${RESULTS_DIR}/t2d_munged.sumstats.gz,${RESULTS_DIR}/cad_munged.sumstats.gz" \
    --ref-ld-chr "${LD_DIR}" \
    --w-ld-chr   "${LD_DIR}" \
    --out        "${RESULTS_DIR}/rg_t2d_cad"

echo ""
echo "Done. Results:"
echo "  Log:  ${RESULTS_DIR}/rg_t2d_cad.log"
echo ""
echo "Key output (from log):"
grep -A 5 "Summary of Genetic Correlation Results" "${RESULTS_DIR}/rg_t2d_cad.log" || true
