#!/usr/bin/env bash
# 02_munge_sumstats.sh — Munge preprocessed T2D and CAD sumstats for LDSC.
#
# Converts LDSC-ready text files to .sumstats.gz format required by ldsc.py --rg.
module load python/2.7.17
#
# Prerequisites:
#   - ldsc.py in PATH (or conda activate ldsc)
#   - Preprocessed files from 00_preprocess_sumstats.py (supports N traits via config)
#
# Usage:
#   bash scripts/a0_analysis/02_munge_sumstats.sh
#   (or via Snakefile_a0 rule munge_sumstats)

set -euo pipefail

RESULTS_DIR="results/a0_analysis"

# ---------------------------------------------------------------------------
# Munge T2D
# ---------------------------------------------------------------------------
echo "Munging T2D sumstats..."
python /sc/arion/projects/psychgen/projects/prs/sample_overlap/software/ldsc/munge_sumstats.py \
    --sumstats "${RESULTS_DIR}/t2d_for_munge.txt" \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --signed-sumstats BETA,0 \
    --p P \
    --se SE \
    --n-col N \
    --out "${RESULTS_DIR}/t2d_munged"

echo "  Done: ${RESULTS_DIR}/t2d_munged.sumstats.gz"

# ---------------------------------------------------------------------------
# Munge CAD
# ---------------------------------------------------------------------------
echo "Munging CAD sumstats..."
python /sc/arion/projects/psychgen/projects/prs/sample_overlap/software/ldsc/munge_sumstats.py \
    --sumstats "${RESULTS_DIR}/cad_for_munge.txt" \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --signed-sumstats BETA,0 \
    --p P \
    --se SE \
    --n-col N \
    --out "${RESULTS_DIR}/cad_munged"

echo "  Done: ${RESULTS_DIR}/cad_munged.sumstats.gz"
echo ""
echo "Munge complete. Log files:"
echo "  ${RESULTS_DIR}/t2d_munged.log"
echo "  ${RESULTS_DIR}/cad_munged.log"
