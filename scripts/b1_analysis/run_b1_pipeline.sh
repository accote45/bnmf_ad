#!/bin/bash
#BSUB -J b1_pipeline_k10
#BSUB -P acc_paul_oreilly
#BSUB -q premium
#BSUB -n 4
#BSUB -R "rusage[mem=32000]"
#BSUB -W 4:00
#BSUB -o logs/lsf/b1_pipeline_k10_%J.out
#BSUB -e logs/lsf/b1_pipeline_k10_%J.err

set -euo pipefail

PROJECT_ROOT="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
cd "${PROJECT_ROOT}"

CONFIG="config/b1_config.yaml"

echo "=== B1 Pipeline Rerun (K=10) ==="
echo "Started: $(date)"
echo "Config: ${CONFIG}"
echo ""

echo "--- Step 1/4: Compute Cluster PRS ---"
Rscript scripts/b1_analysis/01_compute_cluster_prs.R --config "${CONFIG}"
echo ""

echo "--- Step 2/4: Cluster PRS Association ---"
Rscript scripts/b1_analysis/02_cluster_prs_association.R --config "${CONFIG}"
echo ""

echo "--- Step 3/4: Forest Plots ---"
Rscript scripts/b1_analysis/03_b1_forest_plots.R --config "${CONFIG}"
echo ""

echo "--- Step 4/5: Quantile PRS Analysis (regressions) ---"
Rscript scripts/b1_analysis/05_quantile_prs_analysis.R --config "${CONFIG}"
echo ""

echo "--- Step 5/5: Quantile PRS Forest Plots ---"
Rscript scripts/b1_analysis/05b_quantile_forest_plots.R --config "${CONFIG}"
echo ""

echo "=== B1 Pipeline Complete ==="
echo "Finished: $(date)"
