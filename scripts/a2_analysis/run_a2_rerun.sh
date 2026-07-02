#!/usr/bin/env bash
# Driver: rerun A2 GTEx, Tabula Sapiens, and catlas enrichment with the new
# META cluster labels. Runs sequentially to stay under the low ulimit -u (256)
# and caps thread spawning. Logs per-script under logs/a2_rerun/.
set -u
cd /sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic

module load gcc/14.2.0 R/4.2.0 2>/dev/null

# Keep thread fan-out small (data.table / OpenMP / BLAS) given ulimit -u = 256
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export R_DATATABLE_NUM_THREADS=2

LOGDIR=logs/a2_rerun
mkdir -p "$LOGDIR"

run() {
  local name=$1 script=$2
  echo "=== [$(date +%H:%M:%S)] START $name ==="
  if Rscript "$script" > "$LOGDIR/$name.log" 2>&1; then
    echo "=== [$(date +%H:%M:%S)] DONE  $name ==="
  else
    echo "=== [$(date +%H:%M:%S)] FAIL  $name (see $LOGDIR/$name.log) ==="
  fi
}

run gtex             scripts/a2_analysis/a2_gtex_heatmaps.R
run tabula_sapiens   scripts/a2_analysis/a2_tabula_sapiens_heatmaps.R
run catlas_pancreas  scripts/a2_analysis/a2_catlas_pancreas_heatmaps.R
run catlas_esophagus scripts/a2_analysis/a2_catlas_esophagus_heatmaps.R
run catlas_wss       scripts/a2_analysis/a2_catlas_wss_summary.R

echo "=== ALL DONE [$(date +%H:%M:%S)] ==="
