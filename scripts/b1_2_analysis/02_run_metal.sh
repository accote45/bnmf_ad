#!/bin/bash
# 02_run_metal.sh
# Execute METAL IVW meta-analysis of T2D + CAD EUR GWAS.
#
# Usage:
#   bash scripts/b1_2_analysis/02_run_metal.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
METAL_DIR="${PROJECT_ROOT}/results/b1_2_analysis/metal"

echo "=== B1.2 Step 2: Run METAL ==="
echo "  METAL dir: ${METAL_DIR}"

module load metal/2018-08-28

cd "${METAL_DIR}"
metal < run_metal.sh

# Post-process: rename output
METAL_OUT=$(ls -t meta_t2d_cad*.tbl 2>/dev/null | head -1)
if [ -z "${METAL_OUT}" ]; then
  echo "ERROR: METAL output not found"
  exit 1
fi

mv "${METAL_OUT}" meta_t2d_cad.txt
echo "  METAL output: ${METAL_DIR}/meta_t2d_cad.txt"

# Also rename the info file if present
METAL_INFO=$(ls -t meta_t2d_cad*.tbl.info 2>/dev/null | head -1)
if [ -n "${METAL_INFO}" ]; then
  mv "${METAL_INFO}" meta_t2d_cad.info
fi

# Quick summary
echo ""
echo "  Line count: $(wc -l < meta_t2d_cad.txt)"
echo "  Header: $(head -1 meta_t2d_cad.txt)"
echo "  Sample row: $(sed -n '2p' meta_t2d_cad.txt)"

echo ""
echo "=== Step 2 complete ==="
