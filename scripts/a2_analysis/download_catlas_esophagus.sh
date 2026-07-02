#!/bin/bash
# Download CATLAS esophageal cCRE BED files and liftOver from hg38 → hg19
#
# Usage:
#   bash scripts/a2_analysis/download_catlas_esophagus.sh

set -euo pipefail

BASE_DIR="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
HG38_DIR="${BASE_DIR}/data/catlas/esophagus/hg38"
HG19_DIR="${BASE_DIR}/data/catlas/esophagus/hg19"
CHAIN_FILE="${BASE_DIR}/reference/hg38ToHg19.over.chain.gz"

CATLAS_URL="https://catlas.org/humanenhancer/data/cCREs"

FILES=(
  "Esophageal_Epithelial_Cell.bed"
  "Smooth_Muscle_Esophageal_Mucosal.bed"
  "Smooth_Muscle_Esophageal_Muscularis_1.bed"
  "Smooth_Muscle_Esophageal_Muscularis_2.bed"
  "Smooth_Muscle_Esophageal_Muscularis_3.bed"
)

mkdir -p "${HG38_DIR}" "${HG19_DIR}"

# --- Step 1: Check chain file ---
if [ ! -f "${CHAIN_FILE}" ]; then
  echo "Downloading hg38ToHg19 chain file..."
  wget -q -O "${CHAIN_FILE}" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
  echo "  Saved: ${CHAIN_FILE}"
else
  echo "Chain file already exists: ${CHAIN_FILE}"
fi

# --- Step 2: Download BED files ---
echo ""
echo "Downloading ${#FILES[@]} esophageal cCRE BED files..."
for f in "${FILES[@]}"; do
  out="${HG38_DIR}/${f}"
  if [ -f "${out}" ]; then
    echo "  Already exists: ${f}"
  else
    echo "  Downloading: ${f}"
    wget -q -O "${out}" "${CATLAS_URL}/${f}"
  fi
done

# --- Step 3: LiftOver hg38 → hg19 ---
module load liftover 2>/dev/null || true

echo ""
echo "LiftOver hg38 → hg19..."
for f in "${FILES[@]}"; do
  input="${HG38_DIR}/${f}"
  output="${HG19_DIR}/${f}"
  unmapped="${HG19_DIR}/${f%.bed}_unmapped.bed"

  if [ -f "${output}" ]; then
    echo "  Already lifted: ${f}"
    continue
  fi

  liftOver "${input}" "${CHAIN_FILE}" "${output}" "${unmapped}"
  n_in=$(wc -l < "${input}")
  n_out=$(wc -l < "${output}")
  n_fail=$(grep -c "^[^#]" "${unmapped}" 2>/dev/null || echo 0)
  echo "  ${f}: ${n_in} → ${n_out} mapped, ${n_fail} unmapped"
done

echo ""
echo "Done. hg19 BED files in: ${HG19_DIR}"
