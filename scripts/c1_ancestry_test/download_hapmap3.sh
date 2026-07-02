#!/bin/bash
# download_hapmap3.sh
# Download the LDpred2/bigsnpr HapMap3 SNP list (GRCh37) and convert to
# the project's canonical VAR_ID format (CHR_POS_sortedA1_sortedA2).
#
# Source: Prive et al. (2022) figshare — "map_hm3_plus" variant list (RDS format)
# ~1.05M biallelic autosomal SNPs well-imputed across ancestries.
#
# Usage:
#   bash scripts/c1_ancestry_test/download_hapmap3.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
HM3_DIR="${PROJECT_ROOT}/reference/hapmap3"
mkdir -p "${HM3_DIR}"

RAW_FILE="${HM3_DIR}/map_hm3_plus.rds"
OUTPUT_FILE="${HM3_DIR}/hapmap3_snps.tsv"

echo "=== Downloading HapMap3 SNP list ==="
echo "  Output dir: ${HM3_DIR}"

# --- Download from figshare ---
if [ -f "${OUTPUT_FILE}" ]; then
  echo "  Output already exists: ${OUTPUT_FILE}"
  echo "  Skipping download. Delete the file to re-download."
  N=$(wc -l < "${OUTPUT_FILE}")
  echo "  Existing file: ${N} lines"
  exit 0
fi

echo "  Downloading RDS from figshare (ndownloader)..."
curl -L -o "${RAW_FILE}" \
  "https://ndownloader.figshare.com/files/25503788"

echo "  Downloaded: $(ls -lh "${RAW_FILE}" | awk '{print $5}')"

# --- Convert to VAR_ID format ---
echo "  Converting to VAR_ID format..."

module load R/4.2.0 2>/dev/null || true

Rscript -e '
library(data.table)

hm3 <- as.data.table(readRDS("'"${RAW_FILE}"'"))
cat(sprintf("  Raw HapMap3 rows: %d\n", nrow(hm3)))
cat(sprintf("  Columns: %s\n", paste(colnames(hm3), collapse = ", ")))

# Restrict to autosomes
hm3 <- hm3[chr %in% 1:22]
cat(sprintf("  After autosome filter: %d\n", nrow(hm3)))

# Build VAR_ID with alphabetically sorted alleles (project convention)
hm3[, a0_str := as.character(a0)]
hm3[, a1_str := as.character(a1)]
hm3[, allele_min := pmin(a0_str, a1_str)]
hm3[, allele_max := pmax(a0_str, a1_str)]
hm3[, VAR_ID := paste(chr, pos, allele_min, allele_max, sep = "_")]

# Keep useful columns
out <- hm3[, .(VAR_ID, rsid, chr, pos, a0 = a0_str, a1 = a1_str)]

# Remove any duplicates
n_before <- nrow(out)
out <- out[!duplicated(VAR_ID)]
cat(sprintf("  After dedup: %d (removed %d)\n", nrow(out), n_before - nrow(out)))

fwrite(out, "'"${OUTPUT_FILE}"'", sep = "\t")
cat(sprintf("\n  Written: %s (%d SNPs)\n", "'"${OUTPUT_FILE}"'", nrow(out)))
'

# Clean up raw file
rm -f "${RAW_FILE}"

echo ""
echo "=== HapMap3 download complete ==="
echo "  File: ${OUTPUT_FILE}"
echo "  Use in config: hapmap3_snp_file: \"reference/hapmap3/hapmap3_snps.tsv\""
