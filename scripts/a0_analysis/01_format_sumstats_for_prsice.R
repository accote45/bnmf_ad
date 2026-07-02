#!/usr/bin/env Rscript
# 01_format_sumstats_for_prsice.R
# Convert harmonized summary statistics to PRSice-compatible format.
#
# PRSice expects: SNP (rsID), A1, A2, P, BETA
# Our harmonized files have: VAR_ID, RSID, Effect_Allele, P_VALUE, BETA, SE, N, MAF, EAF
#
# Special handling for Mahajan T2D: RSID column has CHR:POS format,
# not rs-IDs. We map CHR:POS -> rsID using the imputed .bim files.
#
# Usage:
#   Rscript scripts/a0_analysis/01_format_sumstats_for_prsice.R
#   Rscript scripts/a0_analysis/01_format_sumstats_for_prsice.R \
#     --bim-dir /path/to/imputed/qced_data \
#     --out-dir results/a0_analysis/prs_ct/sumstats_prsice

library(data.table)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)

sumstats_dir <- "sumstats/harmonized"
bim_dir      <- "/sc/arion/projects/paul_oreilly/data/Biobanks/UKB/imputed/qced_data"
out_dir      <- "results/a0_analysis/prs_ct/sumstats_prsice"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--sumstats-dir") {
    sumstats_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--bim-dir") {
    bim_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--out-dir") {
    out_dir <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Trait definitions ---
traits <- list(
  T2D    = "Mahajan_NatureGenetics_2022.t2d.EUR.GRCh37.processed.txt.gz",
  CAD    = "Tcheandjieu_NatureMed_2023.CAD.EUR.GRCh37.processed.txt.gz",
  Angina = "Verma_Science_2024.Angina.EUR.GRCh37.processed.txt.gz",
  MI     = "Verma_Science_2024.MI.EUR.GRCh37.processed.txt.gz",
  Stroke = "Mishra_Nature_2022.stroke.EUR.GRCh37.processed.txt.gz",
  PAD    = "Verma_Science_2024.PAD.EUR.GRCh37.processed.txt.gz"
)

cat("=== Format Summary Statistics for PRSice ===\n")
cat(sprintf("  Sumstats dir: %s\n", sumstats_dir))
cat(sprintf("  BIM dir:      %s\n", bim_dir))
cat(sprintf("  Output dir:   %s\n\n", out_dir))

# --- Build CHR:POS -> rsID lookup from .bim files (for Mahajan T2D) ---
cat("--- Building CHR:BP -> rsID lookup from imputed .bim files ---\n")

bim_list <- lapply(1:22, function(chr) {
  bim_file <- file.path(bim_dir, sprintf("chr%d_qced.bim", chr))
  dt <- fread(bim_file, header = FALSE,
              col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  dt[, .(CHR, BP, SNP)]
})
bim_lookup <- rbindlist(bim_list)
# Create CHR:BP key for matching
bim_lookup[, chrpos := paste(CHR, BP, sep = ":")]
# Remove duplicates (keep first occurrence)
bim_lookup <- bim_lookup[!duplicated(chrpos)]
cat(sprintf("  Loaded %d variants from .bim files\n\n", nrow(bim_lookup)))

# --- Helper: extract A2 (non-effect allele) from VAR_ID ---
# VAR_ID format: CHR_POS_A_B (sorted alleles)
extract_a2 <- function(var_id, effect_allele) {
  parts <- tstrsplit(var_id, "_", fixed = TRUE)
  allele_a <- parts[[3]]
  allele_b <- parts[[4]]
  fifelse(effect_allele == allele_a, allele_b, allele_a)
}

# --- Process each trait ---
for (trait in names(traits)) {
  cat(sprintf("--- Processing %s ---\n", trait))

  infile <- file.path(sumstats_dir, traits[[trait]])
  cat(sprintf("  Input: %s\n", infile))

  dt <- fread(infile)
  cat(sprintf("  Variants read: %d\n", nrow(dt)))

  # Extract A2 from VAR_ID
  dt[, A2 := extract_a2(VAR_ID, Effect_Allele)]

  # Rename to PRSice convention
  setnames(dt, c("Effect_Allele", "P_VALUE"), c("A1", "P"))

  if (trait == "T2D") {
    # Mahajan RSID has CHR:POS format — map to rs-IDs via .bim
    cat("  Mapping CHR:POS -> rsID via .bim lookup...\n")
    dt[, chrpos := RSID]  # RSID column is already CHR:POS
    dt <- merge(dt, bim_lookup[, .(chrpos, SNP)], by = "chrpos", all.x = FALSE)
    cat(sprintf("  Variants after rsID mapping: %d\n", nrow(dt)))
  } else {
    dt[, SNP := RSID]
  }

  # Filter: require valid SNP, A1, A2, P, BETA
  dt <- dt[!is.na(SNP) & SNP != "" & !is.na(P) & !is.na(BETA)]
  # Remove SNPs without rs-IDs (e.g., indels with non-rs identifiers)
  dt <- dt[grepl("^rs", SNP)]
  cat(sprintf("  Variants after QC: %d\n", nrow(dt)))

  # Select PRSice columns
  out_dt <- dt[, .(SNP, A1, A2, P, BETA)]

  outfile <- file.path(out_dir, sprintf("%s_prsice.txt", trait))
  fwrite(out_dt, outfile, sep = "\t")
  cat(sprintf("  Output: %s (%d variants)\n\n", outfile, nrow(out_dt)))
}

cat("=== Summary statistics formatting complete ===\n")
