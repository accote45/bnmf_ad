#!/usr/bin/env Rscript
# 01_format_gwas.R
# Convert PLINK2 .glm.linear output to the harmonized format expected
# by the existing bNMF pipeline (read_gwas / build_z_matrix).
#
# PLINK2 output columns:
#   #CHROM  POS  ID  REF  ALT  A1  ...  BETA  SE  ...  P  ...  OBS_CT
#
# Pipeline expected format:
#   VAR_ID  Effect_Allele  P_VALUE  BETA  SE  N  MAF
#
# VAR_ID = CHR_POS_sortedA1_sortedA2 (alphabetical allele sorting)
#
# Usage:
#   Rscript 01_format_gwas.R --input gwas.glm.linear --output formatted.txt.gz --trait LDL

library(data.table)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file  <- NULL
output_file <- NULL
trait_name  <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--input") {
    input_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output") {
    output_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--trait") {
    trait_name <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(input_file) || is.null(output_file)) {
  stop("Usage: Rscript 01_format_gwas.R --input <file> --output <file> [--trait <name>]")
}

cat(sprintf("Formatting PLINK2 output: %s\n", input_file))

# --- Read PLINK2 output ---
dt <- fread(input_file)

cat(sprintf("  Raw rows: %d\n", nrow(dt)))
cat(sprintf("  Columns: %s\n", paste(colnames(dt), collapse = ", ")))

# --- Filter to additive model only ---
if ("TEST" %in% colnames(dt)) {
  dt <- dt[TEST == "ADD"]
  cat(sprintf("  After ADD filter: %d\n", nrow(dt)))
}

# --- Build VAR_ID with sorted alleles ---
# Use #CHROM and POS for coordinates; REF and ALT for alleles
chrom_col <- intersect(colnames(dt), c("#CHROM", "CHROM", "CHR"))[1]
if (is.na(chrom_col)) stop("Cannot find chromosome column")

dt[, CHR := get(chrom_col)]
# Remove "chr" prefix if present
dt[, CHR := sub("^chr", "", CHR)]

# Sort alleles alphabetically (same convention as harmonize_sumstats.py)
dt[, allele1 := pmin(REF, ALT)]
dt[, allele2 := pmax(REF, ALT)]
dt[, VAR_ID := paste(CHR, POS, allele1, allele2, sep = "_")]

# --- Map to pipeline format ---
# Detect A1 frequency column (PLINK2 names it A1_FREQ or ALT_FREQS)
freq_col <- intersect(colnames(dt), c("A1_FREQ", "ALT_FREQS", "FREQ"))[1]
p_col    <- intersect(colnames(dt), c("P", "P_VALUE"))[1]
n_col    <- intersect(colnames(dt), c("OBS_CT", "N", "NMISS"))[1]

if (is.na(p_col)) stop("Cannot find P-value column")
if (is.na(n_col)) stop("Cannot find sample size column")

formatted <- dt[, .(
  VAR_ID       = VAR_ID,
  Effect_Allele = A1,
  P_VALUE      = get(p_col),
  BETA         = BETA,
  SE           = SE,
  N            = get(n_col)
)]

# Add MAF if frequency column exists
if (!is.na(freq_col)) {
  formatted[, MAF := ifelse(dt[[freq_col]] <= 0.5, dt[[freq_col]], 1 - dt[[freq_col]])]
} else {
  formatted[, MAF := NA_real_]
}

# --- Remove rows with missing key values ---
before_n <- nrow(formatted)
formatted <- formatted[!is.na(P_VALUE) & !is.na(BETA) & !is.na(SE)]
cat(sprintf("  After removing NA rows: %d (removed %d)\n",
            nrow(formatted), before_n - nrow(formatted)))

# --- Remove duplicates by VAR_ID (keep lowest P) ---
formatted <- formatted[order(P_VALUE)]
formatted <- formatted[!duplicated(VAR_ID)]
cat(sprintf("  After dedup by VAR_ID: %d\n", nrow(formatted)))

# --- Write output ---
fwrite(formatted, output_file, sep = "\t")
cat(sprintf("  Written: %s (%d variants)\n", output_file, nrow(formatted)))
