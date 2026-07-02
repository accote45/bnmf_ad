#!/usr/bin/env Rscript
# 00_subsample_ukb.R
# Read ADMIXTURE Q-matrix, define AFR (AFR > threshold) and EUR
# (EUR > threshold), count AFR N, draw non-overlapping EUR subsamples
# of size AFR_N, and write PLINK --keep files.
#
# Usage:
#   Rscript scripts/c1_ancestry_test/00_subsample_ukb.R
#   Rscript scripts/c1_ancestry_test/00_subsample_ukb.R --config config/c1_config.yaml

library(data.table)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/c1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

# --- Load config ---
cfg <- yaml::read_yaml(config_path)

admixture_file  <- cfg$ukb$admixture_file
afr_threshold   <- cfg$ukb$afr_threshold
eur_threshold   <- cfg$ukb$eur_threshold
n_subsamples    <- cfg$subsampling$n_subsamples
seed            <- cfg$subsampling$seed
results_dir     <- cfg$results_dir

# --- Create output directories ---
sub_dir <- file.path(results_dir, "subsamples")
afr_dir <- file.path(results_dir, "afr")
dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(afr_dir, recursive = TRUE, showWarnings = FALSE)

# --- Read admixture Q-matrix ---
cat("Reading admixture file:", admixture_file, "\n")
anc <- fread(admixture_file)

cat("\nColumn names:\n")
cat(paste(" ", colnames(anc), collapse = "\n"), "\n")

# Validate required columns
required_cols <- c("FID", "IID", "EUR", "AFR")
missing <- setdiff(required_cols, colnames(anc))
if (length(missing) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
}

# --- Extract EUR and AFR individuals by admixture proportion ---
cat(sprintf("\nAFR threshold: >%.2f\n", afr_threshold))
cat(sprintf("EUR threshold: >%.2f\n", eur_threshold))

eur <- anc[EUR > eur_threshold]
afr <- anc[AFR > afr_threshold]

eur_n <- nrow(eur)
afr_n <- nrow(afr)

cat(sprintf("\nEUR individuals: %d\n", eur_n))
cat(sprintf("AFR individuals: %d\n", afr_n))
cat(sprintf("EUR subsample size (= AFR N): %d\n", afr_n))
cat(sprintf("Number of subsamples: %d\n", n_subsamples))
cat(sprintf("Total EUR needed: %d\n", n_subsamples * afr_n))

if (eur_n < n_subsamples * afr_n) {
  stop(sprintf(
    "Not enough EUR individuals for %d non-overlapping subsamples of size %d. Have %d, need %d.",
    n_subsamples, afr_n, eur_n, n_subsamples * afr_n
  ))
}

# --- Subsample EUR ---
set.seed(seed)
shuffled_eur <- eur[sample(.N)]

cat("\nWriting EUR subsample keep files...\n")
for (i in seq_len(n_subsamples)) {
  start_idx <- (i - 1) * afr_n + 1
  end_idx   <- i * afr_n
  sub_ids   <- shuffled_eur[start_idx:end_idx]
  keep_df   <- sub_ids[, .(FID, IID)]

  # Create per-subsample directory
  sub_out_dir <- file.path(sub_dir, sprintf("sub_%02d", i))
  dir.create(sub_out_dir, recursive = TRUE, showWarnings = FALSE)

  keep_file <- file.path(sub_out_dir, sprintf("subsample_%02d.keep", i))
  fwrite(keep_df, keep_file, sep = "\t", col.names = FALSE)
  cat(sprintf("  sub_%02d: %d individuals -> %s\n", i, nrow(keep_df), keep_file))
}

# --- Write AFR keep file ---
cat("\nWriting AFR keep file...\n")
afr_keep <- afr[, .(FID, IID)]
afr_keep_file <- file.path(afr_dir, "afr_all.keep")
fwrite(afr_keep, afr_keep_file, sep = "\t", col.names = FALSE)
cat(sprintf("  AFR: %d individuals -> %s\n", nrow(afr_keep), afr_keep_file))

# --- Save metadata ---
metadata <- data.table(
  parameter = c("admixture_file", "afr_threshold", "eur_threshold",
                "n_subsamples", "eur_total", "afr_total",
                "subsample_size", "seed"),
  value = c(admixture_file, as.character(afr_threshold),
            as.character(eur_threshold), as.character(n_subsamples),
            as.character(eur_n), as.character(afr_n),
            as.character(afr_n), as.character(seed))
)
metadata_file <- file.path(results_dir, "subsample_metadata.tsv")
fwrite(metadata, metadata_file, sep = "\t")
cat(sprintf("\nMetadata saved: %s\n", metadata_file))

cat("\n=== Subsampling complete ===\n")
