#!/usr/bin/env Rscript
# 00_split_ukb_samples.R
# Classify UKB individuals by ancestry via ADMIXTURE .5.Q file,
# split EUR into 80/20 train/validation, and write PLINK --keep files.
#
# Usage:
#   Rscript scripts/b1_analysis/00_split_ukb_samples.R
#   Rscript scripts/b1_analysis/00_split_ukb_samples.R --config config/b1_config.yaml

library(data.table)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)

admixture_file  <- cfg$ancestry$admixture_file
eur_threshold   <- cfg$ancestry$eur_threshold
afr_threshold   <- cfg$ancestry$afr_threshold
eas_threshold   <- cfg$ancestry$eas_threshold
sas_threshold   <- cfg$ancestry$sas_threshold
withdrawn_file  <- cfg$ancestry$withdrawn_file
train_fraction  <- cfg$split$train_fraction
seed            <- cfg$split$seed
results_dir     <- cfg$results_dir

cat("=== B1 Step 0: Split UKB Samples ===\n")
cat(sprintf("  ADMIXTURE file: %s\n", admixture_file))
cat(sprintf("  EUR threshold:  >%.2f\n", eur_threshold))
cat(sprintf("  AFR threshold:  >%.2f\n", afr_threshold))
cat(sprintf("  EAS threshold:  >%.2f\n", eas_threshold))
cat(sprintf("  SAS threshold:  >%.2f\n", sas_threshold))
cat(sprintf("  Train fraction: %.2f\n", train_fraction))
cat(sprintf("  Seed:           %d\n", seed))

# --- Output directory ---
out_dir <- if (!is.null(cfg$samples_dir)) cfg$samples_dir else file.path(results_dir, "samples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Read admixture ---
anc <- fread(admixture_file)
cat(sprintf("\nAdmixture: %d individuals, columns: %s\n",
            nrow(anc), paste(colnames(anc), collapse = ", ")))

# --- Remove withdrawn ---
if (file.exists(withdrawn_file)) {
  withdrawn <- fread(withdrawn_file, header = FALSE)$V1
  n_before <- nrow(anc)
  anc <- anc[!IID %in% withdrawn]
  cat(sprintf("Removed %d withdrawn individuals (%d remaining)\n",
              n_before - nrow(anc), nrow(anc)))
}

# --- Classify by ancestry ---
eur <- anc[EUR > eur_threshold]
afr <- anc[AFR > afr_threshold]
eas <- anc[EAS > eas_threshold]
sas <- anc[SAS > sas_threshold]

# Exclude GWAS derivation samples if specified (for split-UKB sensitivity)
if (!is.null(cfg$exclude_samples) && file.exists(cfg$exclude_samples)) {
  excl_ids <- fread(cfg$exclude_samples, header = FALSE)$V1
  n_before <- nrow(eur)
  eur <- eur[!IID %in% excl_ids]
  cat(sprintf("Excluded %d GWAS derivation samples from EUR (%d remaining)\n",
              n_before - nrow(eur), nrow(eur)))
}

cat(sprintf("\nAncestry counts:\n"))
cat(sprintf("  EUR (>%.2f): %d\n", eur_threshold, nrow(eur)))
cat(sprintf("  AFR (>%.2f): %d\n", afr_threshold, nrow(afr)))
cat(sprintf("  EAS (>%.2f): %d\n", eas_threshold, nrow(eas)))
cat(sprintf("  SAS (>%.2f): %d\n", sas_threshold, nrow(sas)))

# --- Split EUR into train/validation ---
set.seed(seed)
n_eur <- nrow(eur)
n_train <- round(n_eur * train_fraction)
train_idx <- sample(seq_len(n_eur), size = n_train, replace = FALSE)
val_idx <- setdiff(seq_len(n_eur), train_idx)

eur_train <- eur[train_idx, .(FID, IID)]
eur_val   <- eur[val_idx,   .(FID, IID)]

cat(sprintf("\nEUR split:\n"))
cat(sprintf("  Train:      %d (%.1f%%)\n", nrow(eur_train), 100 * nrow(eur_train) / n_eur))
cat(sprintf("  Validation: %d (%.1f%%)\n", nrow(eur_val), 100 * nrow(eur_val) / n_eur))

# --- Write keep files ---
write_keep <- function(dt, filename) {
  path <- file.path(out_dir, filename)
  fwrite(dt[, .(FID, IID)], path, sep = "\t", col.names = FALSE)
  cat(sprintf("  Wrote: %s (%d individuals)\n", path, nrow(dt)))
}

cat("\n--- Writing keep files ---\n")
write_keep(eur_train, "eur_train.keep")
write_keep(eur_val,   "eur_validation.keep")
write_keep(afr[, .(FID, IID)], "afr_validation.keep")
write_keep(eas[, .(FID, IID)], "eas_validation.keep")
write_keep(sas[, .(FID, IID)], "sas_validation.keep")

# --- Write summary ---
summary_lines <- c(
  sprintf("=== B1 Sample Split Summary ==="),
  sprintf("Date: %s", Sys.time()),
  sprintf("Seed: %d", seed),
  sprintf("ADMIXTURE file: %s", admixture_file),
  "",
  sprintf("EUR (>%.2f): %d total", eur_threshold, n_eur),
  sprintf("  Train: %d (%.1f%%)", nrow(eur_train), 100 * nrow(eur_train) / n_eur),
  sprintf("  Validation: %d (%.1f%%)", nrow(eur_val), 100 * nrow(eur_val) / n_eur),
  sprintf("AFR (>%.2f): %d", afr_threshold, nrow(afr)),
  sprintf("EAS (>%.2f): %d", eas_threshold, nrow(eas)),
  sprintf("SAS (>%.2f): %d", sas_threshold, nrow(sas))
)
summary_path <- file.path(out_dir, "sample_summary.txt")
writeLines(summary_lines, summary_path)
cat(sprintf("\nSummary: %s\n", summary_path))
cat("=== Done ===\n")
