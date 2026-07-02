#!/usr/bin/env Rscript
# 01_compute_cluster_prs.R
# Compute cluster-specific PRS for all UKB ancestry groups using bNMF
# W matrix weights and max(|BETA_T2D|, |BETA_CAD|) per variant.
# Runs plink2 --score per chromosome per ancestry group, merges into one file.
#
# Usage:
#   Rscript scripts/b1_analysis/01_compute_cluster_prs.R
#   Rscript scripts/b1_analysis/01_compute_cluster_prs.R --config config/b1_config.yaml

library(data.table)
library(yaml)

# Shared scorer + BETA-lookup helpers (also used by 01b_compute_genomewide_prs.R)
source("scripts/b1_analysis/prs_scoring_helpers.R")

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}
hard_assignment <- "--hard-assignment" %in% args

cfg <- read_yaml(config_path)

w_matrix_path    <- cfg$a1_results$w_matrix
filt_var_path    <- cfg$a1_results$filtered_variants
# bNMF input GWAS lists (config/a1_config.yaml ref_gwas.META). Per-variant BETA
# is the strongest-signal effect across each trait's list, so cluster PRS,
# genome-wide PRS, and bNMF clusters all draw effect sizes from the same GWAS.
t2d_files        <- unlist(cfg$gwas$t2d_files)
cad_files        <- unlist(cfg$gwas$cad_files)
geno_prefix      <- cfg$genotypes$prefix
results_dir      <- cfg$results_dir

plink2_bin <- "/hpc/packages/minerva-centos7/plink2/2.3/plink2"

cat("=== B1 Step 1: Compute Cluster PRS ===\n")
cat(sprintf("  W matrix:     %s\n", w_matrix_path))
cat(sprintf("  Filtered var: %s\n", filt_var_path))
cat(sprintf("  T2D GWAS (%d): %s\n", length(t2d_files), paste(basename(t2d_files), collapse = ", ")))
cat(sprintf("  CAD GWAS (%d): %s\n", length(cad_files), paste(basename(cad_files), collapse = ", ")))
cat(sprintf("  Genotypes:    %s\n", geno_prefix))

# --- Output directories ---
sample_dir <- if (!is.null(cfg$samples_dir)) cfg$samples_dir else file.path(results_dir, "samples")
if (hard_assignment) {
  results_dir <- file.path(results_dir, "prs_hard_assignment")
  cat("  ** HARD ASSIGNMENT MODE **\n")
  cat(sprintf("  Output dir: %s\n", results_dir))
}
prs_dir <- file.path(results_dir, "prs")
dir.create(prs_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load W matrix ---
w_dt <- fread(w_matrix_path)
cluster_cols <- setdiff(colnames(w_dt), "VAR_ID")
cat(sprintf("\nW matrix: %d variants x %d clusters (%s)\n",
            nrow(w_dt), length(cluster_cols), paste(cluster_cols, collapse = ", ")))

# --- Load filtered variants (for RSID + Effect_Allele) ---
filt_dt <- fread(filt_var_path)
cat(sprintf("Filtered variants: %d rows\n", nrow(filt_dt)))

# --- Strongest-signal BETA per variant across each trait's GWAS list ---
# (strongest_beta lives in prs_scoring_helpers.R; restrict to the W-matrix
# variant universe so we only read the relevant rows.)
cat("\n--- Loading GWAS for BETA lookup (strongest-signal across list) ---\n")
universe_ids <- w_dt$VAR_ID
t2d <- strongest_beta(universe_ids, t2d_files)[, .(VAR_ID, BETA_T2D = BETA)]
cat(sprintf("  T2D: %d/%d universe variants with BETA\n", nrow(t2d), length(universe_ids)))

cad <- strongest_beta(universe_ids, cad_files)[, .(VAR_ID, BETA_CAD = BETA)]
cat(sprintf("  CAD: %d/%d universe variants with BETA\n", nrow(cad), length(universe_ids)))

# --- Merge W matrix with variant info and GWAS betas ---
merged <- merge(w_dt, filt_dt[, .(VAR_ID, RSID, Effect_Allele)], by = "VAR_ID")
merged <- merge(merged, t2d, by = "VAR_ID", all.x = TRUE)
merged <- merge(merged, cad, by = "VAR_ID", all.x = TRUE)

cat(sprintf("\nMerged: %d variants\n", nrow(merged)))
cat(sprintf("  With T2D BETA: %d\n", sum(!is.na(merged$BETA_T2D))))
cat(sprintf("  With CAD BETA: %d\n", sum(!is.na(merged$BETA_CAD))))

# --- Select BETA with max |BETA| ---
merged[, BETA := fifelse(
  is.na(BETA_T2D) & is.na(BETA_CAD), NA_real_,
  fifelse(is.na(BETA_T2D), BETA_CAD,
  fifelse(is.na(BETA_CAD), BETA_T2D,
  fifelse(abs(BETA_T2D) >= abs(BETA_CAD), BETA_T2D, BETA_CAD)))
)]

merged[, beta_source := fifelse(
  is.na(BETA_T2D) & is.na(BETA_CAD), "none",
  fifelse(is.na(BETA_T2D), "CAD",
  fifelse(is.na(BETA_CAD), "T2D",
  fifelse(abs(BETA_T2D) >= abs(BETA_CAD), "T2D", "CAD")))
)]

cat(sprintf("  BETA source: T2D=%d, CAD=%d, none=%d\n",
            sum(merged$beta_source == "T2D"),
            sum(merged$beta_source == "CAD"),
            sum(merged$beta_source == "none")))

# Drop variants with no BETA
merged <- merged[!is.na(BETA)]
cat(sprintf("  Variants with BETA: %d\n", nrow(merged)))

# --- Compute weighted scores: BETA * W_k ---
for (k in cluster_cols) {
  merged[[k]] <- merged$BETA * merged[[k]]
}

# --- Hard assignment: zero out non-max-W clusters per variant ---
if (hard_assignment) {
  cat("\n--- Applying hard assignment (max-W cluster only) ---\n")
  w_scores <- as.matrix(merged[, ..cluster_cols])
  max_idx <- max.col(abs(w_scores), ties.method = "first")
  keep_mask <- matrix(0, nrow = nrow(merged), ncol = length(cluster_cols))
  keep_mask[cbind(seq_len(nrow(merged)), max_idx)] <- 1
  for (j in seq_along(cluster_cols)) {
    set(merged, j = cluster_cols[j], value = merged[[cluster_cols[j]]] * keep_mask[, j])
  }
  assign_counts <- table(factor(cluster_cols[max_idx], levels = cluster_cols))
  cat("  Variants assigned per cluster:\n")
  print(assign_counts)

  assign_dt <- data.table(VAR_ID = merged$VAR_ID, assigned_cluster = cluster_cols[max_idx])
  fwrite(assign_dt, file.path(prs_dir, "variant_cluster_assignments.tsv"), sep = "\t")
  cat(sprintf("  Saved: %s\n", file.path(prs_dir, "variant_cluster_assignments.tsv")))
}

# --- Filter to non-zero weights ---
keep <- rowSums(abs(merged[, ..cluster_cols])) > 0
merged <- merged[keep]
cat(sprintf("  Variants with non-zero weights: %d\n", nrow(merged)))

# --- Parse CHR and POS from VAR_ID for per-chromosome matching ---
merged[, c("v_chr", "v_pos") := tstrsplit(VAR_ID, "_", keep = 1:2)]
merged[, v_chr := as.integer(v_chr)]
merged[, v_pos := as.integer(v_pos)]

# --- Save variant-to-beta mapping for reference ---
ref_dt <- merged[, .(VAR_ID, RSID, Effect_Allele, BETA, beta_source,
                      BETA_T2D, BETA_CAD)]
fwrite(ref_dt, file.path(prs_dir, "variant_beta_mapping.tsv"), sep = "\t")

# --- Per-chromosome score files (CHR:POS match to imputed .bim) ---
cat("\n--- Creating per-chromosome score files ---\n")
chr_score_paths <- write_chr_score_files(merged, cluster_cols, geno_prefix,
                                         prs_dir, tag = "score")

# --- Run plink2 --score per ancestry group per chromosome ---
cat("\n--- Running plink2 --score ---\n")

keep_files <- list(
  eur_train      = file.path(sample_dir, "eur_train.keep"),
  eur_validation = file.path(sample_dir, "eur_validation.keep"),
  afr_validation = file.path(sample_dir, "afr_validation.keep"),
  eas_validation = file.path(sample_dir, "eas_validation.keep"),
  sas_validation = file.path(sample_dir, "sas_validation.keep")
)

all_prs <- score_by_group(chr_score_paths, length(cluster_cols),
                          paste0("PRS_", cluster_cols),
                          keep_files, geno_prefix, plink2_bin, prs_dir)

cat("\n--- Combining all groups ---\n")
cat(sprintf("Total individuals: %d\n", nrow(all_prs)))
cat(sprintf("Groups: %s\n", paste(unique(all_prs$group), collapse = ", ")))

# --- Write output ---
out_path <- file.path(prs_dir, "cluster_prs_all.tsv")
fwrite(all_prs, out_path, sep = "\t")
cat(sprintf("\nOutput: %s\n", out_path))

# --- Summary stats ---
prs_cols <- paste0("PRS_", cluster_cols)
cat("\n--- PRS summary (all samples) ---\n")
for (pc in prs_cols) {
  vals <- all_prs[[pc]]
  cat(sprintf("  %s: mean=%.6f, sd=%.6f, min=%.6f, max=%.6f, NAs=%d\n",
              pc, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
              min(vals, na.rm = TRUE), max(vals, na.rm = TRUE),
              sum(is.na(vals))))
}

# Per-group summary
cat("\n--- PRS summary per group ---\n")
for (g in unique(all_prs$group)) {
  g_dt <- all_prs[group == g]
  cat(sprintf("  %s (n=%d):\n", g, nrow(g_dt)))
  for (pc in prs_cols) {
    vals <- g_dt[[pc]]
    cat(sprintf("    %s: mean=%.6f, sd=%.6f\n",
                pc, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE)))
  }
}

# --- Clean up per-chromosome logs ---
# Keep score_chr*.tsv, *_chr*.sscore, and *.cachekey files so subsequent runs
# can skip plink2 scoring for chromosomes whose inputs are unchanged (see the
# score cache above). Only the (tiny, disposable) plink2 logs are removed.
chr_logs <- list.files(prs_dir, pattern = "_chr\\d+\\.log$", full.names = TRUE)
if (length(chr_logs) > 0) file.remove(chr_logs)

cat("\n=== Done ===\n")
