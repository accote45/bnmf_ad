#!/usr/bin/env Rscript
# 00_compute_unweighted_prs.R
# Compute UNWEIGHTED cluster-specific PRS for all UKB ancestry groups.
# Each T2D risk-increasing allele counts as +1 in its assigned cluster
# (hard assignment to max-W_k cluster).
#
# Usage:
#   Rscript scripts/b2_analysis/00_compute_unweighted_prs.R
#   Rscript scripts/b2_analysis/00_compute_unweighted_prs.R --config config/b2_1_config.yaml

library(data.table)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b2_1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)

# Read paths from a supplementary b1 config for genotype and GWAS info
b1_cfg <- read_yaml("config/b1_config.yaml")

w_matrix_path    <- b1_cfg$a1_results$w_matrix
filt_var_path    <- b1_cfg$a1_results$filtered_variants
t2d_gwas_path    <- b1_cfg$gwas$t2d
geno_prefix      <- b1_cfg$genotypes$prefix
results_dir      <- cfg$results_dir

plink2_bin <- "/hpc/packages/minerva-centos7/plink2/2.3/plink2"

cat("=== B2.1 Step 0: Compute Unweighted Cluster PRS ===\n")
cat(sprintf("  W matrix:     %s\n", w_matrix_path))
cat(sprintf("  Filtered var: %s\n", filt_var_path))
cat(sprintf("  T2D GWAS:     %s\n", t2d_gwas_path))
cat(sprintf("  Genotypes:    %s\n", geno_prefix))

# --- Output directories ---
prs_dir <- file.path(results_dir, "prs")
sample_dir <- file.path(b1_cfg$results_dir, "samples")
dir.create(prs_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load W matrix ---
w_dt <- fread(w_matrix_path)
cluster_cols <- setdiff(colnames(w_dt), "VAR_ID")
cat(sprintf("\nW matrix: %d variants x %d clusters (%s)\n",
            nrow(w_dt), length(cluster_cols), paste(cluster_cols, collapse = ", ")))

# --- Load filtered variants (for RSID + Effect_Allele) ---
filt_dt <- fread(filt_var_path)
cat(sprintf("Filtered variants: %d rows\n", nrow(filt_dt)))

# --- Load T2D GWAS for BETA to determine risk direction ---
cat("\n--- Loading T2D GWAS for risk allele alignment ---\n")
t2d <- fread(t2d_gwas_path, select = c("VAR_ID", "BETA"))
setnames(t2d, "BETA", "BETA_T2D")
cat(sprintf("  T2D: %d variants\n", nrow(t2d)))

# --- Merge W matrix with variant info and T2D BETA ---
merged <- merge(w_dt, filt_dt[, .(VAR_ID, RSID, Effect_Allele)], by = "VAR_ID")
merged <- merge(merged, t2d, by = "VAR_ID", all.x = TRUE)

cat(sprintf("\nMerged: %d variants\n", nrow(merged)))
cat(sprintf("  With T2D BETA: %d\n", sum(!is.na(merged$BETA_T2D))))

# Drop variants with no T2D BETA
n_before <- nrow(merged)
merged <- merged[!is.na(BETA_T2D)]
cat(sprintf("  Dropped %d variants with no T2D BETA\n", n_before - nrow(merged)))

# --- Hard-assign each variant to its max-W_k cluster ---
cat("\n--- Hard cluster assignment (max W_k) ---\n")

# For each variant, find cluster with highest W_k
w_mat <- as.matrix(merged[, ..cluster_cols])
max_cluster_idx <- apply(w_mat, 1, which.max)

# Build binary assignment matrix
for (j in seq_along(cluster_cols)) {
  set(merged, j = cluster_cols[j], value = as.numeric(max_cluster_idx == j))
}

# Report variants per cluster
cat("  Variants per cluster:\n")
for (j in seq_along(cluster_cols)) {
  n_in <- sum(merged[[cluster_cols[j]]] == 1)
  cat(sprintf("    %s: %d variants\n", cluster_cols[j], n_in))
}

# --- Align allele to T2D risk-increasing direction ---
cat("\n--- Aligning alleles to T2D risk-increasing direction ---\n")

# Parse REF and ALT from VAR_ID (format: CHR_POS_REF_ALT)
merged[, c("v_chr", "v_pos", "v_ref", "v_alt") := tstrsplit(VAR_ID, "_", keep = 1:4)]
merged[, v_chr := as.integer(v_chr)]
merged[, v_pos := as.integer(v_pos)]

# Determine risk allele: if BETA_T2D > 0, Effect_Allele is risk-increasing; else flip
merged[, other_allele := fifelse(Effect_Allele == v_ref, v_alt, v_ref)]
merged[, risk_allele := fifelse(BETA_T2D > 0, Effect_Allele, other_allele)]

n_flipped <- sum(merged$BETA_T2D < 0)
cat(sprintf("  Alleles kept as-is (BETA_T2D > 0): %d\n", sum(merged$BETA_T2D > 0)))
cat(sprintf("  Alleles flipped (BETA_T2D < 0): %d\n", n_flipped))

# --- Save diagnostic: variant -> cluster assignment ---
diag_dt <- merged[, .(VAR_ID, RSID, Effect_Allele, risk_allele, BETA_T2D,
                       assigned_cluster = cluster_cols[max_cluster_idx])]
fwrite(diag_dt, file.path(prs_dir, "variant_cluster_assignment.tsv"), sep = "\t")
cat(sprintf("\n  Saved diagnostic: %s\n", file.path(prs_dir, "variant_cluster_assignment.tsv")))

# --- Per-chromosome score files ---
cat("\n--- Creating per-chromosome score files ---\n")

chromosomes <- 1:22
chr_score_paths <- list()
total_matched <- 0

for (chr in chromosomes) {
  chr_vars <- merged[v_chr == chr]
  if (nrow(chr_vars) == 0) next

  # Read imputed bim for this chr
  bim_path <- sprintf("%s/chr%d.bim", dirname(geno_prefix), chr)
  if (!file.exists(bim_path)) {
    bim_path <- paste0(gsub("\\{chr\\}", chr, geno_prefix), ".bim")
  }
  bim <- fread(bim_path,
               col.names = c("CHR", "SNP", "CM", "POS", "A1_bim", "A2_bim"))

  # Match by position
  chr_vars[, pos_key := v_pos]
  bim[, pos_key := POS]
  matched <- merge(chr_vars, bim[, .(pos_key, SNP)], by = "pos_key")

  if (nrow(matched) == 0) {
    cat(sprintf("  chr%d: 0/%d matched\n", chr, nrow(chr_vars)))
    next
  }

  total_matched <- total_matched + nrow(matched)

  # Write score file: ID, A1 (risk allele), K1-K6 (0 or 1)
  score_dt <- matched[, c("SNP", "risk_allele", cluster_cols), with = FALSE]
  setnames(score_dt, c("SNP", "risk_allele"), c("ID", "A1"))

  chr_score_path <- file.path(prs_dir, sprintf("score_chr%d.tsv", chr))
  fwrite(score_dt, chr_score_path, sep = "\t")
  chr_score_paths[[as.character(chr)]] <- chr_score_path

  cat(sprintf("  chr%d: %d/%d matched\n", chr, nrow(matched), nrow(chr_vars)))
}

cat(sprintf("\nTotal matched: %d / %d variants across %d chromosomes\n",
            total_matched, nrow(merged), length(chr_score_paths)))

# --- Run plink2 --score per ancestry group per chromosome ---
cat("\n--- Running plink2 --score ---\n")

keep_files <- list(
  eur_train      = file.path(sample_dir, "eur_train.keep"),
  eur_validation = file.path(sample_dir, "eur_validation.keep"),
  afr_validation = file.path(sample_dir, "afr_validation.keep"),
  eas_validation = file.path(sample_dir, "eas_validation.keep"),
  sas_validation = file.path(sample_dir, "sas_validation.keep")
)

score_col_nums <- paste(seq(3, 2 + length(cluster_cols)), collapse = ",")

all_prs_list <- list()

for (group_name in names(keep_files)) {
  keep_path <- keep_files[[group_name]]
  if (!file.exists(keep_path)) {
    cat(sprintf("  SKIP %s: keep file not found\n", group_name))
    next
  }

  cat(sprintf("\n  === %s ===\n", group_name))

  chr_sscore_files <- c()

  for (chr in names(chr_score_paths)) {
    chr_score_path <- chr_score_paths[[chr]]
    chr_geno <- sprintf("%s/chr%s", dirname(geno_prefix), chr)
    out_prefix <- file.path(prs_dir, sprintf("%s_chr%s", group_name, chr))

    cmd <- sprintf(
      "OMP_NUM_THREADS=1 %s --bfile %s --keep %s --score %s 1 2 header cols=+scoresums ignore-dup-ids --score-col-nums %s --threads 1 --out %s 2>&1",
      plink2_bin, chr_geno, keep_path, chr_score_path, score_col_nums, out_prefix
    )

    ret <- system(cmd, intern = TRUE)
    sscore_path <- paste0(out_prefix, ".sscore")
    if (file.exists(sscore_path)) {
      chr_sscore_files <- c(chr_sscore_files, sscore_path)
    } else {
      cat(sprintf("    WARNING: chr%s sscore not found\n", chr))
    }
  }

  cat(sprintf("    %d chromosome sscore files\n", length(chr_sscore_files)))

  if (length(chr_sscore_files) == 0) next

  # Merge per-chromosome scores by summing, joined by FID/IID
  group_dt <- NULL
  for (f in chr_sscore_files) {
    dt <- fread(f)
    setnames(dt, "#FID", "FID", skip_absent = TRUE)
    score_sum_cols <- grep("^SCORE.*_SUM$", colnames(dt), value = TRUE)

    chr_dt <- dt[, c("FID", "IID", score_sum_cols), with = FALSE]

    if (is.null(group_dt)) {
      group_dt <- chr_dt
    } else {
      group_dt <- merge(group_dt, chr_dt, by = c("FID", "IID"), all = TRUE,
                        suffixes = c("", ".new"))
      for (sc in score_sum_cols) {
        new_col <- paste0(sc, ".new")
        if (new_col %in% colnames(group_dt)) {
          group_dt[[sc]] <- fifelse(is.na(group_dt[[sc]]), 0, group_dt[[sc]]) +
                            fifelse(is.na(group_dt[[new_col]]), 0, group_dt[[new_col]])
          group_dt[[new_col]] <- NULL
        }
      }
    }
  }

  # Rename score columns to PRS names
  score_sum_cols <- grep("^SCORE.*_SUM$", colnames(group_dt), value = TRUE)
  prs_names <- paste0("PRS_", cluster_cols)
  setnames(group_dt, score_sum_cols, prs_names)

  group_dt[, group := group_name]
  all_prs_list[[group_name]] <- group_dt

  cat(sprintf("    %d individuals scored\n", nrow(group_dt)))
}

# --- Combine all groups ---
cat("\n--- Combining all groups ---\n")
all_prs <- rbindlist(all_prs_list)
cat(sprintf("Total individuals: %d\n", nrow(all_prs)))
cat(sprintf("Groups: %s\n", paste(unique(all_prs$group), collapse = ", ")))

# --- Write output ---
out_path <- file.path(prs_dir, "cluster_prs_unweighted.tsv")
fwrite(all_prs, out_path, sep = "\t")
cat(sprintf("\nOutput: %s\n", out_path))

# --- Summary stats ---
prs_cols <- paste0("PRS_", cluster_cols)
cat("\n--- PRS summary (all samples) ---\n")
for (pc in prs_cols) {
  vals <- all_prs[[pc]]
  cat(sprintf("  %s: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f, NAs=%d\n",
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
    cat(sprintf("    %s: mean=%.4f, sd=%.4f\n",
                pc, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE)))
  }
}

# --- Clean up per-chromosome sscore files ---
chr_files <- list.files(prs_dir, pattern = "_chr\\d+\\.sscore$", full.names = TRUE)
if (length(chr_files) > 0) file.remove(chr_files)
chr_logs <- list.files(prs_dir, pattern = "_chr\\d+\\.log$", full.names = TRUE)
if (length(chr_logs) > 0) file.remove(chr_logs)
chr_scores <- list.files(prs_dir, pattern = "^score_chr\\d+\\.tsv$", full.names = TRUE)
if (length(chr_scores) > 0) file.remove(chr_scores)

cat("\n=== Done ===\n")
