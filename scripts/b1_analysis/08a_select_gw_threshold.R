#!/usr/bin/env Rscript
# 08a_select_gw_threshold.R
# Select the best genome-wide (C+T) PRS p-value threshold for T2D and CAD by
# model fit (Nagelkerke R2) on the EUR TRAINING split, so the chosen threshold
# is locked in before evaluation on eur_validation in 08_joint_vs_gw_comparison.R.
#
# Why eur_train: selecting the threshold on the same data used to report
# performance (eur_validation) would optimistically bias the GW PRS. Picking it
# on eur_train avoids that leak.
#
# Inputs:  PRSice .all_score files (all 8 Pt_ thresholds, all individuals)
#          phenotype + covariate files, and the group column from cluster_prs_all.tsv
# Outputs: gw_threshold_selection.csv (per-threshold fit table; best flagged)
#
# Usage:
#   Rscript scripts/b1_analysis/08a_select_gw_threshold.R --config config/b1_config.yaml

library(data.table)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Traits to select thresholds for, and where their PRSice output lives
prsice_dir <- "results/a0_analysis/prs_ct/prsice_output"
traits <- list(
  T2D = file.path(prsice_dir, "T2D"),
  CAD = file.path(prsice_dir, "CAD")
)

# Covariate terms (identical to 08_joint_vs_gw_comparison.R)
covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))

# Group used for threshold selection
SELECTION_GROUP <- "eur_train"

cat("=== B1 Step 8a: Select GW PRS threshold by model fit (Nagelkerke R2) ===\n")
cat(sprintf("  Config:           %s\n", config_path))
cat(sprintf("  Selection group:  %s\n", SELECTION_GROUP))
cat(sprintf("  Covariates:       %s\n", paste(covar_terms, collapse = ", ")))
cat(sprintf("  Traits:           %s\n", paste(names(traits), collapse = ", ")))
cat(sprintf("  Output:           %s\n", file.path(results_dir, "gw_threshold_selection.csv")))

# Nagelkerke R2 (identical to 08_joint_vs_gw_comparison.R)
nagelkerke_r2 <- function(full_model, null_model, n) {
  ll_full <- as.numeric(logLik(full_model))
  ll_null <- as.numeric(logLik(null_model))
  cox_snell <- 1 - exp((2 / n) * (ll_null - ll_full))
  max_r2 <- 1 - exp((2 / n) * ll_null)
  cox_snell / max_r2
}

# ------------------------------------------------------------------
# Load shared inputs: phenotypes, covariates, and group membership
# ------------------------------------------------------------------
cat("\n--- Loading shared inputs ---\n")
pheno <- fread(cfg$phenotypes$phenotype_file)
covar <- fread(cfg$phenotypes$covariate_file)

# group column comes from cluster_prs_all.tsv (same source as script 08)
group_dt <- fread(file.path(results_dir, "prs", "cluster_prs_all.tsv"),
                  select = c("FID", "IID", "group"))
train_ids <- group_dt[group == SELECTION_GROUP, .(FID, IID)]
cat(sprintf("  %s individuals: %d\n", SELECTION_GROUP, nrow(train_ids)))

# ------------------------------------------------------------------
# Per-trait, per-threshold model fit
# ------------------------------------------------------------------
results_list <- list()

for (trait in names(traits)) {
  cat(sprintf("\n--- %s ---\n", trait))
  prefix <- traits[[trait]]

  score   <- fread(paste0(prefix, ".all_score"))
  num_snp <- fread(paste0(prefix, ".prsice"))   # Threshold, Num_SNP per threshold

  thr_cols <- grep("^Pt_", colnames(score), value = TRUE)

  # Assemble modelling table: outcome + all threshold scores + covariates,
  # restricted to the selection group
  dt <- merge(score, pheno[, c("FID", "IID", trait), with = FALSE],
              by = c("FID", "IID"))
  dt <- merge(dt, covar, by = c("FID", "IID"))
  dt <- merge(dt, train_ids, by = c("FID", "IID"))   # restrict to eur_train

  if ("Batch" %in% colnames(dt)) dt[, Batch := as.factor(Batch)]

  for (thr in thr_cols) {
    md <- dt[, c(trait, thr, covar_terms), with = FALSE]
    # Rename to a syntactically valid name; threshold labels like "Pt_5e-08"
    # contain a hyphen and break R formula parsing
    setnames(md, thr, "PRS")
    md <- na.omit(md)
    n <- nrow(md)
    n_cases <- sum(md[[trait]] == 1)

    # Standardize the PRS so betas are per-SD (comparable across thresholds)
    md[["PRS"]] <- scale(md[["PRS"]])[, 1]

    f_null <- as.formula(paste(trait, "~", paste(covar_terms, collapse = " + ")))
    f_full <- as.formula(paste(trait, "~", paste(c("PRS", covar_terms), collapse = " + ")))

    fit_null <- glm(f_null, data = md, family = binomial)
    fit_full <- glm(f_full, data = md, family = binomial)

    r2 <- nagelkerke_r2(fit_full, fit_null, n)
    coefs <- summary(fit_full)$coefficients
    prs_beta <- coefs["PRS", "Estimate"]
    prs_p    <- coefs["PRS", "Pr(>|z|)"]

    # Look up SNP count for this threshold (.prsice Threshold has no "Pt_" prefix)
    thr_label <- sub("^Pt_", "", thr)
    ns <- num_snp[Threshold == thr_label, Num_SNP]
    ns <- if (length(ns) == 1) ns else NA_integer_

    cat(sprintf("  %-10s N=%d cases=%d SNPs=%s  R2=%.5f  beta=%.4f  p=%.2e\n",
                thr, n, n_cases, ifelse(is.na(ns), "NA", as.character(ns)),
                r2, prs_beta, prs_p))

    results_list[[length(results_list) + 1]] <- data.table(
      trait = trait,
      threshold = thr,
      num_snp = ns,
      n = n,
      n_cases = n_cases,
      nagelkerke_r2 = r2,
      prs_beta_per_sd = prs_beta,
      prs_p = prs_p
    )
  }
}

results <- rbindlist(results_list)

# Flag the best threshold per trait (max Nagelkerke R2)
results[, best := nagelkerke_r2 == max(nagelkerke_r2), by = trait]

out_csv <- file.path(results_dir, "gw_threshold_selection.csv")
fwrite(results, out_csv)
cat(sprintf("\n--- Saved: %s ---\n", out_csv))

# Report selected thresholds
cat("\n=== Selected thresholds (best Nagelkerke R2 on eur_train) ===\n")
best <- results[best == TRUE]
for (i in seq_len(nrow(best))) {
  cat(sprintf("  %s: %s  (R2=%.5f, SNPs=%s)\n",
              best$trait[i], best$threshold[i], best$nagelkerke_r2[i],
              ifelse(is.na(best$num_snp[i]), "NA", as.character(best$num_snp[i]))))
}
cat("\n=== Done ===\n")
