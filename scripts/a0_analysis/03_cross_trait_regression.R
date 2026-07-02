#!/usr/bin/env Rscript
# 03_cross_trait_regression.R
# Cross-trait logistic regressions: T2D ~ CVD PRS and CVD ~ T2D PRS,
# adjusted for age, age^2, sex, PC1-PC10, and batch.
#
# 10 models total:
#   T2D ~ CAD_PRS, T2D ~ Angina_PRS, T2D ~ MI_PRS, T2D ~ Stroke_PRS, T2D ~ PAD_PRS
#   CAD ~ T2D_PRS, Angina ~ T2D_PRS, MI ~ T2D_PRS, Stroke ~ T2D_PRS, PAD ~ T2D_PRS
#
# PRS at p < 0.05 threshold (primary), with sensitivity at other thresholds.
#
# Usage:
#   Rscript scripts/a0_analysis/03_cross_trait_regression.R
#   Rscript scripts/a0_analysis/03_cross_trait_regression.R \
#     --prs-dir results/a0_analysis/prs_ct/prsice_output \
#     --pheno results/a0_analysis/prs_ct/phenotypes_combined.txt \
#     --covar results/a0_analysis/prs_ct/covariates.txt \
#     --out-dir results/a0_analysis/prs_ct

library(tidyverse)
library(data.table)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)

prs_dir    <- "results/a0_analysis/prs_ct/prsice_output"
pheno_file <- "results/a0_analysis/prs_ct/phenotypes_combined.txt"
covar_file <- "results/a0_analysis/prs_ct/covariates.txt"
out_dir    <- "results/a0_analysis/prs_ct"
threshold  <- 0.05  # primary p-value threshold

i <- 1
while (i <= length(args)) {
  if (args[i] == "--prs-dir") {
    prs_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pheno") {
    pheno_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--covar") {
    covar_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--out-dir") {
    out_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--threshold") {
    threshold <- as.numeric(args[i + 1]); i <- i + 2
  } else {
    i <- i + 1
  }
}

cat("=== Cross-Trait PRS Regression ===\n")
cat(sprintf("  PRS dir:     %s\n", prs_dir))
cat(sprintf("  Pheno file:  %s\n", pheno_file))
cat(sprintf("  Covar file:  %s\n", covar_file))
cat(sprintf("  Output dir:  %s\n", out_dir))
cat(sprintf("  P threshold: %g\n\n", threshold))

# --- Nagelkerke R-squared ---
nagelkerke_r2 <- function(full_model, null_model, n) {
  ll_full <- as.numeric(logLik(full_model))
  ll_null <- as.numeric(logLik(null_model))
  cox_snell <- 1 - exp((2 / n) * (ll_null - ll_full))
  max_r2 <- 1 - exp((2 / n) * ll_null)
  cox_snell / max_r2
}

# --- Load data ---
cat("--- Loading phenotypes ---\n")
pheno <- fread(pheno_file)
cat(sprintf("  %d individuals, %d columns\n", nrow(pheno), ncol(pheno)))

cat("--- Loading covariates ---\n")
covar <- fread(covar_file)
cat(sprintf("  %d individuals\n", nrow(covar)))

# --- Load PRS scores ---
# PRSice --fastscore --no-regress outputs: FID, IID, then one column per threshold
# Column names are the thresholds (e.g., "5e-08", "1e-05", "0.001", "0.01", "0.05", "0.1", "0.5", "1")
cat("\n--- Loading PRS scores ---\n")

traits <- c("T2D", "CAD", "Angina", "MI", "Stroke", "PAD")
prs_list <- list()

for (trait in traits) {
  score_file <- file.path(prs_dir, sprintf("%s.all_score", trait))
  if (!file.exists(score_file)) {
    cat(sprintf("  WARNING: %s not found, skipping %s\n", score_file, trait))
    next
  }
  dt <- fread(score_file)
  cat(sprintf("  %s: %d individuals, thresholds: %s\n",
              trait, nrow(dt), paste(colnames(dt)[-(1:2)], collapse = ", ")))

  # Find the column closest to our target threshold
  # PRSice columns may be named "Pt_0.05" or plain "0.05"
  thresh_cols <- colnames(dt)[-(1:2)]
  thresh_vals <- as.numeric(sub("^Pt_", "", thresh_cols))
  best_idx <- which.min(abs(thresh_vals - threshold))
  best_col <- thresh_cols[best_idx]
  cat(sprintf("    Using threshold: %s\n", best_col))

  prs_col_name <- paste0(trait, "_PRS")
  prs_dt <- dt[, .(FID, IID)]
  prs_dt[[prs_col_name]] <- dt[[best_col]]

  # Standardize PRS (mean=0, sd=1) for interpretable OR
  prs_dt[[prs_col_name]] <- scale(prs_dt[[prs_col_name]])[, 1]

  prs_list[[trait]] <- prs_dt
}

if (length(prs_list) == 0) {
  stop("No PRS score files found. Run PRSice (step 02) first.")
}

# --- Merge all data ---
cat("\n--- Merging data ---\n")

# Start with phenotypes + covariates
merged <- inner_join(as_tibble(pheno), as_tibble(covar), by = c("FID", "IID"))

# Add each PRS
for (trait in names(prs_list)) {
  merged <- inner_join(merged, as_tibble(prs_list[[trait]]), by = c("FID", "IID"))
}

cat(sprintf("  Merged dataset: %d individuals\n", nrow(merged)))

# Convert Batch to factor
merged <- merged %>% mutate(Batch = as.factor(Batch))

# --- Define regression models ---
cvd_traits <- c("CAD", "Angina", "MI", "Stroke", "PAD")

models <- bind_rows(
  # T2D as outcome, CVD PRS as predictor
  tibble(outcome = "T2D", predictor_prs = paste0(cvd_traits, "_PRS"),
         predictor_trait = cvd_traits),
  # CVD as outcome, T2D PRS as predictor
  tibble(outcome = cvd_traits, predictor_prs = "T2D_PRS",
         predictor_trait = "T2D")
)

# --- Covariate formula (shared across all models) ---
covar_terms <- c("age", "age2", "sex",
                 paste0("PC", 1:10),
                 "Batch")

# --- Run regressions ---
cat("\n--- Running logistic regressions ---\n")

results <- models %>%
  pmap_dfr(function(outcome, predictor_prs, predictor_trait) {
    # Check that both columns exist
    if (!outcome %in% colnames(merged) || !predictor_prs %in% colnames(merged)) {
      cat(sprintf("  SKIP: %s ~ %s (missing data)\n", outcome, predictor_prs))
      return(tibble())
    }

    # Build formula
    rhs <- paste(c(predictor_prs, covar_terms), collapse = " + ")
    formula_full <- as.formula(paste(outcome, "~", rhs))
    formula_null <- as.formula(paste(outcome, "~", paste(covar_terms, collapse = " + ")))

    # Subset to non-missing
    model_data <- merged %>%
      select(all_of(c(outcome, predictor_prs, covar_terms))) %>%
      drop_na()

    n_cases <- sum(model_data[[outcome]] == 1)
    n_controls <- sum(model_data[[outcome]] == 0)

    cat(sprintf("  %s ~ %s: %d cases, %d controls\n",
                outcome, predictor_prs, n_cases, n_controls))

    # Fit models
    fit_full <- glm(formula_full, data = model_data, family = binomial)
    fit_null <- glm(formula_null, data = model_data, family = binomial)

    # Extract PRS coefficient
    coef_summary <- summary(fit_full)$coefficients
    prs_row <- coef_summary[predictor_prs, ]

    or <- exp(prs_row["Estimate"])
    ci_lower <- exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"])
    ci_upper <- exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"])
    p_val <- prs_row["Pr(>|z|)"]

    # Nagelkerke R-squared (incremental: full vs null)
    r2 <- nagelkerke_r2(fit_full, fit_null, nrow(model_data))

    cat(sprintf("    OR=%.3f [%.3f-%.3f], P=%.2e, R2=%.4f\n",
                or, ci_lower, ci_upper, p_val, r2))

    tibble(
      outcome = outcome,
      predictor_prs = predictor_prs,
      predictor_trait = predictor_trait,
      OR = round(or, 4),
      CI_lower = round(ci_lower, 4),
      CI_upper = round(ci_upper, 4),
      p_value = p_val,
      nagelkerke_r2 = round(r2, 6),
      n_cases = n_cases,
      n_controls = n_controls,
      n_total = n_cases + n_controls,
      prs_threshold = threshold
    )
  })

# --- Write results ---
cat("\n--- Writing results ---\n")

out_file <- file.path(out_dir, "cross_trait_regression_results.csv")
write_csv(results, out_file)
cat(sprintf("  Results: %s (%d models)\n", out_file, nrow(results)))

# Print summary table
cat("\n=== Results Summary ===\n")
results %>%
  select(outcome, predictor_trait, OR, CI_lower, CI_upper, p_value, nagelkerke_r2,
         n_cases, n_controls) %>%
  print(n = Inf, width = Inf)

cat("\n=== Cross-trait regression complete ===\n")
