#!/usr/bin/env Rscript
# 02_cluster_prs_association.R
# Test associations between bNMF cluster PRS and T2D/CAD/synthetic outcomes.
# Runs logistic regressions on train set, evaluates on EUR val + AFR/EAS/SAS val.
#
# Usage:
#   Rscript scripts/b1_analysis/02_cluster_prs_association.R
#   Rscript scripts/b1_analysis/02_cluster_prs_association.R --config config/b1_config.yaml

library(data.table)
library(yaml)
library(pROC)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}
hard_assignment <- "--hard-assignment" %in% args

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
if (hard_assignment) {
  results_dir <- file.path(results_dir, "prs_hard_assignment")
  cat("  ** HARD ASSIGNMENT MODE **\n")
  cat(sprintf("  Output dir: %s\n", results_dir))
}
pheno_file  <- cfg$phenotypes$phenotype_file
covar_file  <- cfg$phenotypes$covariate_file

cat("=== B1 Step 2: Cluster PRS Association ===\n")
cat(sprintf("  PRS file:   %s\n", file.path(results_dir, "prs/cluster_prs_all.tsv")))
cat(sprintf("  Phenotypes: %s\n", pheno_file))
cat(sprintf("  Covariates: %s\n", covar_file))

# --- Nagelkerke R-squared ---
nagelkerke_r2 <- function(full_model, null_model, n) {
  ll_full <- as.numeric(logLik(full_model))
  ll_null <- as.numeric(logLik(null_model))
  cox_snell <- 1 - exp((2 / n) * (ll_null - ll_full))
  max_r2 <- 1 - exp((2 / n) * ll_null)
  cox_snell / max_r2
}

# --- Load data ---
cat("\n--- Loading data ---\n")
prs_all <- fread(file.path(results_dir, "prs/cluster_prs_all.tsv"))
pheno   <- fread(pheno_file)
covar   <- fread(covar_file)

cat(sprintf("  PRS: %d individuals, columns: %s\n",
            nrow(prs_all), paste(colnames(prs_all), collapse = ", ")))

# Detect cluster PRS columns
prs_cols <- grep("^PRS_K", colnames(prs_all), value = TRUE)
cat(sprintf("  Cluster PRS columns: %s\n", paste(prs_cols, collapse = ", ")))

# --- Merge phenotypes + covariates ---
merged <- merge(prs_all, pheno, by = c("FID", "IID"), all.x = TRUE)
merged <- merge(merged, covar, by = c("FID", "IID"), all.x = TRUE)
cat(sprintf("  Merged: %d individuals\n", nrow(merged)))

# --- Create synthetic outcome (T2D OR CAD) ---
merged[, Synthetic := fifelse(T2D == 1 & CAD == 1, 1L,
                      fifelse(T2D == 0 | CAD == 0, 0L, NA_integer_))]

cat(sprintf("  Synthetic cases: %d, controls: %d, NA: %d\n",
            sum(merged$Synthetic == 1, na.rm = TRUE),
            sum(merged$Synthetic == 0, na.rm = TRUE),
            sum(is.na(merged$Synthetic))))

# --- Define groups ---
groups <- list(
  eur_train      = merged[group == "eur_train"],
  eur_validation = merged[group == "eur_validation"],
  afr_validation = merged[group == "afr_validation"],
  eas_validation = merged[group == "eas_validation"],
  sas_validation = merged[group == "sas_validation"]
)

for (g in names(groups)) {
  cat(sprintf("  %s: %d individuals\n", g, nrow(groups[[g]])))
}

# --- Covariate terms ---
covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))

# --- Define models ---
outcomes <- c("T2D", "CAD", "Synthetic")

# --- Run regressions for one group ---
run_association <- function(data, group_name) {
  cat(sprintf("\n=== Running models on: %s (%d individuals) ===\n",
              group_name, nrow(data)))

  # Convert Batch to factor if present
  if ("Batch" %in% colnames(data)) {
    data[, Batch := as.factor(Batch)]
  }

  results_list <- list()

  for (outcome in outcomes) {
    if (!outcome %in% colnames(data)) {
      cat(sprintf("  SKIP %s: column not found\n", outcome))
      next
    }

    # --- Individual models: outcome ~ PRS_K{i} + covariates ---
    for (prs_col in prs_cols) {
      model_data <- data[, c(outcome, prs_col, covar_terms), with = FALSE]
      model_data <- na.omit(model_data)

      n_cases <- sum(model_data[[outcome]] == 1)
      n_controls <- sum(model_data[[outcome]] == 0)

      if (n_cases < 10 || n_controls < 10) {
        cat(sprintf("  SKIP %s ~ %s: too few cases (%d) or controls (%d)\n",
                    outcome, prs_col, n_cases, n_controls))
        next
      }

      # Standardize PRS
      model_data[[prs_col]] <- scale(model_data[[prs_col]])[, 1]

      rhs <- paste(c(prs_col, covar_terms), collapse = " + ")
      formula_full <- as.formula(paste(outcome, "~", rhs))
      formula_null <- as.formula(paste(outcome, "~", paste(covar_terms, collapse = " + ")))

      fit_full <- tryCatch(
        glm(formula_full, data = model_data, family = binomial),
        error = function(e) { cat(sprintf("  ERROR %s ~ %s: %s\n", outcome, prs_col, e$message)); NULL }
      )
      if (is.null(fit_full)) next

      fit_null <- glm(formula_null, data = model_data, family = binomial)

      coef_summ <- summary(fit_full)$coefficients
      prs_row <- coef_summ[prs_col, ]

      or <- exp(prs_row["Estimate"])
      ci_lower <- exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"])
      ci_upper <- exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"])
      p_val <- prs_row["Pr(>|z|)"]
      r2 <- nagelkerke_r2(fit_full, fit_null, nrow(model_data))

      # AUC
      pred_probs <- predict(fit_full, type = "response")
      auc_val <- tryCatch(
        as.numeric(auc(roc(model_data[[outcome]], pred_probs, quiet = TRUE))),
        error = function(e) NA_real_
      )

      cat(sprintf("  %s ~ %s: OR=%.3f [%.3f-%.3f], P=%.2e, R2=%.4f, AUC=%.3f (%d/%d)\n",
                  outcome, prs_col, or, ci_lower, ci_upper, p_val, r2, auc_val,
                  n_cases, n_controls))

      results_list[[length(results_list) + 1]] <- data.table(
        group = group_name, outcome = outcome, model_type = "individual",
        predictor = prs_col,
        OR = round(or, 4), CI_lower = round(ci_lower, 4), CI_upper = round(ci_upper, 4),
        p_value = p_val, nagelkerke_r2 = round(r2, 6), AUC = round(auc_val, 4),
        n_cases = n_cases, n_controls = n_controls, n_total = n_cases + n_controls
      )
    }

    # --- Joint model: outcome ~ PRS_K1 + ... + PRS_K6 + covariates ---
    model_data <- data[, c(outcome, prs_cols, covar_terms), with = FALSE]
    model_data <- na.omit(model_data)

    n_cases <- sum(model_data[[outcome]] == 1)
    n_controls <- sum(model_data[[outcome]] == 0)

    if (n_cases < 10 || n_controls < 10) {
      cat(sprintf("  SKIP %s ~ joint: too few cases/controls\n", outcome))
      next
    }

    # Standardize all PRS columns
    for (pc in prs_cols) {
      model_data[[pc]] <- scale(model_data[[pc]])[, 1]
    }

    rhs <- paste(c(prs_cols, covar_terms), collapse = " + ")
    formula_full <- as.formula(paste(outcome, "~", rhs))
    formula_null <- as.formula(paste(outcome, "~", paste(covar_terms, collapse = " + ")))

    fit_full <- tryCatch(
      glm(formula_full, data = model_data, family = binomial),
      error = function(e) { cat(sprintf("  ERROR %s ~ joint: %s\n", outcome, e$message)); NULL }
    )
    if (is.null(fit_full)) next

    fit_null <- glm(formula_null, data = model_data, family = binomial)

    coef_summ <- summary(fit_full)$coefficients
    r2 <- nagelkerke_r2(fit_full, fit_null, nrow(model_data))

    # AUC for joint model
    pred_probs <- predict(fit_full, type = "response")
    auc_val <- tryCatch(
      as.numeric(auc(roc(model_data[[outcome]], pred_probs, quiet = TRUE))),
      error = function(e) NA_real_
    )

    cat(sprintf("  %s ~ JOINT: R2=%.4f, AUC=%.3f (%d/%d)\n",
                outcome, r2, auc_val, n_cases, n_controls))

    # Extract each PRS coefficient from joint model
    for (pc in prs_cols) {
      if (pc %in% rownames(coef_summ)) {
        prs_row <- coef_summ[pc, ]
        or <- exp(prs_row["Estimate"])
        ci_lower <- exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"])
        ci_upper <- exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"])
        p_val <- prs_row["Pr(>|z|)"]
      } else {
        or <- ci_lower <- ci_upper <- p_val <- NA_real_
      }

      results_list[[length(results_list) + 1]] <- data.table(
        group = group_name, outcome = outcome, model_type = "joint",
        predictor = pc,
        OR = round(or, 4), CI_lower = round(ci_lower, 4), CI_upper = round(ci_upper, 4),
        p_value = p_val, nagelkerke_r2 = round(r2, 6), AUC = round(auc_val, 4),
        n_cases = n_cases, n_controls = n_controls, n_total = n_cases + n_controls
      )
    }
  }

  rbindlist(results_list)
}

# --- Run on all groups ---
all_results <- rbindlist(lapply(names(groups), function(g) {
  run_association(groups[[g]], g)
}))

# --- Write results ---
cat("\n--- Writing results ---\n")

for (g in unique(all_results$group)) {
  out_file <- file.path(results_dir, sprintf("association_results_%s.csv", g))
  fwrite(all_results[group == g], out_file)
  cat(sprintf("  %s: %d rows -> %s\n", g, nrow(all_results[group == g]), out_file))
}

# Combined file
combined_path <- file.path(results_dir, "association_results_all.csv")
fwrite(all_results, combined_path)
cat(sprintf("  Combined: %d rows -> %s\n", nrow(all_results), combined_path))

# --- Print summary table ---
cat("\n=== Results Summary ===\n")
print(all_results[, .(group, outcome, model_type, predictor,
                       OR, CI_lower, CI_upper, p_value, nagelkerke_r2, AUC,
                       n_cases, n_controls)],
      nrow = Inf)

cat("\n=== Done ===\n")
