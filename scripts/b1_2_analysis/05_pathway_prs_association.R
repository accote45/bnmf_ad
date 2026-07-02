#!/usr/bin/env Rscript
# 05_pathway_prs_association.R
# Compute OR/SD for each pathway PRS against T2D, CAD, and Synthetic outcomes.
# Train on EUR train, evaluate on EUR validation (matching b1 approach).
#
# Usage:
#   Rscript scripts/b1_2_analysis/05_pathway_prs_association.R
#   Rscript scripts/b1_2_analysis/05_pathway_prs_association.R --config config/b1_2_config.yaml

library(data.table)
library(yaml)
library(pROC)
library(parallel)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_2_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
prset_dir <- file.path(results_dir, "prset")
n_cores <- cfg$prset$threads
n_top <- cfg$analysis$n_top_pathways

cat("=== B1.2 Step 5: Pathway PRS Association Testing ===\n")
cat(sprintf("  Config: %s\n", config_path))
cat(sprintf("  Cores: %d\n", n_cores))
cat(sprintf("  Top pathways: %d\n", n_top))

# --- Nagelkerke R-squared ---
nagelkerke_r2 <- function(full_model, null_model, n) {
  ll_full <- as.numeric(logLik(full_model))
  ll_null <- as.numeric(logLik(null_model))
  cox_snell <- 1 - exp((2 / n) * (ll_null - ll_full))
  max_r2 <- 1 - exp((2 / n) * ll_null)
  cox_snell / max_r2
}

# --- Select top pathways from summary ---
cat("\n--- Selecting top pathways by PRS.R2 ---\n")
summary_file <- file.path(prset_dir, "prset_output.summary")
pw_summary <- fread(summary_file)
max_snps <- cfg$analysis$max_pathway_snps
pw_summary <- pw_summary[Set != "Base"]
if (!is.null(max_snps)) {
  pw_summary <- pw_summary[Num_SNP <= max_snps]
  cat(sprintf("  Filtered to pathways with <= %d SNPs: %d remain\n", max_snps, nrow(pw_summary)))
}
pw_summary <- pw_summary[order(-PRS.R2)]
top_pathways <- pw_summary$Set[1:n_top]

cat(sprintf("  Top %d pathways by R²:\n", n_top))
for (i in seq_along(top_pathways)) {
  row <- pw_summary[Set == top_pathways[i]]
  cat(sprintf("    %d. %s (R²=%.6f, P=%.2e, %d SNPs)\n",
              i, top_pathways[i], row$PRS.R2, row$P, row$Num_SNP))
}

# --- Load PRSet .best file (top pathways only) ---
cat("\n--- Loading PRSet .best file (selective columns) ---\n")
best_file <- file.path(prset_dir, "prset_output.best")
prs <- fread(best_file, select = c("FID", "IID", top_pathways), nThread = 2)
cat(sprintf("  PRSet .best: %d individuals, %d columns\n", nrow(prs), ncol(prs)))

pathway_cols <- top_pathways
cat(sprintf("  Pathway PRS columns: %d\n", length(pathway_cols)))

# --- Load phenotypes and covariates ---
cat("\n--- Loading phenotypes and covariates ---\n")
pheno <- fread(cfg$phenotypes$phenotype_file)
covar <- fread(cfg$phenotypes$covariate_file)

# Create Synthetic
pheno[, Synthetic := fifelse(T2D == 1 | CAD == 1, 1L,
                    fifelse(T2D == 0 & CAD == 0, 0L, NA_integer_))]

# --- Merge all data ---
dat <- merge(prs, pheno[, .(FID, IID, T2D, CAD, Synthetic)], by = c("FID", "IID"))
dat <- merge(dat, covar, by = c("FID", "IID"))
cat(sprintf("  Merged: %d individuals\n", nrow(dat)))

# --- Split into train and validation ---
cat("\n--- Splitting into train/validation ---\n")
train_ids <- fread(cfg$samples$eur_train, header = FALSE)$V1
val_ids <- fread(cfg$samples$eur_validation, header = FALSE)$V1

train <- dat[IID %in% train_ids]
val <- dat[IID %in% val_ids]
cat(sprintf("  Train: %d individuals\n", nrow(train)))
cat(sprintf("  Validation: %d individuals\n", nrow(val)))

# Convert Batch to factor
if ("Batch" %in% colnames(train)) {
  train[, Batch := as.factor(Batch)]
  val[, Batch := as.factor(Batch)]
}

# --- Define covariate terms ---
covar_terms <- c("age", "age2", "sex", paste0("PC", 1:10))

# --- Define outcomes ---
outcomes <- c("T2D", "CAD", "Synthetic")

# --- Run association tests ---
cat("\n--- Running association tests ---\n")
cat(sprintf("  %d pathways x %d outcomes = %d models\n",
            length(pathway_cols), length(outcomes),
            length(pathway_cols) * length(outcomes)))

run_one_pathway_outcome <- function(pw_col, outcome, train_data, val_data) {
  # Extract needed columns
  train_sub <- train_data[, c(outcome, pw_col, covar_terms), with = FALSE]
  train_sub <- na.omit(train_sub)

  val_sub <- val_data[, c(outcome, pw_col, covar_terms), with = FALSE]
  val_sub <- na.omit(val_sub)

  n_cases_train <- sum(train_sub[[outcome]] == 1)
  n_controls_train <- sum(train_sub[[outcome]] == 0)
  n_cases_val <- sum(val_sub[[outcome]] == 1)
  n_controls_val <- sum(val_sub[[outcome]] == 0)

  if (n_cases_train < 10 || n_controls_train < 10 ||
      n_cases_val < 10 || n_controls_val < 10) {
    return(NULL)
  }

  # Check variance
  if (sd(train_sub[[pw_col]], na.rm = TRUE) == 0) return(NULL)

  # Standardize PRS on training set
  train_mean <- mean(train_sub[[pw_col]], na.rm = TRUE)
  train_sd <- sd(train_sub[[pw_col]], na.rm = TRUE)
  if (train_sd == 0) return(NULL)

  train_sub[[pw_col]] <- (train_sub[[pw_col]] - train_mean) / train_sd
  val_sub[[pw_col]] <- (val_sub[[pw_col]] - train_mean) / train_sd

  rhs <- paste(c(pw_col, covar_terms), collapse = " + ")
  formula_full <- as.formula(paste(outcome, "~", rhs))
  formula_null <- as.formula(paste(outcome, "~", paste(covar_terms, collapse = " + ")))

  # Fit on training
  fit_train <- tryCatch(
    glm(formula_full, data = train_sub, family = binomial),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(glm(formula_full, data = train_sub, family = binomial))
  )
  if (is.null(fit_train)) return(NULL)

  # Extract training OR
  coef_summ <- summary(fit_train)$coefficients
  if (!pw_col %in% rownames(coef_summ)) return(NULL)

  prs_row <- coef_summ[pw_col, ]
  or_train <- exp(prs_row["Estimate"])
  ci_lower_train <- exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"])
  ci_upper_train <- exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"])
  p_val_train <- prs_row["Pr(>|z|)"]

  fit_null_train <- tryCatch(
    glm(formula_null, data = train_sub, family = binomial),
    error = function(e) NULL
  )
  r2_train <- if (!is.null(fit_null_train)) nagelkerke_r2(fit_train, fit_null_train, nrow(train_sub)) else NA_real_

  # Evaluate on validation
  fit_val <- tryCatch(
    glm(formula_full, data = val_sub, family = binomial),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(glm(formula_full, data = val_sub, family = binomial))
  )

  if (!is.null(fit_val)) {
    coef_val <- summary(fit_val)$coefficients
    if (pw_col %in% rownames(coef_val)) {
      prs_row_val <- coef_val[pw_col, ]
      or_val <- exp(prs_row_val["Estimate"])
      ci_lower_val <- exp(prs_row_val["Estimate"] - 1.96 * prs_row_val["Std. Error"])
      ci_upper_val <- exp(prs_row_val["Estimate"] + 1.96 * prs_row_val["Std. Error"])
      p_val_val <- prs_row_val["Pr(>|z|)"]
    } else {
      or_val <- ci_lower_val <- ci_upper_val <- p_val_val <- NA_real_
    }

    fit_null_val <- tryCatch(
      glm(formula_null, data = val_sub, family = binomial),
      error = function(e) NULL
    )
    r2_val <- if (!is.null(fit_null_val)) nagelkerke_r2(fit_val, fit_null_val, nrow(val_sub)) else NA_real_

    auc_val <- tryCatch({
      pred_probs <- predict(fit_val, type = "response")
      as.numeric(auc(roc(val_sub[[outcome]], pred_probs, quiet = TRUE)))
    }, error = function(e) NA_real_)
  } else {
    or_val <- ci_lower_val <- ci_upper_val <- p_val_val <- r2_val <- auc_val <- NA_real_
  }

  data.table(
    pathway = pw_col,
    outcome = outcome,
    OR_train = round(or_train, 4),
    CI_lower_train = round(ci_lower_train, 4),
    CI_upper_train = round(ci_upper_train, 4),
    p_value_train = p_val_train,
    r2_train = round(r2_train, 6),
    OR_val = round(or_val, 4),
    CI_lower_val = round(ci_lower_val, 4),
    CI_upper_val = round(ci_upper_val, 4),
    p_value_val = p_val_val,
    r2_val = round(r2_val, 6),
    AUC_val = round(auc_val, 4),
    n_cases_train = n_cases_train,
    n_controls_train = n_controls_train,
    n_cases_val = n_cases_val,
    n_controls_val = n_controls_val
  )
}

# Build task list
tasks <- expand.grid(pathway = pathway_cols, outcome = outcomes,
                     stringsAsFactors = FALSE)
cat(sprintf("  Total tasks: %d\n", nrow(tasks)))

# Run in parallel
cat("  Running regressions...\n")
t0 <- Sys.time()

results <- mclapply(seq_len(nrow(tasks)), function(i) {
  run_one_pathway_outcome(
    tasks$pathway[i], tasks$outcome[i], train, val
  )
}, mc.cores = min(n_cores, 4))

# Filter NULLs and combine
results <- results[!sapply(results, is.null)]
results_dt <- rbindlist(results)

t1 <- Sys.time()
cat(sprintf("  Completed %d models in %.1f minutes\n",
            nrow(results_dt), as.numeric(difftime(t1, t0, units = "mins"))))

# --- FDR correction ---
cat("\n--- Applying FDR correction ---\n")
for (outcome_name in outcomes) {
  idx <- results_dt$outcome == outcome_name
  results_dt[idx, fdr_val := p.adjust(p_value_val, method = "BH")]
  results_dt[idx, fdr_train := p.adjust(p_value_train, method = "BH")]
}

# --- Write results ---
out_file <- file.path(results_dir, "pathway_association_results.csv")
fwrite(results_dt, out_file)
cat(sprintf("\n  Saved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Print top pathways by Synthetic validation OR ---
cat(sprintf("\n--- Top %d pathways by Synthetic validation OR ---\n", n_top))
synth_results <- results_dt[outcome == "Synthetic"][order(-OR_val)]
if (nrow(synth_results) > 0) {
  print(synth_results[1:min(n_top, nrow(synth_results)),
    .(pathway, OR_val, CI_lower_val, CI_upper_val, p_value_val, fdr_val)])
}

cat("\n=== Step 5 complete ===\n")
