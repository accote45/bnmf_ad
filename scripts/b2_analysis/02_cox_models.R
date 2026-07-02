#!/usr/bin/env Rscript
# 02_cox_models.R
# Cox proportional hazards models for time-to-comorbidity.
# Build on EUR training, evaluate on EUR validation.
#
# Three directions:
#   CAD-first -> T2D (event = subsequent T2D diagnosis)
#   T2D-first -> CAD (event = subsequent CAD diagnosis)
#   Combined  (event = developing comorbidity regardless of order)
#
# For each direction, fit individual models (one per cluster PRS).
# PRS is continuous, standardized on training set. HR per SD reported.
#
# Usage:
#   Rscript scripts/b2_analysis/02_cox_models.R --config config/b2_config.yaml
#   Rscript scripts/b2_analysis/02_cox_models.R --config config/b2_1_config.yaml

library(data.table)
library(yaml)
library(survival)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") {
    config_path <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}
if (is.null(config_path)) stop("--config is required. Usage: Rscript 02_cox_models.R --config <config.yaml>")

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
clusters <- cfg$analysis$clusters
prs_cols <- paste0("PRS_", clusters)

cat("=== B2 Step 2: Cox Proportional Hazards Models ===\n")
cat(sprintf("  Config: %s\n\n", config_path))

# --- Load data ---
cat("--- Loading data ---\n")

# Survival phenotypes
surv <- fread(file.path(results_dir, "survival_phenotypes.tsv"))
cat(sprintf("  Survival phenotypes: %d individuals\n", nrow(surv)))

# PRS scores
prs <- fread(cfg$b1_results$prs_file)
cat(sprintf("  PRS: %d individuals\n", nrow(prs)))

# Covariates
covar <- fread(cfg$phenotypes$covariate_file)
cat(sprintf("  Covariates: %d individuals\n", nrow(covar)))

# Sample splits
eur_train_ids <- fread(file.path(cfg$b1_results$samples_dir, "eur_train.keep"))$V1
eur_val_ids <- fread(file.path(cfg$b1_results$samples_dir, "eur_validation.keep"))$V1

# --- Merge all data ---
cat("\n--- Merging data ---\n")
dt <- merge(surv, prs[, c("IID", prs_cols), with = FALSE], by = "IID")
dt <- merge(dt, covar[, .(IID, sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)],
            by = "IID")
cat(sprintf("  Merged dataset: %d individuals\n", nrow(dt)))

# Compute age² at first dx
dt[, age_at_first_dx2 := age_at_first_dx^2]

# Split into train and validation
train <- dt[IID %in% eur_train_ids]
val   <- dt[IID %in% eur_val_ids]
cat(sprintf("  EUR train: %d\n", nrow(train)))
cat(sprintf("  EUR validation: %d\n", nrow(val)))

# --- Define covariate terms ---
covar_terms <- c("age_at_first_dx", "age_at_first_dx2", "sex",
                 paste0("PC", 1:10))

# --- Run Cox models ---
cat("\n--- Running Cox models ---\n")

directions <- c("CAD", "T2D", "Combined")  # first_dx value or "Combined"
dir_labels <- c("CAD_to_T2D", "T2D_to_CAD", "Combined")

results_list <- list()

for (d in seq_along(directions)) {
  dx <- directions[d]
  dir_label <- dir_labels[d]

  if (dx == "Combined") {
    cat("\n  Direction: Combined (any comorbidity)\n")
    train_sub <- copy(train)
    val_sub   <- copy(val)
  } else {
    cat(sprintf("\n  Direction: %s-first -> %s\n",
                dx, ifelse(dx == "CAD", "T2D", "CAD")))
    train_sub <- train[first_dx == dx]
    val_sub   <- val[first_dx == dx]
  }

  cat(sprintf("    Train: %d individuals, %d events (%.1f%%)\n",
              nrow(train_sub), sum(train_sub$event),
              100 * mean(train_sub$event)))
  cat(sprintf("    Validation: %d individuals, %d events (%.1f%%)\n",
              nrow(val_sub), sum(val_sub$event),
              100 * mean(val_sub$event)))

  for (k in seq_along(clusters)) {
    prs_col <- prs_cols[k]
    cluster_name <- clusters[k]

    # Standardize PRS on training set
    train_mean <- mean(train_sub[[prs_col]], na.rm = TRUE)
    train_sd   <- sd(train_sub[[prs_col]], na.rm = TRUE)

    train_sub[, prs_z := (get(prs_col) - train_mean) / train_sd]
    val_sub[, prs_z := (get(prs_col) - train_mean) / train_sd]

    # Build formula
    form <- as.formula(paste(
      "Surv(time_years, event) ~",
      paste(c("prs_z", covar_terms), collapse = " + ")
    ))

    # Fit on training
    fit_train <- tryCatch(
      coxph(form, data = train_sub),
      error = function(e) { cat(sprintf("    ERROR %s: %s\n", cluster_name, e$message)); NULL },
      warning = function(w) {
        cat(sprintf("    WARNING %s: %s\n", cluster_name, w$message))
        suppressWarnings(coxph(form, data = train_sub))
      }
    )

    if (is.null(fit_train)) next

    # Extract results from training fit
    s <- summary(fit_train)
    hr <- exp(coef(fit_train)["prs_z"])
    ci <- exp(confint(fit_train)["prs_z", ])
    pval <- s$coefficients["prs_z", "Pr(>|z|)"]
    conc_train <- s$concordance["C"]

    # PH assumption test
    ph_test <- tryCatch({
      zph <- cox.zph(fit_train)
      zph$table["prs_z", "p"]
    }, error = function(e) NA_real_)

    # Evaluate on validation: predict linear predictor and compute concordance
    # Negate lp because concordance() expects higher predictor = longer survival,
    # but coxph lp has higher = higher risk (shorter survival)
    conc_val <- tryCatch({
      lp_val <- predict(fit_train, newdata = val_sub, type = "lp")
      surv_val <- Surv(val_sub$time_years, val_sub$event)
      1 - concordance(surv_val ~ lp_val)$concordance
    }, error = function(e) NA_real_)

    results_list[[length(results_list) + 1]] <- data.table(
      direction = dir_label,
      cluster = cluster_name,
      HR = round(hr, 4),
      CI_lower = round(ci[1], 4),
      CI_upper = round(ci[2], 4),
      p_value = pval,
      c_index_train = round(conc_train, 4),
      c_index_val = round(conc_val, 4),
      ph_test_p = round(ph_test, 4),
      n_events_train = sum(train_sub$event),
      n_total_train = nrow(train_sub),
      n_events_val = sum(val_sub$event),
      n_total_val = nrow(val_sub)
    )

    cat(sprintf("    %s: HR=%.3f (%.3f-%.3f), p=%.2e, C_train=%.3f, C_val=%.3f, PH_p=%.3f\n",
                cluster_name, hr, ci[1], ci[2], pval,
                conc_train, ifelse(is.na(conc_val), 0, conc_val),
                ifelse(is.na(ph_test), 0, ph_test)))

    # Clean up temp column
    train_sub[, prs_z := NULL]
    val_sub[, prs_z := NULL]
  }
}

# --- Combine and save results ---
results_dt <- rbindlist(results_list)

out_file <- file.path(results_dir, "cox_results.csv")
fwrite(results_dt, out_file)
cat(sprintf("\n  Saved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Print formatted summary ---
cat("\n--- Summary Table ---\n")
cat(sprintf("%-15s %-8s %8s %15s %12s %10s %10s\n",
            "Direction", "Cluster", "HR", "95% CI", "p-value", "C_train", "C_val"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (r in seq_len(nrow(results_dt))) {
  row <- results_dt[r]
  cat(sprintf("%-15s %-8s %8.3f (%6.3f-%6.3f) %12.2e %10.3f %10.3f\n",
              row$direction, row$cluster, row$HR,
              row$CI_lower, row$CI_upper, row$p_value,
              row$c_index_train, ifelse(is.na(row$c_index_val), 0, row$c_index_val)))
}

cat("\n=== Step 2 complete ===\n")
