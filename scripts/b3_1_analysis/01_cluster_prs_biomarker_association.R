#!/usr/bin/env Rscript
# 01_cluster_prs_biomarker_association.R
# Test associations between bNMF cluster PRS and continuous biomarker traits.
# Runs linear regressions (individual + joint) on EUR validation cohort.
#
# Usage:
#   Rscript scripts/b3_1_analysis/01_cluster_prs_biomarker_association.R
#   Rscript scripts/b3_1_analysis/01_cluster_prs_biomarker_association.R --config config/b3_1_config.yaml

library(data.table)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b3_1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

hard_assignment <- "--hard-assignment" %in% args

cfg <- read_yaml(config_path)
phenotype_dir <- cfg$results_dir
results_dir <- cfg$results_dir
primary_group <- cfg$analysis$primary_group

if (hard_assignment) {
  cfg$b1_results$prs_file <- "results/b1_analysis/prs_hard_assignment/prs/cluster_prs_all.tsv"
  results_dir <- "results/b3_1_analysis/prs_hard_assignment"
  cat("  ** HARD ASSIGNMENT MODE **\n")
  cat(sprintf("  PRS file: %s\n", cfg$b1_results$prs_file))
  cat(sprintf("  Output dir: %s\n\n", results_dir))
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B3 Step 1: Cluster PRS Biomarker Association ===\n")
cat(sprintf("  Config:  %s\n", config_path))
cat(sprintf("  Group:   %s\n", primary_group))

# --- Load data ---
cat("\n--- Loading data ---\n")

prs_all <- fread(cfg$b1_results$prs_file)
biomarkers <- fread(file.path(phenotype_dir, "biomarker_phenotypes.tsv"))
covar <- fread(cfg$phenotypes$covariate_file)

cat(sprintf("  PRS: %d individuals\n", nrow(prs_all)))
cat(sprintf("  Biomarkers: %d individuals\n", nrow(biomarkers)))
cat(sprintf("  Covariates: %d individuals\n", nrow(covar)))

# Detect cluster PRS columns
prs_cols <- grep("^PRS_K", colnames(prs_all), value = TRUE)
cat(sprintf("  Cluster PRS columns: %s\n", paste(prs_cols, collapse = ", ")))

# Drop age/sex from biomarkers (already in covariate file) to avoid merge collision
biomarkers[, c("age", "sex") := NULL]

# --- Merge ---
merged <- merge(prs_all, covar, by = c("FID", "IID"), all.x = TRUE)
merged <- merge(merged, biomarkers, by = c("FID", "IID"), all.x = TRUE)
cat(sprintf("  Merged: %d individuals\n", nrow(merged)))

# Filter to primary group
dat <- merged[group == primary_group]
cat(sprintf("  %s: %d individuals\n", primary_group, nrow(dat)))

# --- Define traits and covariates ---
outcome_traits <- c("Lpa", "CRP", "HbA1c_corrected", "glucose_corrected",
                    "eGFR", "non_HDL_chol", "triglycerides_corrected",
                    "LDL_corrected", "WHtR", "waist_circumference",
                    "BMI", "SBP_corrected", "PREVENT_ASCVD")

trait_display_map <- c(
  Lpa = "Lpa",
  CRP = "CRP",
  HbA1c_corrected = "HbA1c",
  glucose_corrected = "Fasting Glucose",
  eGFR = "eGFR",
  non_HDL_chol = "non-HDL Cholesterol",
  triglycerides_corrected = "Triglycerides",
  LDL_corrected = "LDL",
  WHtR = "WHtR",
  waist_circumference = "Waist Circumference",
  BMI = "BMI",
  SBP_corrected = "SBP",
  PREVENT_ASCVD = "PREVENT-ASCVD 10yr"
)

# Traits to log-transform (from config)
log_traits <- c("Lpa", "CRP", "triglycerides_corrected", "glucose_corrected")

covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))

# Ensure age2 exists
if (!"age2" %in% colnames(dat)) {
  dat[, age2 := as.numeric(age)^2]
}

# Convert Batch to factor for use as covariate
dat[, Batch := as.factor(Batch)]

# Pre-compute log-transformed columns
for (trait in outcome_traits) {
  if (trait %in% log_traits) {
    log_col <- paste0("log_", trait)
    dat[, (log_col) := log(get(trait))]
    dat[is.infinite(get(log_col)), (log_col) := NA]
  }
}

# --- Partial R-squared for linear models ---
partial_r2 <- function(full_model, null_model) {
  ss_res_full <- sum(residuals(full_model)^2)
  ss_res_null <- sum(residuals(null_model)^2)
  (ss_res_null - ss_res_full) / ss_res_null
}

# --- Sex-stratified analysis groups ---
sex_groups <- list(All = NULL, Male = 1L, Female = 0L)
covar_terms_all <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))
covar_terms_sex <- c("age", "age2", "Batch", paste0("PC", 1:10))

all_individual <- list()
all_joint <- list()

for (sg_name in names(sex_groups)) {
  sg_val <- sex_groups[[sg_name]]

  if (is.null(sg_val)) {
    sg_dat <- dat
    sg_covar <- covar_terms_all
  } else {
    sg_dat <- dat[sex == sg_val]
    sg_covar <- covar_terms_sex
  }

  cat(sprintf("\n=== Sex group: %s (N=%d) ===\n", sg_name, nrow(sg_dat)))

  # ============================================================
  # Individual models: trait ~ PRS_Ki + covariates
  # ============================================================
  cat("\n--- Running individual models ---\n")

  results_list <- list()

  for (trait in outcome_traits) {
    trait_col <- if (trait %in% log_traits) paste0("log_", trait) else trait

    for (prs_col in prs_cols) {
      model_data <- sg_dat[, c(trait_col, prs_col, sg_covar), with = FALSE]
      model_data <- na.omit(model_data)

      if (nrow(model_data) < 50) {
        cat(sprintf("  SKIP %s ~ %s: N=%d < 50\n", trait, prs_col, nrow(model_data)))
        next
      }

      model_data[[prs_col]] <- scale(model_data[[prs_col]])[, 1]

      rhs_full <- paste(c(prs_col, sg_covar), collapse = " + ")
      rhs_null <- paste(sg_covar, collapse = " + ")
      formula_full <- as.formula(paste(trait_col, "~", rhs_full))
      formula_null <- as.formula(paste(trait_col, "~", rhs_null))

      fit_full <- tryCatch(
        lm(formula_full, data = model_data),
        error = function(e) { cat(sprintf("  ERROR %s ~ %s: %s\n", trait, prs_col, e$message)); NULL }
      )
      if (is.null(fit_full)) next

      fit_null <- lm(formula_null, data = model_data)

      coef_summ <- summary(fit_full)$coefficients
      if (!prs_col %in% rownames(coef_summ)) next

      prs_row <- coef_summ[prs_col, ]
      beta <- prs_row["Estimate"]
      se <- prs_row["Std. Error"]
      ci_lower <- beta - 1.96 * se
      ci_upper <- beta + 1.96 * se
      p_val <- prs_row["Pr(>|t|)"]
      pr2 <- partial_r2(fit_full, fit_null)

      cat(sprintf("  [%s] %s ~ %s: beta=%.4f [%.4f, %.4f], P=%.2e (N=%d)\n",
                  sg_name, trait, prs_col, beta, ci_lower, ci_upper, p_val, nrow(model_data)))

      results_list[[length(results_list) + 1]] <- data.table(
        group = primary_group,
        sex_group = sg_name,
        trait = trait,
        trait_display = trait_display_map[trait],
        model_type = "individual",
        predictor = prs_col,
        beta = round(beta, 6),
        se = round(se, 6),
        ci_lower = round(ci_lower, 6),
        ci_upper = round(ci_upper, 6),
        p_value = p_val,
        partial_r2 = round(pr2, 8),
        n_obs = nrow(model_data)
      )
    }
  }

  sg_individual <- rbindlist(results_list)
  cat(sprintf("\n  Individual models (%s): %d results\n", sg_name, nrow(sg_individual)))

  # ============================================================
  # Joint models: trait ~ PRS_K1 + ... + PRS_K9 + covariates
  # ============================================================
  cat("\n--- Running joint models ---\n")

  joint_list <- list()

  for (trait in outcome_traits) {
    trait_col <- if (trait %in% log_traits) paste0("log_", trait) else trait

    model_data <- sg_dat[, c(trait_col, prs_cols, sg_covar), with = FALSE]
    model_data <- na.omit(model_data)

    if (nrow(model_data) < 50) {
      cat(sprintf("  SKIP %s ~ joint: N=%d < 50\n", trait, nrow(model_data)))
      next
    }

    for (pc in prs_cols) {
      model_data[[pc]] <- scale(model_data[[pc]])[, 1]
    }

    rhs_full <- paste(c(prs_cols, sg_covar), collapse = " + ")
    rhs_null <- paste(sg_covar, collapse = " + ")
    formula_full <- as.formula(paste(trait_col, "~", rhs_full))
    formula_null <- as.formula(paste(trait_col, "~", rhs_null))

    fit_full <- tryCatch(
      lm(formula_full, data = model_data),
      error = function(e) { cat(sprintf("  ERROR %s ~ joint: %s\n", trait, e$message)); NULL }
    )
    if (is.null(fit_full)) next

    fit_null <- lm(formula_null, data = model_data)
    pr2_overall <- partial_r2(fit_full, fit_null)

    coef_summ <- summary(fit_full)$coefficients

    cat(sprintf("  [%s] %s ~ JOINT: overall R2=%.6f (N=%d)\n",
                sg_name, trait, pr2_overall, nrow(model_data)))

    for (pc in prs_cols) {
      if (pc %in% rownames(coef_summ)) {
        prs_row <- coef_summ[pc, ]
        beta <- prs_row["Estimate"]
        se <- prs_row["Std. Error"]
        ci_lower <- beta - 1.96 * se
        ci_upper <- beta + 1.96 * se
        p_val <- prs_row["Pr(>|t|)"]
      } else {
        beta <- se <- ci_lower <- ci_upper <- p_val <- NA_real_
      }

      joint_list[[length(joint_list) + 1]] <- data.table(
        group = primary_group,
        sex_group = sg_name,
        trait = trait,
        trait_display = trait_display_map[trait],
        model_type = "joint",
        predictor = pc,
        beta = round(beta, 6),
        se = round(se, 6),
        ci_lower = round(ci_lower, 6),
        ci_upper = round(ci_upper, 6),
        p_value = p_val,
        partial_r2 = round(pr2_overall, 8),
        n_obs = nrow(model_data)
      )
    }
  }

  sg_joint <- rbindlist(joint_list)
  cat(sprintf("  Joint models (%s): %d results\n", sg_name, nrow(sg_joint)))

  all_individual[[sg_name]] <- sg_individual
  all_joint[[sg_name]] <- sg_joint
}

individual_results <- rbindlist(all_individual)
joint_results <- rbindlist(all_joint)

# ============================================================
# Multiple testing correction (per sex_group)
# ============================================================
cat("\n--- Multiple testing correction ---\n")

bonf_thresh <- 0.05 / cfg$analysis$bonferroni_n_tests

individual_results[, fdr_q := p.adjust(p_value, method = "BH"), by = sex_group]
individual_results[, bonferroni_sig := p_value < bonf_thresh]

for (sg in names(sex_groups)) {
  sg_res <- individual_results[sex_group == sg]
  n_fdr <- sum(sg_res$fdr_q < cfg$analysis$fdr_threshold, na.rm = TRUE)
  n_bonf <- sum(sg_res$bonferroni_sig, na.rm = TRUE)
  cat(sprintf("  %s: FDR < %.2f: %d / %d | Bonferroni: %d / %d\n",
              sg, cfg$analysis$fdr_threshold, n_fdr, nrow(sg_res), n_bonf, nrow(sg_res)))
}

joint_results[, fdr_q := p.adjust(p_value, method = "BH"), by = sex_group]
joint_results[, bonferroni_sig := p_value < bonf_thresh]

# ============================================================
# Write results
# ============================================================
cat("\n--- Writing results ---\n")

indiv_file <- file.path(results_dir, "biomarker_association_individual.csv")
fwrite(individual_results, indiv_file)
cat(sprintf("  Individual: %s (%d rows)\n", indiv_file, nrow(individual_results)))

joint_file <- file.path(results_dir, "biomarker_association_joint.csv")
fwrite(joint_results, joint_file)
cat(sprintf("  Joint: %s (%d rows)\n", joint_file, nrow(joint_results)))

all_results <- rbindlist(list(individual_results, joint_results), fill = TRUE)
combined_file <- file.path(results_dir, "biomarker_association_all.csv")
fwrite(all_results, combined_file)
cat(sprintf("  Combined: %s (%d rows)\n", combined_file, nrow(all_results)))

# --- Print summary ---
cat("\n=== FDR-significant results (individual models, All) ===\n")
sig_results <- individual_results[sex_group == "All" & fdr_q < cfg$analysis$fdr_threshold]
if (nrow(sig_results) > 0) {
  sig_results <- sig_results[order(fdr_q)]
  print(sig_results[, .(trait_display, predictor, beta, ci_lower, ci_upper,
                         p_value, fdr_q, partial_r2, n_obs)],
        nrow = Inf)
} else {
  cat("  No FDR-significant results.\n")
}

cat("\n=== Done ===\n")
