#!/usr/bin/env Rscript
# 05_quantile_prs_analysis.R
# Quantile-based cluster PRS associations: OR for top 20%, 10%, 5%, 1%
# vs the rest. Writes quantile_prs_associations.csv.
#
# Forest plots are drawn separately by 05b_quantile_forest_plots.R (which reads
# the CSV written here), so figures can be tweaked without rerunning regressions.
#
# Usage: Rscript scripts/b1_analysis/05_quantile_prs_analysis.R

library(data.table)
library(yaml)

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

cat("=== B1 Step 5: Quantile PRS Analysis ===\n")

# ============================================================
# Load and merge data (same as 03_b1_forest_plots.R)
# ============================================================
cat("\n--- Loading data ---\n")

prs_all <- fread(file.path(results_dir, "prs/cluster_prs_all.tsv"))
prs_cols <- grep("^PRS_K", colnames(prs_all), value = TRUE)
n_clusters <- length(prs_cols)
cat(sprintf("  PRS: %d individuals, %d clusters\n", nrow(prs_all), n_clusters))

# Genome-wide PRS (bNMF-matched; built by 01b_compute_genomewide_prs.R).
# Carries the sample `group`, which prs_all already provides — drop the dup.
gw_all <- fread(file.path(results_dir, "prs/genomewide_prs_all.tsv"))
gw_all[, group := NULL]

pheno <- fread(cfg$phenotypes$phenotype_file)
covar <- fread(cfg$phenotypes$covariate_file)

dat <- merge(prs_all, gw_all, by = c("FID", "IID"))
dat <- merge(dat, pheno, by = c("FID", "IID"))
dat <- merge(dat, covar, by = c("FID", "IID"))
dat[, Synthetic := fifelse(T2D == 1 & CAD == 1, 1L,
                   fifelse(T2D == 0 | CAD == 0, 0L, NA_integer_))]
cat(sprintf("  Merged: %d individuals\n", nrow(dat)))

covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))

# ============================================================
# Quantile associations
# ============================================================
cat("\n--- Computing quantile associations ---\n")

quantile_fracs <- c(0.20, 0.10, 0.05, 0.01)
quantile_names <- c("top_20pct", "top_10pct", "top_5pct", "top_1pct")
outcomes <- c("T2D", "CAD", "Synthetic")
val_groups <- c("eur_validation", "afr_validation", "eas_validation", "sas_validation")

# All PRS columns to test: cluster PRS + genome-wide (trait-specific + the
# supplementary combined-variant scores over the shared bNMF universe).
all_prs <- c(prs_cols, "GW_T2D", "GW_CAD", "GW_T2D_combined", "GW_CAD_combined")
all_prs_labels <- c(gsub("PRS_", "", prs_cols), "GW T2D", "GW CAD",
                    "GW T2D (combined variants)", "GW CAD (combined variants)")

results_list <- list()

for (grp in val_groups) {
  grp_data <- dat[group == grp]
  cat(sprintf("\n  %s (n=%d)\n", grp, nrow(grp_data)))

  for (pi in seq_along(all_prs)) {
    prs_col <- all_prs[pi]
    prs_label <- all_prs_labels[pi]

    for (qi in seq_along(quantile_fracs)) {
      frac <- quantile_fracs[qi]
      qname <- quantile_names[qi]

      # Compute threshold within this group
      thresh_val <- quantile(grp_data[[prs_col]], probs = 1 - frac, na.rm = TRUE)
      grp_data[, top_ind := fifelse(get(prs_col) >= thresh_val, 1L, 0L)]
      n_top <- sum(grp_data$top_ind == 1, na.rm = TRUE)

      for (outcome in outcomes) {
        model_data <- grp_data[, c(outcome, "top_ind", covar_terms), with = FALSE]
        model_data <- na.omit(model_data)

        n_cases <- sum(model_data[[outcome]] == 1)
        n_controls <- sum(model_data[[outcome]] == 0)
        if (n_cases < 10 || n_controls < 10) next

        n_cases_top <- sum(model_data[[outcome]] == 1 & model_data$top_ind == 1)
        n_controls_top <- sum(model_data[[outcome]] == 0 & model_data$top_ind == 1)

        fit <- tryCatch(
          glm(as.formula(paste(outcome, "~ top_ind +",
              paste(covar_terms, collapse = " + "))),
              data = model_data, family = binomial),
          error = function(e) NULL
        )
        if (is.null(fit)) next

        coef_summ <- summary(fit)$coefficients
        if (!"top_ind" %in% rownames(coef_summ)) next

        est <- coef_summ["top_ind", "Estimate"]
        se <- coef_summ["top_ind", "Std. Error"]
        p_val <- coef_summ["top_ind", "Pr(>|z|)"]

        results_list[[length(results_list) + 1]] <- data.table(
          group = grp,
          cluster = prs_label,
          outcome = outcome,
          quantile = qname,
          OR = round(exp(est), 4),
          CI_lower = round(exp(est - 1.96 * se), 4),
          CI_upper = round(exp(est + 1.96 * se), 4),
          p_value = p_val,
          n_top = n_top,
          n_cases_top = n_cases_top,
          n_controls_top = n_controls_top
        )
      }
    }
  }
}

grp_data[, top_ind := NULL]

results <- rbindlist(results_list)

# Add cluster labels
cluster_labels <- if (!is.null(cfg$cluster_labels)) unlist(cfg$cluster_labels) else setNames(paste0("K", 1:n_clusters), paste0("K", 1:n_clusters))
label_map <- c(cluster_labels, "GW T2D" = "GW T2D", "GW CAD" = "GW CAD",
               "GW T2D (combined variants)" = "GW T2D (combined variants)",
               "GW CAD (combined variants)" = "GW CAD (combined variants)")
results[, cluster_label := label_map[cluster]]

fwrite(results, file.path(results_dir, "quantile_prs_associations.csv"))
cat(sprintf("\n--- Saved: quantile_prs_associations.csv (%d rows) ---\n", nrow(results)))

# Print EUR validation top 10% summary
cat("\n  EUR validation top 10%:\n")
print(results[group == "eur_validation" & quantile == "top_10pct",
              .(cluster_label, outcome, OR, CI_lower, CI_upper, p_value)])

cat("\n=== Done (regressions) ===\n")
cat("  Forest plots: Rscript scripts/b1_analysis/05b_quantile_forest_plots.R\n")
