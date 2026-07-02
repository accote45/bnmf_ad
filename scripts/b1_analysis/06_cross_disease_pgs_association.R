#!/usr/bin/env Rscript
# 06_cross_disease_pgs_association.R
# Test cluster PGS associations with cross-disease cardiometabolic risk:
#   - T2D among CAD cases
#   - CAD among T2D cases
# Runs continuous (per-SD) and top 10% (binary indicator) models.
# Generates CSV results + supplementary Excel table.
#
# Usage:
#   Rscript scripts/b1_analysis/06_cross_disease_pgs_association.R --hard-assignment

library(data.table)
library(yaml)
library(openxlsx)

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
}

cluster_labels <- unlist(cfg$cluster_labels)
cluster_keys <- paste0("K", 1:10)

pheno_file <- cfg$phenotypes$phenotype_file
covar_file <- cfg$phenotypes$covariate_file

cat("=== B1 Step 6: Cross-Disease PGS Association ===\n")
cat(sprintf("  Results dir: %s\n", results_dir))
cat(sprintf("  Phenotypes:  %s\n", pheno_file))
cat(sprintf("  Covariates:  %s\n", covar_file))

# --- Load and merge data ---
cat("\n--- Loading data ---\n")

prs_all <- fread(file.path(results_dir, "prs/cluster_prs_all.tsv"))
pheno   <- fread(pheno_file)
covar   <- fread(covar_file)

prs_cols <- grep("^PRS_K", colnames(prs_all), value = TRUE)
cat(sprintf("  PRS: %d individuals, %d clusters\n", nrow(prs_all), length(prs_cols)))

merged <- merge(prs_all, pheno, by = c("FID", "IID"), all.x = TRUE)
merged <- merge(merged, covar, by = c("FID", "IID"), all.x = TRUE)
cat(sprintf("  Merged: %d individuals\n", nrow(merged)))

# --- Analysis parameters ---
covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))

cross_disease_tests <- list(
  list(restrict_col = "CAD", restrict_val = 1L,
       outcome = "T2D", label = "T2D among CAD cases"),
  list(restrict_col = "T2D", restrict_val = 1L,
       outcome = "CAD", label = "CAD among T2D cases")
)

val_groups <- c("eur_validation", "afr_validation",
                "eas_validation", "sas_validation")

# --- Formatting helpers ---
format_or <- function(x) sprintf("%.2f", x)

format_ci <- function(lower, upper) sprintf("%.2f–%.2f", lower, upper)

format_pval <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, sprintf("%.2e", p),
      ifelse(p < 0.01, sprintf("%.3f", p),
        sprintf("%.2f", p))))
}

# --- Run cross-disease associations ---
run_cross_disease <- function(data, group_name) {
  cat(sprintf("\n=== %s (%d individuals) ===\n", group_name, nrow(data)))

  data[, Batch := droplevels(as.factor(Batch))]

  results_list <- list()

  for (test in cross_disease_tests) {
    restricted <- data[get(test$restrict_col) == test$restrict_val]
    restricted[, Batch := droplevels(Batch)]

    n_restricted <- nrow(restricted)
    n_outcome_cases <- sum(restricted[[test$outcome]] == 1, na.rm = TRUE)
    n_outcome_controls <- sum(restricted[[test$outcome]] == 0, na.rm = TRUE)

    cat(sprintf("  %s: n=%d (cases=%d, controls=%d)\n",
                test$label, n_restricted, n_outcome_cases, n_outcome_controls))

    if (n_outcome_cases < 10 || n_outcome_controls < 10) {
      cat(sprintf("  SKIP %s: too few cases (%d) or controls (%d)\n",
                  test$label, n_outcome_cases, n_outcome_controls))
      next
    }

    for (prs_col in prs_cols) {
      cluster_id <- gsub("PRS_", "", prs_col)

      # --- Continuous (per-SD) model ---
      model_data <- restricted[, c(test$outcome, prs_col, covar_terms), with = FALSE]
      model_data <- na.omit(model_data)
      model_data[[prs_col]] <- scale(model_data[[prs_col]])[, 1]

      n_cases <- sum(model_data[[test$outcome]] == 1)
      n_controls <- sum(model_data[[test$outcome]] == 0)

      if (n_cases >= 10 && n_controls >= 10) {
        fit <- tryCatch(
          glm(as.formula(paste(test$outcome, "~", paste(c(prs_col, covar_terms), collapse = " + "))),
              data = model_data, family = binomial),
          error = function(e) { cat(sprintf("  ERROR %s ~ %s (cont): %s\n", test$label, prs_col, e$message)); NULL }
        )
        if (!is.null(fit)) {
          cs <- summary(fit)$coefficients
          if (prs_col %in% rownames(cs)) {
            prs_row <- cs[prs_col, ]
            results_list[[length(results_list) + 1]] <- data.table(
              group = group_name, restrict_disease = test$restrict_col,
              outcome = test$outcome, label = test$label,
              analysis_type = "continuous", predictor = prs_col, cluster = cluster_id,
              OR = round(exp(prs_row["Estimate"]), 4),
              CI_lower = round(exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"]), 4),
              CI_upper = round(exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"]), 4),
              p_value = prs_row["Pr(>|z|)"],
              n_restricted = n_restricted, n_cases = n_cases, n_controls = n_controls
            )
            cat(sprintf("    %s ~ %s (per-SD): OR=%.3f [%.3f-%.3f] P=%.2e\n",
                        test$label, cluster_id,
                        exp(prs_row["Estimate"]),
                        exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"]),
                        exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"]),
                        prs_row["Pr(>|z|)"]))
          }
        }
      }

      # --- Top 10% model ---
      model_data_q <- copy(restricted[, c(test$outcome, prs_col, covar_terms), with = FALSE])
      model_data_q <- na.omit(model_data_q)

      thresh <- quantile(model_data_q[[prs_col]], probs = 0.90, na.rm = TRUE)
      model_data_q[, top_10pct := fifelse(get(prs_col) >= thresh, 1L, 0L)]

      n_top <- sum(model_data_q$top_10pct == 1)
      n_cases_top <- sum(model_data_q[[test$outcome]] == 1 & model_data_q$top_10pct == 1)
      n_cases <- sum(model_data_q[[test$outcome]] == 1)
      n_controls <- sum(model_data_q[[test$outcome]] == 0)

      if (n_cases >= 10 && n_controls >= 10 && n_cases_top >= 2) {
        fit_q <- tryCatch(
          glm(as.formula(paste(test$outcome, "~ top_10pct +", paste(covar_terms, collapse = " + "))),
              data = model_data_q, family = binomial),
          error = function(e) { cat(sprintf("  ERROR %s ~ %s (top10): %s\n", test$label, prs_col, e$message)); NULL }
        )
        if (!is.null(fit_q)) {
          cs <- summary(fit_q)$coefficients
          if ("top_10pct" %in% rownames(cs)) {
            prs_row <- cs["top_10pct", ]
            results_list[[length(results_list) + 1]] <- data.table(
              group = group_name, restrict_disease = test$restrict_col,
              outcome = test$outcome, label = test$label,
              analysis_type = "top_10pct", predictor = prs_col, cluster = cluster_id,
              OR = round(exp(prs_row["Estimate"]), 4),
              CI_lower = round(exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"]), 4),
              CI_upper = round(exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"]), 4),
              p_value = prs_row["Pr(>|z|)"],
              n_restricted = n_restricted, n_cases = n_cases, n_controls = n_controls
            )
            cat(sprintf("    %s ~ %s (top10%%): OR=%.3f [%.3f-%.3f] P=%.2e (n_top=%d, cases_top=%d)\n",
                        test$label, cluster_id,
                        exp(prs_row["Estimate"]),
                        exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"]),
                        exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"]),
                        prs_row["Pr(>|z|)"], n_top, n_cases_top))
          }
        }
      } else {
        cat(sprintf("    SKIP %s ~ %s (top10%%): insufficient (cases=%d, controls=%d, cases_top=%d)\n",
                    test$label, cluster_id, n_cases, n_controls, n_cases_top))
      }
    }
  }

  rbindlist(results_list)
}

# --- Run on all validation groups ---
cat("\n--- Running cross-disease associations ---\n")

all_results <- rbindlist(lapply(val_groups, function(g) {
  gdata <- merged[group == g]
  if (nrow(gdata) == 0) {
    cat(sprintf("  SKIP %s: no individuals\n", g))
    return(data.table())
  }
  run_cross_disease(gdata, g)
}))

# --- Save CSV ---
cat("\n--- Saving CSV results ---\n")
csv_path <- file.path(results_dir, "cross_disease_pgs_associations.csv")
fwrite(all_results, csv_path)
cat(sprintf("  %d rows -> %s\n", nrow(all_results), csv_path))

# --- Build supplementary Excel table (EUR validation only) ---
cat("\n--- Generating supplementary table ---\n")

eur_results <- all_results[group == "eur_validation"]
if (nrow(eur_results) == 0) {
  cat("  No EUR validation results — skipping Excel table.\n")
  quit(status = 0)
}

out_dir <- "results/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reshape_section <- function(dt) {
  label_order <- c("T2D among CAD cases", "CAD among T2D cases")
  dt[, label := factor(label, levels = label_order)]
  setorder(dt, label)

  out_list <- vector("list", length(label_order))
  for (i in seq_along(label_order)) {
    sub <- dt[label == label_order[i]]
    row_data <- list(Outcome = label_order[i])
    for (k in cluster_keys) {
      krow <- sub[cluster == k]
      if (nrow(krow) == 1) {
        row_data[[paste0(k, "_OR")]] <- format_or(krow$OR)
        row_data[[paste0(k, "_CI")]] <- format_ci(krow$CI_lower, krow$CI_upper)
        row_data[[paste0(k, "_P")]]  <- format_pval(krow$p_value)
      } else {
        row_data[[paste0(k, "_OR")]] <- ""
        row_data[[paste0(k, "_CI")]] <- ""
        row_data[[paste0(k, "_P")]]  <- ""
      }
    }
    out_list[[i]] <- as.data.table(row_data)
  }
  rbindlist(out_list)
}

wide_cont <- reshape_section(eur_results[analysis_type == "continuous"])
wide_top10 <- reshape_section(eur_results[analysis_type == "top_10pct"])

total_cols <- 1L + length(cluster_keys) * 3L  # Outcome + 10 clusters * 3

wb <- createWorkbook()
addWorksheet(wb, "Data")

title_style  <- createStyle(textDecoration = "bold", fontSize = 11)
section_style <- createStyle(textDecoration = "bold", fontSize = 10, fgFill = "#D9D9D9")
cluster_style <- createStyle(textDecoration = "bold", halign = "center", border = "bottom")
sub_style    <- createStyle(textDecoration = "bold", halign = "center", border = "bottom")

# Row 1: title
title_text <- "Supplementary Table X: Association of cluster polygenic scores with cross-disease cardiometabolic risk in the European validation cohort"
writeData(wb, 1, x = title_text, startRow = 1, startCol = 1)
addStyle(wb, 1, title_style, rows = 1, cols = 1)
mergeCells(wb, 1, cols = 1:total_cols, rows = 1)

# --- Section A: Continuous (per-SD) ---
# Row 3: section header
writeData(wb, 1, x = "A. Full PGS distribution (per-SD odds ratios)", startRow = 3, startCol = 1)
addStyle(wb, 1, section_style, rows = 3, cols = 1:total_cols, gridExpand = TRUE)
mergeCells(wb, 1, cols = 1:total_cols, rows = 3)

# Row 4-5: cluster headers + sub-headers
writeData(wb, 1, x = "Outcome", startRow = 4, startCol = 1)
addStyle(wb, 1, sub_style, rows = 4, cols = 1)
mergeCells(wb, 1, cols = 1, rows = 4:5)

for (k_idx in seq_along(cluster_keys)) {
  k <- cluster_keys[k_idx]
  col_start <- 1L + (k_idx - 1L) * 3L + 1L
  writeData(wb, 1, x = cluster_labels[k], startRow = 4, startCol = col_start)
  addStyle(wb, 1, cluster_style, rows = 4, cols = col_start)
  mergeCells(wb, 1, cols = col_start:(col_start + 2L), rows = 4)

  writeData(wb, 1, x = "OR/SD", startRow = 5, startCol = col_start)
  addStyle(wb, 1, sub_style, rows = 5, cols = col_start)
  writeData(wb, 1, x = "95% CI", startRow = 5, startCol = col_start + 1L)
  addStyle(wb, 1, sub_style, rows = 5, cols = col_start + 1L)
  writeData(wb, 1, x = "P", startRow = 5, startCol = col_start + 2L)
  addStyle(wb, 1, sub_style, rows = 5, cols = col_start + 2L)
}

# Rows 6-7: continuous data
writeData(wb, 1, x = wide_cont, startRow = 6, startCol = 1, colNames = FALSE)

# --- Section B: Top 10% ---
# Row 9: section header
writeData(wb, 1, x = "B. Top 10% of PGS distribution", startRow = 9, startCol = 1)
addStyle(wb, 1, section_style, rows = 9, cols = 1:total_cols, gridExpand = TRUE)
mergeCells(wb, 1, cols = 1:total_cols, rows = 9)

# Row 10-11: cluster headers + sub-headers
writeData(wb, 1, x = "Outcome", startRow = 10, startCol = 1)
addStyle(wb, 1, sub_style, rows = 10, cols = 1)
mergeCells(wb, 1, cols = 1, rows = 10:11)

for (k_idx in seq_along(cluster_keys)) {
  k <- cluster_keys[k_idx]
  col_start <- 1L + (k_idx - 1L) * 3L + 1L
  writeData(wb, 1, x = cluster_labels[k], startRow = 10, startCol = col_start)
  addStyle(wb, 1, cluster_style, rows = 10, cols = col_start)
  mergeCells(wb, 1, cols = col_start:(col_start + 2L), rows = 10)

  writeData(wb, 1, x = "OR", startRow = 11, startCol = col_start)
  addStyle(wb, 1, sub_style, rows = 11, cols = col_start)
  writeData(wb, 1, x = "95% CI", startRow = 11, startCol = col_start + 1L)
  addStyle(wb, 1, sub_style, rows = 11, cols = col_start + 1L)
  writeData(wb, 1, x = "P", startRow = 11, startCol = col_start + 2L)
  addStyle(wb, 1, sub_style, rows = 11, cols = col_start + 2L)
}

# Rows 12-13: top 10% data
writeData(wb, 1, x = wide_top10, startRow = 12, startCol = 1, colNames = FALSE)

# Column widths
setColWidths(wb, 1, cols = 1, widths = 22)
or_cols <- 1L + seq(1, 30, by = 3)
ci_cols <- 1L + seq(2, 30, by = 3)
p_cols  <- 1L + seq(3, 30, by = 3)
setColWidths(wb, 1, cols = or_cols, widths = 8)
setColWidths(wb, 1, cols = ci_cols, widths = 14)
setColWidths(wb, 1, cols = p_cols, widths = 12)

# Legend
legend_lines <- c(
  "Legend:",
  "OR/SD = odds ratio per one standard deviation increase in cluster PGS.",
  "OR = odds ratio for top 10% vs. remaining 90% of the cluster PGS distribution.",
  "95% CI = 95% confidence interval.",
  paste0("P = P-value from logistic regression adjusted for age, age², sex, ",
         "genotyping batch, and PC1–10."),
  paste0("\"T2D among CAD cases\" = association of cluster PGS with T2D restricted ",
         "to individuals with prevalent CAD."),
  paste0("\"CAD among T2D cases\" = association of cluster PGS with CAD restricted ",
         "to individuals with prevalent T2D."),
  paste0("For the top 10% analysis, decile cutoffs were computed within each ",
         "disease-restricted subset."),
  "",
  "Abbreviations:",
  paste0("PGS = polygenic score; T2D = type 2 diabetes; CAD = coronary artery disease; ",
         "OR = odds ratio; CI = confidence interval; SD = standard deviation; ",
         "ALP = alkaline phosphatase; SHBG = sex hormone-binding globulin; ",
         "Lpa = lipoprotein(a).")
)

legend_start <- 15L
legend_style <- createStyle(fontSize = 9, wrapText = TRUE)
bold_legend  <- createStyle(textDecoration = "bold", fontSize = 9)

for (i in seq_along(legend_lines)) {
  r <- legend_start + i - 1L
  writeData(wb, 1, x = legend_lines[i], startRow = r, startCol = 1)
  sty <- if (i %in% c(1, 10)) bold_legend else legend_style
  addStyle(wb, 1, style = sty, rows = r, cols = 1)
  mergeCells(wb, 1, cols = 1:total_cols, rows = r)
}

out_path <- file.path(out_dir, "cluster_pgs_cross_disease_associations.xlsx")
saveWorkbook(wb, out_path, overwrite = TRUE)
cat(sprintf("  Excel table -> %s\n", out_path))

# --- Summary ---
cat("\n=== Summary ===\n")
cat(sprintf("  CSV: %s (%d rows)\n", csv_path, nrow(all_results)))
cat(sprintf("  Excel: %s\n", out_path))
if (nrow(eur_results) > 0) {
  cat("\n  EUR validation results:\n")
  print(eur_results[, .(label, analysis_type, cluster, OR, CI_lower, CI_upper, p_value,
                        n_restricted, n_cases, n_controls)])
}

cat("\n=== Done ===\n")
