#!/usr/bin/env Rscript
# make_supptable_cluster_pgs.R
# Generate supplementary Excel tables for individual cluster PGS associations
# with T2D, CAD, and T2D/CAD composite outcomes.

library(data.table)
library(openxlsx)
library(yaml)

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

cfg <- read_yaml("config/b1_config.yaml")
cluster_labels <- unlist(cfg$cluster_labels)
cluster_keys <- paste0("K", 1:10)

outcome_display <- c(T2D = "T2D", CAD = "CAD", Synthetic = "T2D/CAD")
outcome_order <- c("T2D", "CAD", "Synthetic")
ancestry_display <- c(afr_validation = "AFR", eas_validation = "EAS", sas_validation = "SAS")
ancestry_order <- c("afr_validation", "eas_validation", "sas_validation")

# Default to the soft (probabilistic) cluster-PGS outputs; --hard-assignment
# switches to the hard-assignment subdirectory (mirrors the b1 analysis scripts).
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if ("--hard-assignment" %in% args) "results/b1_analysis/prs_hard_assignment" else "results/b1_analysis"
out_dir <- "results/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Formatting helpers -------------------------------------------------------

format_or <- function(x) sprintf("%.2f", x)

format_ci <- function(lower, upper) sprintf("%.2f–%.2f", lower, upper)

format_pval <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, sprintf("%.2e", p),
      ifelse(p < 0.01, sprintf("%.3f", p),
        sprintf("%.2f", p))))
}

# --- Reshape long -> wide -----------------------------------------------------

reshape_to_wide <- function(dt, include_ancestry = FALSE) {
  dt <- copy(dt)
  dt[, outcome := factor(outcome, levels = outcome_order)]

  if (include_ancestry) {
    dt[, anc := factor(ancestry_display[group], levels = c("AFR", "EAS", "SAS"))]
    setorder(dt, anc, outcome)
    rows <- unique(dt[, .(group, outcome, anc)])
  } else {
    setorder(dt, outcome)
    rows <- unique(dt[, .(outcome)])
  }

  out_list <- vector("list", nrow(rows))
  for (i in seq_len(nrow(rows))) {
    row_data <- list()
    if (include_ancestry) {
      sub <- dt[group == rows$group[i] & outcome == rows$outcome[i]]
      row_data[["Ancestry"]] <- as.character(rows$anc[i])
      row_data[["Outcome"]] <- outcome_display[as.character(rows$outcome[i])]
    } else {
      sub <- dt[outcome == rows$outcome[i]]
      row_data[["Outcome"]] <- outcome_display[as.character(rows$outcome[i])]
    }
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

# --- Build workbook -----------------------------------------------------------

build_workbook <- function(wide_dt, title_text, or_col_name, include_ancestry) {
  wb <- createWorkbook()
  addWorksheet(wb, "Data")

  n_label_cols <- if (include_ancestry) 2L else 1L
  total_cols <- n_label_cols + length(cluster_keys) * 3L

  # Row 1: title
  title_style <- createStyle(textDecoration = "bold", fontSize = 11)
  writeData(wb, 1, x = title_text, startRow = 1, startCol = 1)
  addStyle(wb, 1, title_style, rows = 1, cols = 1)
  mergeCells(wb, 1, cols = 1:total_cols, rows = 1)

  # Row 3-4: multi-level headers
  cluster_style <- createStyle(textDecoration = "bold", halign = "center", border = "bottom")
  sub_style <- createStyle(textDecoration = "bold", halign = "center", border = "bottom")

  # Label columns span rows 3-4
  if (include_ancestry) {
    writeData(wb, 1, x = "Ancestry", startRow = 3, startCol = 1)
    addStyle(wb, 1, sub_style, rows = 3, cols = 1)
    mergeCells(wb, 1, cols = 1, rows = 3:4)
    writeData(wb, 1, x = "Outcome", startRow = 3, startCol = 2)
    addStyle(wb, 1, sub_style, rows = 3, cols = 2)
    mergeCells(wb, 1, cols = 2, rows = 3:4)
  } else {
    writeData(wb, 1, x = "Outcome", startRow = 3, startCol = 1)
    addStyle(wb, 1, sub_style, rows = 3, cols = 1)
    mergeCells(wb, 1, cols = 1, rows = 3:4)
  }

  # Cluster headers (row 3, merged across 3 cols) and sub-headers (row 4)
  for (k_idx in seq_along(cluster_keys)) {
    k <- cluster_keys[k_idx]
    col_start <- n_label_cols + (k_idx - 1L) * 3L + 1L
    writeData(wb, 1, x = cluster_labels[k], startRow = 3, startCol = col_start)
    addStyle(wb, 1, cluster_style, rows = 3, cols = col_start)
    mergeCells(wb, 1, cols = col_start:(col_start + 2L), rows = 3)

    writeData(wb, 1, x = or_col_name, startRow = 4, startCol = col_start)
    addStyle(wb, 1, sub_style, rows = 4, cols = col_start)
    writeData(wb, 1, x = "95% CI", startRow = 4, startCol = col_start + 1L)
    addStyle(wb, 1, sub_style, rows = 4, cols = col_start + 1L)
    writeData(wb, 1, x = "P", startRow = 4, startCol = col_start + 2L)
    addStyle(wb, 1, sub_style, rows = 4, cols = col_start + 2L)
  }

  # Data rows
  writeData(wb, 1, x = wide_dt, startRow = 5, startCol = 1, colNames = FALSE)

  # Column widths
  if (include_ancestry) {
    setColWidths(wb, 1, cols = 1, widths = 10)
    setColWidths(wb, 1, cols = 2, widths = 12)
  } else {
    setColWidths(wb, 1, cols = 1, widths = 12)
  }
  or_cols <- n_label_cols + seq(1, 30, by = 3)
  ci_cols <- n_label_cols + seq(2, 30, by = 3)
  p_cols  <- n_label_cols + seq(3, 30, by = 3)
  setColWidths(wb, 1, cols = or_cols, widths = 8)
  setColWidths(wb, 1, cols = ci_cols, widths = 14)
  setColWidths(wb, 1, cols = p_cols, widths = 12)

  # Legend
  is_quantile <- or_col_name == "OR"
  or_desc <- if (is_quantile) {
    "OR = odds ratio for top 10% vs. remaining 90% of the cluster PGS distribution."
  } else {
    "OR/SD = odds ratio per one standard deviation increase in cluster PGS."
  }
  legend_lines <- c(
    "Legend:",
    or_desc,
    "95% CI = 95% confidence interval.",
    "P = P-value from logistic regression adjusted for age, age², sex, genotyping batch, and PC1–10.",
    "T2D/CAD = composite outcome defined as having both T2D and CAD.",
    "",
    "Abbreviations:",
    paste0("PGS = polygenic score; T2D = type 2 diabetes; CAD = coronary artery disease; ",
           "OR = odds ratio; CI = confidence interval; SD = standard deviation; ",
           "AFR = African; EAS = East Asian; SAS = South Asian; EUR = European; ",
           "ALP = alkaline phosphatase; SHBG = sex hormone-binding globulin; ",
           "Lpa = lipoprotein(a).")
  )
  legend_start <- nrow(wide_dt) + 7
  legend_style <- createStyle(fontSize = 9, wrapText = TRUE)
  bold_legend <- createStyle(textDecoration = "bold", fontSize = 9)
  for (i in seq_along(legend_lines)) {
    r <- legend_start + i - 1L
    writeData(wb, 1, x = legend_lines[i], startRow = r, startCol = 1)
    sty <- if (i %in% c(1, 7)) bold_legend else legend_style
    addStyle(wb, 1, style = sty, rows = r, cols = 1)
    mergeCells(wb, 1, cols = 1:total_cols, rows = r)
  }

  wb
}

# --- Load data ----------------------------------------------------------------

cat("Loading data...\n")

eur_cont <- fread(file.path(data_dir, "association_results_eur_validation.csv"))
eur_cont <- eur_cont[model_type == "individual"]
eur_cont[, cluster := gsub("PRS_", "", predictor)]

noneur_cont <- rbindlist(lapply(c("afr", "eas", "sas"), function(a) {
  fread(file.path(data_dir, sprintf("association_results_%s_validation.csv", a)))
}))
noneur_cont <- noneur_cont[model_type == "individual"]
noneur_cont[, cluster := gsub("PRS_", "", predictor)]

quant <- fread(file.path(data_dir, "quantile_prs_associations.csv"))
quant <- quant[quantile == "top_10pct" & grepl("^K\\d+$", cluster)]
eur_quant <- quant[group == "eur_validation"]
noneur_quant <- quant[group %in% paste0(c("afr", "eas", "sas"), "_validation")]

# --- Generate tables ----------------------------------------------------------

cat("Generating EUR validation table...\n")
wide <- reshape_to_wide(eur_cont, include_ancestry = FALSE)
wb <- build_workbook(wide,
  title_text = "Supplementary Table X: Association of cluster polygenic scores with cardiometabolic outcomes in the European validation cohort",
  or_col_name = "OR/SD", include_ancestry = FALSE)
saveWorkbook(wb, file.path(out_dir, "cluster_pgs_eur_validation.xlsx"), overwrite = TRUE)

cat("Generating non-EUR validation table...\n")
wide <- reshape_to_wide(noneur_cont, include_ancestry = TRUE)
wb <- build_workbook(wide,
  title_text = "Supplementary Table X: Association of cluster polygenic scores with cardiometabolic outcomes in non-European validation cohorts",
  or_col_name = "OR/SD", include_ancestry = TRUE)
saveWorkbook(wb, file.path(out_dir, "cluster_pgs_noneur_validation.xlsx"), overwrite = TRUE)

cat("Generating EUR top 10%% table...\n")
wide <- reshape_to_wide(eur_quant, include_ancestry = FALSE)
wb <- build_workbook(wide,
  title_text = "Supplementary Table X: Association of top 10% cluster polygenic scores with cardiometabolic outcomes in the European validation cohort",
  or_col_name = "OR", include_ancestry = FALSE)
saveWorkbook(wb, file.path(out_dir, "cluster_pgs_eur_top10pct.xlsx"), overwrite = TRUE)

cat("Generating non-EUR top 10%% table...\n")
wide <- reshape_to_wide(noneur_quant, include_ancestry = TRUE)
wb <- build_workbook(wide,
  title_text = "Supplementary Table X: Association of top 10% cluster polygenic scores with cardiometabolic outcomes in non-European validation cohorts",
  or_col_name = "OR", include_ancestry = TRUE)
saveWorkbook(wb, file.path(out_dir, "cluster_pgs_noneur_top10pct.xlsx"), overwrite = TRUE)

cat("Done. Files saved to results/supplementary_tables/\n")
