#!/usr/bin/env Rscript
# 03_supplementary_table.R
# Generate Excel supplementary table for B3 cluster PGS biomarker associations.
# Each cluster is expanded into three rows: Male, Female, then All (combined).
#
# Usage:
#   Rscript scripts/b3_1_analysis/03_supplementary_table.R
#   Rscript scripts/b3_1_analysis/03_supplementary_table.R --hard-assignment

library(data.table)
library(openxlsx)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
hard_assignment <- "--hard-assignment" %in% args

cfg <- read_yaml("config/b3_1_config.yaml")
results_dir <- cfg$results_dir
if (hard_assignment) {
  results_dir <- "results/b3_1_analysis/prs_hard_assignment"
  cat("  ** HARD ASSIGNMENT MODE **\n")
}

out_dir <- "results/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data (keep all sex groups) ---
indiv <- fread(file.path(results_dir, "biomarker_association_individual.csv"))
indiv <- indiv[model_type == "individual"]

# --- Cluster labels (ordered K1-K10) ---
cluster_labels <- unlist(cfg$cluster_labels)
k_order <- paste0("K", seq_along(cluster_labels))
cluster_map <- setNames(cluster_labels[k_order], paste0("PRS_", k_order))

# --- Sex strata: Male, Female, then All (combined) per cluster ---
sex_order <- c("Male", "Female", "All")
n_sex <- length(sex_order)

# --- 12 cardiometabolic traits (exclude PREVENT-ASCVD) ---
trait_order <- c(
  "HbA1c_corrected", "glucose_corrected",
  "Lpa", "non_HDL_chol", "triglycerides_corrected", "LDL_corrected",
  "WHtR", "waist_circumference", "BMI",
  "CRP", "eGFR", "SBP_corrected"
)
trait_display <- c(
  HbA1c_corrected = "HbA1c", glucose_corrected = "Fasting Glucose",
  Lpa = "Lpa", non_HDL_chol = "non-HDL Cholesterol",
  triglycerides_corrected = "Triglycerides", LDL_corrected = "LDL",
  WHtR = "WHtR", waist_circumference = "Waist Circumference", BMI = "BMI",
  CRP = "CRP", eGFR = "eGFR", SBP_corrected = "SBP"
)

indiv <- indiv[trait %in% trait_order]

# --- Sample sizes per trait (from the combined "All" group) ---
sample_sizes <- indiv[sex_group == "All" & predictor == paste0("PRS_", k_order[1]),
                      .(n = n_obs[1]), by = trait]
n_map <- setNames(sample_sizes$n, sample_sizes$trait)

# --- Build wide data matrix: rows = clusters x sex (Male/Female/All), cols = (beta,se,p) x traits ---
n_traits   <- length(trait_order)
n_clusters <- length(k_order)
n_rows     <- n_clusters * n_sex

data_mat <- matrix(NA_character_, nrow = n_rows, ncol = n_traits * 3)

for (i in seq_along(k_order)) {
  prs_col <- paste0("PRS_", k_order[i])
  for (s in seq_along(sex_order)) {
    r <- (i - 1) * n_sex + s
    for (j in seq_along(trait_order)) {
      row <- indiv[predictor == prs_col & trait == trait_order[j] &
                     sex_group == sex_order[s]]
      if (nrow(row) == 1) {
        col_base <- (j - 1) * 3
        data_mat[r, col_base + 1] <- formatC(row$beta, format = "f", digits = 4)
        data_mat[r, col_base + 2] <- formatC(row$se, format = "f", digits = 4)
        data_mat[r, col_base + 3] <- formatC(row$p_value, format = "e", digits = 2)
      }
    }
  }
}

# --- Build Excel workbook ---
n_label_cols <- 2L                 # Cluster + Sex
total_cols   <- n_label_cols + n_traits * 3

wb <- createWorkbook()
addWorksheet(wb, "Table")

# Styles
title_style   <- createStyle(fontSize = 12, textDecoration = "bold")
trait_style   <- createStyle(fontSize = 11, textDecoration = "bold", halign = "center",
                             border = "bottom", borderStyle = "thin")
sub_style     <- createStyle(fontSize = 10, textDecoration = "bold", halign = "center",
                             border = "bottom", borderStyle = "thin")
cluster_style <- createStyle(fontSize = 10, textDecoration = "bold", valign = "center")
sex_style     <- createStyle(fontSize = 10, halign = "center")
data_style    <- createStyle(fontSize = 10, halign = "center")
legend_style  <- createStyle(fontSize = 9, wrapText = TRUE)

# Row 1: Title
title_text <- "Supplementary Table X: Validation of individual-level cluster PGS associations with cardiometabolic risk factors, by sex"
writeData(wb, "Table", title_text, startCol = 1, startRow = 1)
addStyle(wb, "Table", title_style, rows = 1, cols = 1)
mergeCells(wb, "Table", cols = 1:total_cols, rows = 1)

# Rows 3-4: Cluster + Sex headers (merged across rows 3:4)
writeData(wb, "Table", "Cluster", startCol = 1, startRow = 3)
mergeCells(wb, "Table", cols = 1, rows = 3:4)
addStyle(wb, "Table", trait_style, rows = 3:4, cols = 1, gridExpand = TRUE)
writeData(wb, "Table", "Sex", startCol = 2, startRow = 3)
mergeCells(wb, "Table", cols = 2, rows = 3:4)
addStyle(wb, "Table", trait_style, rows = 3:4, cols = 2, gridExpand = TRUE)

# Rows 3-4: Trait headers (merged 3 cols) + Beta/SE/P sub-headers
for (j in seq_along(trait_order)) {
  col_start <- n_label_cols + 1 + (j - 1) * 3
  n <- format(n_map[trait_order[j]], big.mark = ",")
  header <- sprintf("%s (n=%s)", trait_display[trait_order[j]], n)
  writeData(wb, "Table", header, startCol = col_start, startRow = 3)
  mergeCells(wb, "Table", cols = col_start:(col_start + 2), rows = 3)
  addStyle(wb, "Table", trait_style, rows = 3, cols = col_start:(col_start + 2),
           gridExpand = TRUE)

  writeData(wb, "Table", "Beta", startCol = col_start, startRow = 4)
  writeData(wb, "Table", "SE", startCol = col_start + 1, startRow = 4)
  writeData(wb, "Table", "P-value", startCol = col_start + 2, startRow = 4)
  addStyle(wb, "Table", sub_style, rows = 4, cols = col_start:(col_start + 2),
           gridExpand = TRUE)
}

# Data: 3 rows (Male/Female/All) per cluster
for (i in seq_along(k_order)) {
  base_row <- 4 + (i - 1) * n_sex + 1            # first data row for this cluster
  # Cluster label spans the 3 sex rows
  writeData(wb, "Table", cluster_map[paste0("PRS_", k_order[i])],
            startCol = 1, startRow = base_row)
  mergeCells(wb, "Table", cols = 1, rows = base_row:(base_row + n_sex - 1))
  addStyle(wb, "Table", cluster_style, rows = base_row:(base_row + n_sex - 1),
           cols = 1, gridExpand = TRUE)

  for (s in seq_along(sex_order)) {
    data_row <- base_row + (s - 1)
    r        <- (i - 1) * n_sex + s
    writeData(wb, "Table", sex_order[s], startCol = 2, startRow = data_row)
    addStyle(wb, "Table", sex_style, rows = data_row, cols = 2)
    for (j in seq(n_traits * 3)) {
      col <- n_label_cols + j
      writeData(wb, "Table", data_mat[r, j], startCol = col, startRow = data_row)
      addStyle(wb, "Table", data_style, rows = data_row, cols = col)
    }
  }
}

# Column widths
setColWidths(wb, "Table", cols = 1, widths = 24)
setColWidths(wb, "Table", cols = 2, widths = 9)
for (j in seq_along(trait_order)) {
  col_start <- n_label_cols + 1 + (j - 1) * 3
  setColWidths(wb, "Table", cols = col_start, widths = 12)
  setColWidths(wb, "Table", cols = col_start + 1, widths = 10)
  setColWidths(wb, "Table", cols = col_start + 2, widths = 12)
}

# Legend
legend_row <- 4 + n_rows + 2
legend_lines <- paste(
  "Abbreviations: PGS, polygenic score; SE, standard error; HbA1c, glycated hemoglobin;",
  "LDL, low-density lipoprotein cholesterol; non-HDL, non-high-density lipoprotein cholesterol;",
  "WHtR, waist-to-height ratio; CRP, C-reactive protein; eGFR, estimated glomerular filtration rate;",
  "SBP, systolic blood pressure.",
  "Each cluster is shown for Male, Female, and All (combined) participants.",
  "Beta coefficients represent the change in trait per 1 SD increase in cluster PGS,",
  "estimated via linear regression adjusting for age, age-squared, sex (in the All model),",
  "genotyping batch, and PC1-10 in the EUR validation cohort.",
  "Sample sizes (n) in the header refer to the combined (All) group; Male/Female are subsets.",
  "Medication corrections were applied to HbA1c, fasting glucose, LDL, total cholesterol,",
  "triglycerides, and SBP prior to analysis.",
  "Lpa, CRP, triglycerides, and fasting glucose were log-transformed before regression."
)
writeData(wb, "Table", legend_lines, startCol = 1, startRow = legend_row)
addStyle(wb, "Table", legend_style, rows = legend_row, cols = 1)
mergeCells(wb, "Table", cols = 1:total_cols, rows = legend_row)

# Save
out_file <- file.path(out_dir, "cluster_pgs_biomarker_associations.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("Saved: %s (%d clusters x %d sex strata = %d data rows)\n",
            out_file, n_clusters, n_sex, n_rows))
