#!/usr/bin/env Rscript
# make_supptable_individual_percentiles.R
# Supplementary table: individual-level cluster PGS percentiles with
# PCE risk, PREVENT risk, and HbA1c for the EUR validation cohort.
#
# Usage:
#   Rscript scripts/b1_analysis/make_supptable_individual_percentiles.R --config config/b1_config.yaml --hard-assignment

library(data.table)
library(openxlsx)
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
}

out_dir <- "results/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Supplementary Table: Individual Cluster PGS Percentiles ===\n")

# ============================================================
# Load data
# ============================================================
cat("\n--- Loading data ---\n")

prs_all <- fread(file.path(results_dir, "prs/cluster_prs_all.tsv"))
biomarker <- fread("results/b3_1_analysis/biomarker_phenotypes.tsv",
                   select = c("FID", "IID", "PCE_ASCVD", "PREVENT_ASCVD", "HbA1c_corrected"))

eur_val <- prs_all[group == "eur_validation"]
cat(sprintf("  EUR validation: %d individuals\n", nrow(eur_val)))

# Remove negative FIDs
n_neg <- sum(eur_val$FID < 0)
eur_val <- eur_val[FID > 0]
cat(sprintf("  Removed %d negative-FID individuals -> %d remaining\n", n_neg, nrow(eur_val)))

cluster_labels <- unlist(cfg$cluster_labels)
prs_cols <- grep("^PRS_K", colnames(eur_val), value = TRUE)
cat(sprintf("  Clusters: %s\n", paste(prs_cols, collapse = ", ")))

# ============================================================
# Merge biomarkers and filter to complete cases
# ============================================================
cat("\n--- Merging biomarkers ---\n")

dat <- merge(eur_val, biomarker, by = c("FID", "IID"), all.x = TRUE)
dat[, HbA1c_pct := 0.09148 * HbA1c_corrected + 2.152]

cat(sprintf("  PCE available: %d / %d\n",
            sum(!is.na(dat$PCE_ASCVD)), nrow(dat)))
cat(sprintf("  PREVENT available: %d / %d\n",
            sum(!is.na(dat$PREVENT_ASCVD)), nrow(dat)))
cat(sprintf("  HbA1c available: %d / %d\n",
            sum(!is.na(dat$HbA1c_pct)), nrow(dat)))

# Keep only individuals with complete PCE, PREVENT, and HbA1c
dat <- dat[!is.na(PCE_ASCVD) & !is.na(PREVENT_ASCVD) & !is.na(HbA1c_pct)]
cat(sprintf("  After requiring complete PCE + PREVENT + HbA1c: %d individuals\n", nrow(dat)))

# ============================================================
# Compute percentiles (within the filtered set)
# ============================================================
cat("\n--- Computing percentiles ---\n")

n <- nrow(dat)
for (pc in prs_cols) {
  pctile_col <- gsub("PRS_", "pctile_", pc)
  dat[, (pctile_col) := frank(get(pc), ties.method = "average") / n * 100]
}

# ============================================================
# Build output table
# ============================================================
cat("--- Building output table ---\n")

# Cluster column order: canonical b1 display order (Glycemic -> Lpa), matching
# scripts/b1_analysis/03_b1_forest_plots.R. Clusters not in the list fall back
# to natural K order at the end.
desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
all_k <- names(cluster_labels)
k_ordered <- c(desired_k_order[desired_k_order %in% all_k],
               setdiff(all_k, desired_k_order))

pctile_cols <- paste0("pctile_", k_ordered)
col_labels <- sprintf("%s (%s)", k_ordered, cluster_labels[k_ordered])

out <- dat[, c("FID", "IID", pctile_cols, "PCE_ASCVD", "PREVENT_ASCVD", "HbA1c_pct"),
           with = FALSE]

setnames(out, pctile_cols, col_labels)
setnames(out, "PCE_ASCVD", "PCE Risk (%)")
setnames(out, "PREVENT_ASCVD", "PREVENT Risk (%)")
setnames(out, "HbA1c_pct", "HbA1c (%)")

setorder(out, -`PCE Risk (%)`)

# Round values
for (cl in col_labels) {
  out[, (cl) := round(get(cl), 1)]
}
out[, `PCE Risk (%)` := round(`PCE Risk (%)`, 2)]
out[, `PREVENT Risk (%)` := round(`PREVENT Risk (%)`, 2)]
out[, `HbA1c (%)` := round(`HbA1c (%)`, 2)]

cat(sprintf("  Output: %d rows, %d columns\n", nrow(out), ncol(out)))

# ============================================================
# Write Excel
# ============================================================
cat("\n--- Writing Excel ---\n")

wb <- createWorkbook()
addWorksheet(wb, "Data")

total_cols <- ncol(out)

title_text <- paste0(
  "Supplementary Table X: Individual-level cluster PGS percentiles and ",
  "cardiovascular risk in the European validation cohort")

title_style <- createStyle(textDecoration = "bold", fontSize = 11)
writeData(wb, 1, x = title_text, startRow = 1, startCol = 1)
addStyle(wb, 1, title_style, rows = 1, cols = 1)
mergeCells(wb, 1, cols = 1:total_cols, rows = 1)

header_style <- createStyle(textDecoration = "bold", halign = "center",
                            border = "bottom", borderStyle = "thin")
for (j in seq_along(colnames(out))) {
  writeData(wb, 1, x = colnames(out)[j], startRow = 3, startCol = j)
  addStyle(wb, 1, header_style, rows = 3, cols = j)
}

writeData(wb, 1, x = out, startRow = 4, startCol = 1, colNames = FALSE)

setColWidths(wb, 1, cols = 1, widths = 12)
setColWidths(wb, 1, cols = 2, widths = 12)
setColWidths(wb, 1, cols = 3:12, widths = 14)
setColWidths(wb, 1, cols = 13, widths = 14)
setColWidths(wb, 1, cols = 14, widths = 16)
setColWidths(wb, 1, cols = 15, widths = 12)

# Legend
legend_lines <- c(
  "Legend:",
  "Each row represents one individual in the European validation cohort.",
  "Cluster PGS percentiles are computed within the included individuals (0 = lowest, 100 = highest).",
  "Individuals with negative FID or missing PCE, PREVENT, or HbA1c values were excluded.",
  "PCE Risk = 10-year atherosclerotic cardiovascular disease risk estimated using the Pooled Cohort Equations (Goff et al., 2014).",
  "PREVENT Risk = 10-year ASCVD risk estimated using the AHA PREVENT equations (Khan et al., 2024).",
  "HbA1c is converted to NGSP units (%) from IFCC units (mmol/mol).",
  "Rows are ordered by PCE Risk in descending order.",
  "",
  "Abbreviations:",
  paste0("PGS = polygenic score; PCE = Pooled Cohort Equations; ",
         "PREVENT = Predicting Risk of cardiovascular disease EVENTs; ",
         "ASCVD = atherosclerotic cardiovascular disease; ",
         "HbA1c = glycated hemoglobin; ",
         "ALP = alkaline phosphatase; SHBG = sex hormone-binding globulin; ",
         "Lpa = lipoprotein(a).")
)

legend_start <- nrow(out) + 6
legend_style <- createStyle(fontSize = 9, wrapText = TRUE)
bold_legend <- createStyle(textDecoration = "bold", fontSize = 9)
for (i in seq_along(legend_lines)) {
  r <- legend_start + i - 1L
  writeData(wb, 1, x = legend_lines[i], startRow = r, startCol = 1)
  sty <- if (i %in% c(1, 9)) bold_legend else legend_style
  addStyle(wb, 1, style = sty, rows = r, cols = 1)
  mergeCells(wb, 1, cols = 1:total_cols, rows = r)
}

out_file <- file.path(out_dir, "cluster_pgs_individual_percentiles.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("  Saved: %s\n", out_file))

cat("\n=== Done ===\n")
