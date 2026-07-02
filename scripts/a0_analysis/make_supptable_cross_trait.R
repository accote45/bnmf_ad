#!/usr/bin/env Rscript
# make_supptable_cross_trait.R
# Generate a supplementary Excel table for the cross-trait PRS regressions
# (T2D ~ CVD PRS and CVD ~ T2D PRS) produced by 03_cross_trait_regression.R.
# Reports N cases, N controls, OR per SD, 95% CI, and P for all 10 models.

library(data.table)
library(openxlsx)

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

in_csv  <- "results/a0_analysis/prs_ct/cross_trait_regression_results.csv"
out_dir <- "results/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_xlsx <- file.path(out_dir, "cross_trait_prs_regression.xlsx")

# --- Load + format ------------------------------------------------------------
res <- fread(in_csv)

# Preserve the analysis order: T2D as outcome (CVD predictors), then CVD outcomes
res[, outcome := factor(outcome, levels = c("T2D", "CAD", "Angina", "MI", "Stroke", "PAD"))]
setorder(res, outcome)
# Keep T2D-outcome block first (predictor != T2D), then CVD-outcome block
res[, block := ifelse(predictor_trait == "T2D", 2L, 1L)]
setorder(res, block, outcome)

prs_threshold <- unique(res$prs_threshold)

tab <- data.table(
  Outcome        = as.character(res$outcome),
  `Predictor PRS` = res$predictor_trait,
  `N cases`      = res$n_cases,
  `N controls`   = res$n_controls,
  `OR/SD`        = round(res$OR, 3),
  `95% CI lower` = round(res$CI_lower, 3),
  `95% CI upper` = round(res$CI_upper, 3),
  `Nagelkerke R2` = signif(res$nagelkerke_r2, 3),
  `P`            = res$p_value
)

cat(sprintf("Loaded %d regressions from %s (PRS threshold p < %s)\n",
            nrow(tab), in_csv, paste(prs_threshold, collapse = ",")))

# --- Build workbook (house style) --------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "Cross-trait PRS")

total_cols <- ncol(tab)

title_text <- paste0(
  "Supplementary Table X: Cross-trait polygenic score associations between ",
  "type 2 diabetes and cardiovascular disease traits in the multi-ancestry analysis cohort")
title_style <- createStyle(textDecoration = "bold", fontSize = 12, wrapText = TRUE,
                           valign = "center")
writeData(wb, 1, x = title_text, startRow = 1, startCol = 1)
addStyle(wb, 1, title_style, rows = 1, cols = 1)
mergeCells(wb, 1, cols = 1:total_cols, rows = 1)
setRowHeights(wb, 1, rows = 1, heights = 42)

# Header row (row 3) and data (row 4+)
header_style <- createStyle(textDecoration = "bold", halign = "center", border = "bottom")
num_p_style  <- createStyle(numFmt = "0.00E+00")
writeData(wb, 1, x = tab, startRow = 3, startCol = 1, headerStyle = header_style)

# Scientific notation for the P column
p_col <- which(colnames(tab) == "P")
addStyle(wb, 1, num_p_style, rows = 4:(3 + nrow(tab)), cols = p_col, gridExpand = TRUE)

# Column widths
setColWidths(wb, 1, cols = 1, widths = 10)        # Outcome
setColWidths(wb, 1, cols = 2, widths = 14)        # Predictor PRS
setColWidths(wb, 1, cols = 3:4, widths = 12)      # N cases / controls
setColWidths(wb, 1, cols = 5:8, widths = 13)      # OR / CI / R2
setColWidths(wb, 1, cols = p_col, widths = 14)    # P

# Legend
legend_lines <- c(
  "Legend:",
  paste0("OR/SD = odds ratio per one standard deviation increase in the predictor ",
         "polygenic score (PRS), at p-value threshold p < ", paste(prs_threshold, collapse = ","), "."),
  "95% CI = 95% confidence interval (lower and upper bounds).",
  "Nagelkerke R2 = incremental pseudo-R2 of the PRS over the covariate-only model.",
  "P = P-value from logistic regression adjusted for age, age2, sex, genotyping batch, and PC1-10.",
  "Cohort = multi-ancestry analysis cohort (478,504 individuals with complete PRS, phenotype, and covariate data).",
  "",
  "Abbreviations:",
  paste0("PRS = polygenic score; SD = standard deviation; OR = odds ratio; CI = confidence interval; ",
         "T2D = type 2 diabetes; CAD = coronary artery disease; MI = myocardial infarction; ",
         "PAD = peripheral artery disease.")
)
legend_start <- nrow(tab) + 5
legend_style <- createStyle(fontSize = 9, wrapText = TRUE)
bold_legend  <- createStyle(textDecoration = "bold", fontSize = 9)
for (i in seq_along(legend_lines)) {
  r <- legend_start + i - 1L
  writeData(wb, 1, x = legend_lines[i], startRow = r, startCol = 1)
  sty <- if (i %in% c(1, 8)) bold_legend else legend_style
  addStyle(wb, 1, style = sty, rows = r, cols = 1)
  mergeCells(wb, 1, cols = 1:total_cols, rows = r)
}

freezePane(wb, 1, firstActiveRow = 4)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("Saved: %s\n", out_xlsx))
