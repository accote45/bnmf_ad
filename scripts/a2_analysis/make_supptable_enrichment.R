#!/usr/bin/env Rscript
#
# Build the A2 enrichment supplementary tables (.xlsx) from the current results:
#   - gtex_tissue_enrichment.xlsx            (GTEx tissue Abs/Rel expression p-values)
#   - tabula_sapiens_single_cell_enrichment.xlsx (TS weighted specificity score / p-value)
#   - catlas_single_cell_enrichment.xlsx     (catlas matched-null Firth enrichment beta / p-value)
#
# Matches the layout of the existing tables: bold title (row 1, size 12), plain
# header (row 2), data rows, a blank row, then a wrapText legend block (size 10).
# Cluster labels use the canonical new META labels.
#
# Usage:
#   module load gcc/14.2.0 R/4.2.0
#   Rscript scripts/a2_analysis/make_supptable_enrichment.R

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
res_dir  <- file.path(base_dir, "results/a2_analysis")
out_dir  <- file.path(base_dir, "results/supplementary_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

EMDASH <- "—"

# Canonical new META cluster labels (keep in sync with the A2 plotting scripts)
cluster_labels <- c(
  K1 = "Lpa", K2 = "Adiponectin", K3 = "Platelet",
  K4 = "SHBG", K5 = "Blood Pressure-Stature", K6 = "Metabolic",
  K7 = "Triglycerides-HDL", K8 = "ALP-LDL", K9 = "Obesity", K10 = "Glycemic"
)
cluster_display <- function(k) sprintf("%s (%s)", k, cluster_labels[k])
cluster_num     <- function(k) as.integer(sub("K", "", k))

# Display formatting (mirrors the heatmap scripts / existing tables)
fmt_tissue_camel <- function(x) gsub("([a-z])([A-Z])", "\\1 \\2", x)  # GTEx camelCase
fmt_underscore   <- function(x) gsub("_", " ", x)                      # TS / catlas

# --- Generic writer: one styled sheet matching the existing supp-table layout ---
build_and_save <- function(out_path, sheet_name, title, df, legend_lines, widths) {
  wb <- createWorkbook()
  addWorksheet(wb, sheet_name)

  # Title (row 1, bold size 12)
  writeData(wb, sheet_name, title, startRow = 1, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet_name, createStyle(textDecoration = "bold", fontSize = 12),
           rows = 1, cols = 1)

  # Header (row 2) + data (row 3+)
  writeData(wb, sheet_name, df, startRow = 2, colNames = TRUE)

  # Legend block: one blank row after the data, then the lines (size 10, wrapped)
  legend_start <- 2 + nrow(df) + 2
  for (i in seq_along(legend_lines)) {
    writeData(wb, sheet_name, legend_lines[i],
              startRow = legend_start + i - 1, startCol = 1, colNames = FALSE)
  }
  addStyle(wb, sheet_name, createStyle(fontSize = 10, wrapText = TRUE),
           rows = legend_start:(legend_start + length(legend_lines) - 1), cols = 1)

  setColWidths(wb, sheet_name, cols = seq_along(widths), widths = widths)
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("  Saved %s (%d data rows)", basename(out_path), nrow(df)))
}

# ============================ 1. GTEx ========================================
message("Building GTEx tissue enrichment table...")
g <- fread(file.path(res_dir, "gtex_ttest_results_nearest.csv"))
gw <- dcast(g, cluster + tissue ~ expr_type, value.var = "p_value")
gw[, `:=`(Cluster = cluster_display(cluster),
          Tissue  = fmt_tissue_camel(tissue),
          knum    = cluster_num(cluster))]
setorder(gw, knum, Tissue)
gtex_df <- gw[, .(Cluster,
                  Tissue,
                  `Relative Expression P-value` = rel,
                  `Absolute Expression P-value` = abs)]
gtex_legend <- c(
  "Legend:",
  "Cluster: bNMF cluster ID and label from the META multi-ancestry decomposition (K = 10).",
  "Tissue: GTEx v8 tissue name (46 tissues total).",
  paste0("Relative Expression P-value: one-sided Welch t-test p-value comparing mean ",
         "tissue-specificity scores of cluster-assigned genes (within ±50 kb of each ",
         "bNMF variant) versus all other protein-coding genes. Tests whether cluster genes ",
         "are more specifically expressed in that tissue than expected."),
  paste0("Absolute Expression P-value: one-sided Welch t-test p-value comparing mean ",
         "log2(TPM + 1) expression of cluster-assigned genes versus other expressed ",
         "protein-coding genes (expressed-gene background). Tests whether cluster genes are ",
         "more highly expressed in that tissue than expected."),
  "TPM: Transcripts Per Million.",
  "bNMF: Bayesian non-negative matrix factorization."
)
build_and_save(
  file.path(out_dir, "gtex_tissue_enrichment.xlsx"),
  "GTEx Analysis",
  "Supplementary Table X: GTEx tissue expression enrichment analysis for cardiometabolic clusters",
  as.data.frame(gtex_df), gtex_legend, widths = c(40, 38, 30, 30)
)

# ----- Shared builder for the single-cell WSS tables (TS + catlas) -----------
build_wss_table <- function(csv_path, out_name, sheet_name, title, legend_lines,
                            tissue_formatter) {
  d <- fread(csv_path)
  d[, `:=`(Cluster   = cluster_display(cluster),
           tdisp     = tissue_formatter(tissue),
           cdisp     = fmt_underscore(cell_type),
           knum      = cluster_num(cluster))]
  d[, `Tissue — Cell Type` := paste(tdisp, cdisp, sep = paste0(" ", EMDASH, " "))]
  setorder(d, knum, tdisp, cdisp)
  out_df <- d[, .(Cluster,
                  `Tissue — Cell Type`,
                  `Weighted Specificity Score` = wss_observed,
                  `P-value` = p_value)]
  build_and_save(file.path(out_dir, out_name), sheet_name, title,
                 as.data.frame(out_df), legend_lines, widths = c(40, 45, 28, 18))
}

# ============================ 2. Tabula Sapiens ==============================
message("Building Tabula Sapiens single-cell enrichment table...")
ts_legend <- c(
  "Legend:",
  "Cluster: bNMF cluster ID and label from the META multi-ancestry decomposition (K = 10).",
  paste0("Tissue — Cell Type: Tabula Sapiens tissue and cell type annotation. Tissues ",
         "include Blood, Fat, Heart, Large Intestine, Liver, Lung, Lymph Node, Mammary, ",
         "Prostate, Skin, Small Intestine, Spleen, Thymus, and Vasculature ",
         "(14 tissues, 307 cell types total)."),
  paste0("Weighted Specificity Score (WSS): For each cluster k and cell type c, ",
         "WSS(k, c) = sum of W[i, k] x specificity[gene_i, c] across all variants i mapped ",
         "to genes present in the tissue. W[i, k] is the continuous bNMF loading weight for ",
         "variant i in cluster k. Specificity scores represent cell-type-specific expression ",
         "relative to other cell types within the same tissue."),
  paste0("P-value: Permutation test p-value (10,000 permutations). The null distribution was ",
         "generated by randomly reassigning gene specificity values to variant positions and ",
         "recomputing the WSS. The p-value is the proportion of permuted WSS values greater ",
         "than or equal to the observed WSS. Per-cluster FDR correction was applied using the ",
         "Benjamini-Hochberg method."),
  "bNMF: Bayesian non-negative matrix factorization.",
  "FDR: False discovery rate.",
  "WSS: Weighted Specificity Score."
)
build_wss_table(
  file.path(res_dir, "tabula_sapiens/tabula_sapiens_wss_results.csv"),
  "tabula_sapiens_single_cell_enrichment.xlsx", "Tabula Sapiens Analysis",
  "Supplementary Table X: Tabula Sapiens single-cell analysis for cardiometabolic clusters",
  ts_legend, tissue_formatter = fmt_underscore
)

# ============================ 3. catlas =====================================
# CATLAS uses the SNP-level matched-null Firth enrichment (the WSS permutation
# columns are replaced by the Firth coefficient beta + its penalized-LR p-value).
# This result set is structurally different from WSS: it covers only the
# Firth-testable set (clusters with >= 10 matched lead SNPs; K1/Lpa excluded),
# and each row carries a real beta + p-value, so it is built directly here
# rather than through build_wss_table().
message("Building catlas single-cell enrichment table (Firth)...")
cat <- fread(file.path(res_dir,
                       "catlas/catlas_firth/catlas_firth_enrichment_results.csv"))
cat[, `:=`(Cluster = cluster_display(cluster),
           knum    = cluster_num(cluster))]
# `tissue` (Pancreas/Esophagus) and `label` (nice cell-type name) are pre-formatted upstream.
cat[, `Tissue — Cell Type` := paste(tissue, label, sep = paste0(" ", EMDASH, " "))]
setorder(cat, knum, tissue, label)
catlas_df <- cat[, .(Cluster,
                     `Tissue — Cell Type`,
                     `Firth Enrichment Beta` = beta,
                     `P-value` = p_value)]

catlas_legend <- c(
  "Legend:",
  paste0("Cluster: bNMF cluster ID and label from the META multi-ancestry decomposition. ",
         "Of the 10 clusters, only those with >= 10 matched lead SNPs are testable; the Lpa ",
         "cluster (K1) is excluded for having too few lead SNPs."),
  paste0("Tissue — Cell Type: CATLAS single-cell ATAC-seq tissue and cell type ",
         "annotation. Tissues include Pancreas and Esophagus (2 tissues, 16 cell types total)."),
  paste0("Firth Enrichment Beta: coefficient (theta) for the cell-type ATAC peak-overlap ",
         "indicator from a matched-null Firth bias-reduced logistic regression, ",
         "logit P(lead) = a0 + a_CDS*CDS + a_5UTR*5'UTR + a_3UTR*3'UTR + theta*(peak overlap). ",
         "Lead SNPs (Y=1) are the bNMF variants assigned to the cluster by maximum loading ",
         "weight; matched null SNPs (Y=0) are 1000 Genomes EUR variants within ±50 kb of a ",
         "cluster lead and not in LD (r2 < 0.05) with any cluster lead. Genic covariates ",
         "(CDS, 5'UTR, 3'UTR) are from GENCODE v19. The beta is a log-odds ratio: a positive ",
         "value means lead SNPs are more likely than matched null SNPs to fall within that ",
         "cell type's open-chromatin peaks."),
  paste0("P-value: penalized likelihood-ratio test p-value for the peak-overlap coefficient ",
         "(theta) from the Firth regression. Significance was assessed with a per-cluster ",
         "Bonferroni correction over the cell types tested."),
  "bNMF: Bayesian non-negative matrix factorization.",
  "LD: linkage disequilibrium.",
  "Firth: Firth's bias-reduced (penalized-likelihood) logistic regression.",
  "cCRE: candidate cis-regulatory element."
)
build_and_save(
  file.path(out_dir, "catlas_single_cell_enrichment.xlsx"),
  "CATLAS Analysis",
  "Supplementary Table X: CATLAS single-cell analysis for cardiometabolic clusters",
  as.data.frame(catlas_df), catlas_legend, widths = c(40, 45, 28, 18)
)

message("\nDone.")
