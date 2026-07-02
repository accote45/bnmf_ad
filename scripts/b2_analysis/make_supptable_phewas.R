#!/usr/bin/env Rscript
# make_supptable_phewas.R
# Build the PheWAS supplementary table (phewas_cluster_pgs_phenome.xlsx) from the
# B2.2 results: one row per 3-character ICD-10 phenotype, with per-cluster
# Beta / SE / P columns under canonical cluster-label headers.
#
# Mirrors b2_2_manhattan_plots.R for the disease-name lookup (UKB code table,
# code_id = 19) and for remapping the cluster labels from config (so the table
# always uses the current canonical labels, overriding any stale labels in the
# CSV). Layout matches the other supplementary tables: bold title (row 1),
# multi-level header (rows 3-4), data (row 5+), wrapText legend below.
#
# Usage:
#   module load gcc/14.2.0 R/4.2.0
#   Rscript scripts/b2_analysis/make_supptable_phewas.R                 # soft (default)
#   Rscript scripts/b2_analysis/make_supptable_phewas.R --hard-assignment

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(RSQLite)
  library(stringr)
  library(yaml)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

args <- commandArgs(trailingOnly = TRUE)
hard_assignment <- "--hard-assignment" %in% args

cfg <- read_yaml("config/b2_2_config.yaml")
cluster_labels <- unlist(cfg$cluster_labels)        # canonical labels, keyed K1..K10
cluster_keys   <- paste0("K", 1:10)

# Soft (default) vs hard-assignment results
if (hard_assignment) {
  results_file <- "results/b2_2_analysis/prs_hard_assignment/phewas_results.csv"
} else {
  results_file <- "results/b2_2_analysis/phewas_results.csv"
}
out_dir  <- "results/supplementary_tables"
out_file <- file.path(out_dir, "phewas_cluster_pgs_phenome.xlsx")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chapter_labels <- c(
  E = "Endocrine",     F = "Mental",       G = "Nervous",
  H = "Eye/Ear",       I = "Circulatory",  J = "Respiratory",
  K = "Digestive",     L = "Skin",         M = "Musculoskeletal",
  N = "Genitourinary"
)

# --- ICD-10 disease names from the UKB database (same source as the manhattan) -
ukb_db <- file.path(base_dir, "../../../data/ukb/phenotype/ukb18177.db")
con    <- dbConnect(SQLite(), ukb_db)
icd10_names <- as.data.table(dbGetQuery(con,
  "SELECT value, meaning FROM code WHERE code_id = 19 AND length(value) = 3"))
dbDisconnect(con)
# Strip the leading code from the meaning (e.g. "E11 Non-insulin-dependent..." -> "...")
icd10_names[, disease_name := str_trim(str_remove(meaning, "^[A-Z]\\d+\\s+"))]
icd10_names <- icd10_names[, .(icd10_code = value, disease_name)]

# --- Load PheWAS results ------------------------------------------------------
res <- fread(results_file)
res[, cluster_label := cluster_labels[cluster]]      # override stale CSV labels
res <- merge(res, icd10_names, by = "icd10_code", all.x = TRUE)

# --- Formatting helpers -------------------------------------------------------
fmt_num <- function(x) ifelse(is.na(x), "", sprintf("%.4f", x))
fmt_p <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, sprintf("%.2e", p),
      ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.3f", p))))
}

# --- Reshape long -> wide (one row per ICD-10 code) ---------------------------
# Row order: by chapter, then ICD-10 code (matches the manhattan ordering).
codes <- unique(res[, .(icd10_code, chapter, disease_name)])
setorder(codes, chapter, icd10_code)

out <- data.table(
  `ICD-10 Code` = codes$icd10_code,
  Disease       = codes$disease_name,
  Chapter       = chapter_labels[codes$chapter]
)
for (k in cluster_keys) {
  sub <- res[cluster == k]
  m   <- sub[match(codes$icd10_code, icd10_code)]
  out[[paste0(k, "_Beta")]] <- fmt_num(m$beta)
  out[[paste0(k, "_SE")]]   <- fmt_num(m$se)
  out[[paste0(k, "_P")]]    <- fmt_p(m$p_value)
}

cat(sprintf("PheWAS table: %d ICD-10 phenotypes x %d clusters\n",
            nrow(out), length(cluster_keys)))

# --- Write Excel (multi-level header matching existing supplements) -----------
n_label_cols <- 3L
total_cols   <- n_label_cols + length(cluster_keys) * 3L

wb <- createWorkbook()
addWorksheet(wb, "PheWAS Results")

title <- "Supplementary Table X: Association of cardiometabolic cluster PGS with phenome"
addStyle(wb, 1, createStyle(textDecoration = "bold", fontSize = 12), rows = 1, cols = 1)
writeData(wb, 1, x = title, startRow = 1, startCol = 1)
mergeCells(wb, 1, cols = 1:total_cols, rows = 1)

hdr  <- createStyle(textDecoration = "bold", halign = "center", border = "bottom", wrapText = TRUE)

# Label columns span rows 3-4
for (j in seq_len(n_label_cols)) {
  lbl <- c("ICD-10 Code", "Disease", "Chapter")[j]
  writeData(wb, 1, x = lbl, startRow = 3, startCol = j)
  addStyle(wb, 1, hdr, rows = 3, cols = j)
  mergeCells(wb, 1, cols = j, rows = 3:4)
}
# Cluster headers (row 3, merged over 3 cols) + Beta/SE/P sub-headers (row 4)
for (ki in seq_along(cluster_keys)) {
  cs <- n_label_cols + (ki - 1L) * 3L + 1L
  writeData(wb, 1, x = cluster_labels[cluster_keys[ki]], startRow = 3, startCol = cs)
  addStyle(wb, 1, hdr, rows = 3, cols = cs)
  mergeCells(wb, 1, cols = cs:(cs + 2L), rows = 3)
  writeData(wb, 1, x = "Beta", startRow = 4, startCol = cs);      addStyle(wb, 1, hdr, rows = 4, cols = cs)
  writeData(wb, 1, x = "SE",   startRow = 4, startCol = cs + 1L); addStyle(wb, 1, hdr, rows = 4, cols = cs + 1L)
  writeData(wb, 1, x = "P",    startRow = 4, startCol = cs + 2L); addStyle(wb, 1, hdr, rows = 4, cols = cs + 2L)
}

writeData(wb, 1, x = out, startRow = 5, startCol = 1, colNames = FALSE)

setColWidths(wb, 1, cols = 1, widths = 12)
setColWidths(wb, 1, cols = 2, widths = 34)
setColWidths(wb, 1, cols = 3, widths = 16)
setColWidths(wb, 1, cols = (n_label_cols + 1):total_cols, widths = 9)

cluster_name_line <- paste0(
  "Cluster names: ",
  paste(sprintf("%s (%s)", cluster_labels[cluster_keys], cluster_keys), collapse = ", "), ".")
legend_lines <- c(
  "Legend:",
  "ICD-10 Code = 3-character ICD-10 diagnosis code (UK Biobank fields 41202 primary + 41204 secondary).",
  "Disease = ICD-10 code description (UK Biobank data coding 19).",
  "Chapter = ICD-10 chapter grouping.",
  "Beta / SE / P = effect size, standard error, and P-value from logistic regression of the binary phenotype on the cluster PGS, adjusted for age, age2, sex, genotyping batch, and PC1-10.",
  cluster_name_line,
  "PGS = polygenic score; SE = standard error."
)
legend_start <- nrow(out) + 7
for (i in seq_along(legend_lines)) {
  r <- legend_start + i - 1L
  writeData(wb, 1, x = legend_lines[i], startRow = r, startCol = 1)
  sty <- if (i == 1) createStyle(textDecoration = "bold", fontSize = 9) else createStyle(fontSize = 9, wrapText = TRUE)
  addStyle(wb, 1, sty, rows = r, cols = 1)
  mergeCells(wb, 1, cols = 1:total_cols, rows = r)
}

saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("Saved: %s\n", out_file))
