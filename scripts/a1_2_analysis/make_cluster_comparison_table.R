#!/usr/bin/env Rscript
#
# make_cluster_comparison_table.R
# Build supplementary table comparing cluster assignments across 3 published studies.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(openxlsx)
})

base_dir <- getwd()

# --- 1. Read published cluster assignment data --------------------------------

cat("Reading published cluster data...\n")

# Suzuki 2024: hard clustering (ST6)
suz_raw <- as.data.table(read_excel(
  file.path(base_dir, "data/published_clusters/Suzuki2024_Nature_SuppTables.xlsx"),
  sheet = "ST6", skip = 2
))
suz_raw <- suz_raw[!is.na(`Index SNV`)]
suz <- suz_raw[, .(
  locus_suz = Locus,
  rsid_suz  = `Index SNV`,
  cluster_suz = `Cluster assignment`
)]
cat(sprintf("  Suzuki: %d variants, %d clusters\n", nrow(suz), uniqueN(suz$cluster_suz)))

# Smith 2024: soft clustering (S5) — assign to max-weight cluster
smith_raw <- as.data.table(read_excel(
  file.path(base_dir, "data/published_clusters/Smith2024_NatMed_SuppTables.xlsx"),
  sheet = "S5", skip = 3
))
# Cluster weight columns are 4-17 (first set of 14 cluster columns)
cluster_cols_smith <- names(smith_raw)[4:17]
# Clean up column names (remove ...N suffix)
clean_names <- gsub("\\.\\.\\.[0-9]+$", "", cluster_cols_smith)
smith_weights <- as.matrix(smith_raw[, ..cluster_cols_smith])
mode(smith_weights) <- "numeric"
max_cluster <- clean_names[max.col(smith_weights, ties.method = "first")]
smith <- data.table(
  locus_smith = smith_raw$locus,
  rsid_smith  = smith_raw$rsID,
  cluster_smith = max_cluster
)
smith <- smith[!is.na(rsid_smith)]
cat(sprintf("  Smith: %d variants, %d clusters\n", nrow(smith), uniqueN(smith$cluster_smith)))

# Pascat 2026: hard clustering (Supplementary Data 3)
pas_raw <- as.data.table(read_excel(
  file.path(base_dir, "data/published_clusters/Pascat2026_NatComms_SuppTables.xlsx"),
  sheet = "Supplementary Data 3", col_names = FALSE, skip = 3
))
pas <- data.table(
  rsid_pas    = pas_raw[[1]],
  locus_pas   = pas_raw[[5]],
  cluster_pas = pas_raw[[4]],
  # Cross-references already in Pascat SD3
  rsid_suz_xref    = pas_raw[[6]],
  ld_r2_suz        = as.numeric(pas_raw[[7]]),
  cluster_suz_xref = pas_raw[[8]],
  rsid_smith_xref    = pas_raw[[9]],
  ld_r2_smith        = as.numeric(pas_raw[[10]]),
  cluster_smith_xref = pas_raw[[11]]
)
pas <- pas[!is.na(rsid_pas)]
cat(sprintf("  Pascat: %d variants, %d clusters\n", nrow(pas), uniqueN(pas$cluster_pas)))

# --- 2. Build unified locus table ---------------------------------------------

cat("\nBuilding unified locus table...\n")

# Strategy: split multi-gene loci into individual genes, then match studies
# by any shared gene. This handles "GCKR" matching "DTNB, KIF3C, GCKR".

split_genes <- function(locus_str) {
  genes <- trimws(unlist(strsplit(locus_str, "[,;/]")))
  toupper(genes[nzchar(genes)])
}

# Build gene→variant lookup for each study
build_gene_map <- function(dt, locus_col, id_col) {
  rows <- list()
  for (i in seq_len(nrow(dt))) {
    genes <- split_genes(dt[[locus_col]][i])
    for (g in genes) rows[[length(rows) + 1]] <- data.table(gene = g, idx = i)
  }
  rbindlist(rows)
}

suz_gmap <- build_gene_map(suz, "locus_suz", "rsid_suz")
smith_gmap <- build_gene_map(smith, "locus_smith", "rsid_smith")
pas_gmap <- build_gene_map(pas, "locus_pas", "rsid_pas")

# Union of all gene symbols
all_genes <- sort(unique(c(suz_gmap$gene, smith_gmap$gene, pas_gmap$gene)))

out_rows <- list()
used_suz <- logical(nrow(suz))
used_smith <- logical(nrow(smith))
used_pas <- logical(nrow(pas))

for (g in all_genes) {
  si <- suz_gmap[gene == g, idx]
  mi <- smith_gmap[gene == g, idx]
  pi <- pas_gmap[gene == g, idx]
  # Skip if all already used
  si <- si[!used_suz[si]]; mi <- mi[!used_smith[mi]]; pi <- pi[!used_pas[pi]]
  if (length(si) == 0 && length(mi) == 0 && length(pi) == 0) next
  # Take first unused variant from each study
  si <- if (length(si)) si[1] else NA_integer_
  mi <- if (length(mi)) mi[1] else NA_integer_
  pi <- if (length(pi)) pi[1] else NA_integer_
  # Mark used
  if (!is.na(si)) used_suz[si] <- TRUE
  if (!is.na(mi)) used_smith[mi] <- TRUE
  if (!is.na(pi)) used_pas[pi] <- TRUE
  # Best locus name
  locus_name <- if (!is.na(si)) suz$locus_suz[si] else
                if (!is.na(mi)) smith$locus_smith[mi] else
                pas$locus_pas[pi]
  out_rows[[length(out_rows) + 1]] <- data.table(
    Locus = locus_name,
    rsid_pas    = if (!is.na(pi)) pas$rsid_pas[pi] else "",
    cluster_pas = if (!is.na(pi)) pas$cluster_pas[pi] else "",
    rsid_smith    = if (!is.na(mi)) smith$rsid_smith[mi] else "",
    cluster_smith = if (!is.na(mi)) smith$cluster_smith[mi] else "",
    rsid_suz    = if (!is.na(si)) suz$rsid_suz[si] else "",
    cluster_suz = if (!is.na(si)) suz$cluster_suz[si] else ""
  )
}

out <- rbindlist(out_rows)
out <- out[order(Locus)]

cat(sprintf("  Total unique loci: %d\n", nrow(out)))
cat(sprintf("  In all 3 studies: %d\n",
    sum(out$rsid_pas != "" & out$rsid_smith != "" & out$rsid_suz != "")))
cat(sprintf("  In 2 studies: %d\n",
    sum((out$rsid_pas != "") + (out$rsid_smith != "") + (out$rsid_suz != "") == 2)))
cat(sprintf("  Pascat only: %d\n",
    sum(out$rsid_pas != "" & out$rsid_smith == "" & out$rsid_suz == "")))
cat(sprintf("  Smith only: %d\n",
    sum(out$rsid_pas == "" & out$rsid_smith != "" & out$rsid_suz == "")))
cat(sprintf("  Suzuki only: %d\n",
    sum(out$rsid_pas == "" & out$rsid_smith == "" & out$rsid_suz != "")))

# --- 3. Write Excel -----------------------------------------------------------

cat("\nWriting Excel file...\n")

wb <- createWorkbook()
addWorksheet(wb, "Cluster Comparison")

# Row 1: Table header (bold)
title_style <- createStyle(textDecoration = "bold", fontSize = 11)
writeData(wb, 1, x = "Supplementary Table X: Comparison of cluster assignments across published T2D genetic clustering studies",
          startRow = 1, startCol = 1)
addStyle(wb, 1, title_style, rows = 1, cols = 1)

# Row 2: blank (skip)

# Row 3: Study-level header (multi-level top)
study_style <- createStyle(textDecoration = "bold", halign = "center", border = "bottom")
writeData(wb, 1, x = "Pascat et al. 2026", startRow = 3, startCol = 2)
addStyle(wb, 1, study_style, rows = 3, cols = 2)
mergeCells(wb, 1, cols = 2:3, rows = 3)

writeData(wb, 1, x = "Smith et al. 2024", startRow = 3, startCol = 4)
addStyle(wb, 1, study_style, rows = 3, cols = 4)
mergeCells(wb, 1, cols = 4:5, rows = 3)

writeData(wb, 1, x = "Suzuki et al. 2024", startRow = 3, startCol = 6)
addStyle(wb, 1, study_style, rows = 3, cols = 6)
mergeCells(wb, 1, cols = 6:7, rows = 3)

# Row 4: Sub-column headers
sub_style <- createStyle(textDecoration = "bold", border = "bottom")
headers <- c("Locus", "Index SNV", "Cluster", "Index SNV", "Cluster", "Index SNV", "Cluster")
for (j in seq_along(headers)) {
  writeData(wb, 1, x = headers[j], startRow = 4, startCol = j)
  addStyle(wb, 1, sub_style, rows = 4, cols = j)
}

# Row 5+: Data
writeData(wb, 1, x = out, startRow = 5, startCol = 1, colNames = FALSE)

# Column widths
setColWidths(wb, 1, cols = 1, widths = 22)
setColWidths(wb, 1, cols = c(2, 4, 6), widths = 16)
setColWidths(wb, 1, cols = c(3, 5, 7), widths = 28)

# Legend below data
legend_row <- nrow(out) + 7
legend_style <- createStyle(fontSize = 9, wrapText = TRUE)
legend_lines <- c(
  "SNV, single nucleotide variant; T2D, type 2 diabetes; CAD, coronary artery disease.",
  "Locus names are matched across studies by gene symbol. Where a study does not include a variant at a given locus, cells are left blank.",
  "Pascat et al. 2026: 5 hard clusters from hierarchical clustering of T2D and blood pressure GWAS loci (1,303 SNVs).",
  "Smith et al. 2024: 14 soft clusters from Bayesian non-negative matrix factorization of T2D GWAS loci (650 SNVs). Each variant is assigned to its maximum-weight cluster.",
  "Suzuki et al. 2024: 8 hard clusters from hierarchical clustering of T2D GWAS loci (1,289 SNVs)."
)
for (i in seq_along(legend_lines)) {
  writeData(wb, 1, x = legend_lines[i], startRow = legend_row + i - 1, startCol = 1)
  addStyle(wb, 1, legend_style, rows = legend_row + i - 1, cols = 1)
}

out_path <- file.path(base_dir, "results/supplementary_tables/comparison_cluster_assignments.xlsx")
saveWorkbook(wb, out_path, overwrite = TRUE)
cat(sprintf("Saved: %s\n", out_path))
cat("Done.\n")
