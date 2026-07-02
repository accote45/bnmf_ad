#!/usr/bin/env Rscript
# make_supptable_cluster_matrices.R
# Build two supplementary tables for the META bNMF clusters:
#   1. cluster_weights_meta.xlsx       — W-matrix variant weights (variants x clusters)
#   2. cluster_trait_loadings_meta.xlsx — H-matrix net trait loadings (traits x clusters)
#
# Cluster columns are labelled with the canonical META cluster names. Variant
# annotation (rsID, Locus) is reused from the cad_t2d_input_variants table so it
# stays identical to that supplement (no re-derivation of gene mapping).
#
# Usage:
#   Rscript scripts/a1_analysis/make_supptable_cluster_matrices.R

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(readxl)
  library(yaml)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

cfg <- read_yaml("config/a1_config.yaml")
results_dir <- cfg$results_dir
out_dir <- "results/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Canonical META cluster labels (K1..K10), from 06/07/09 ---
meta_labels <- c(
  K1 = "Lpa",                    K2 = "Adiponectin",
  K3 = "Platelet",              K4 = "SHBG",
  K5 = "Blood Pressure-Stature", K6 = "Metabolic",
  K7 = "Triglycerides-HDL",      K8 = "ALP-LDL",
  K9 = "Obesity",               K10 = "Glycemic"
)
cluster_name_line <- paste0(
  "Cluster names: ",
  paste(sprintf("%s (%s)", meta_labels, names(meta_labels)), collapse = ", "), ".")

# --- Reuse process_h_matrix() from script 06 (skip its main loop) ---
exprs_06 <- parse(file.path("scripts", "a1_analysis", "06_ancestry_trait_barplots.R"))
for (e in exprs_06) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, envir = globalenv())
  }
}

# --- Reference GWAS study lists (dynamic from config subset) ---
ref_gwas  <- cfg$ref_gwas$META
cad_files <- ref_gwas[grep("^CAD_", names(ref_gwas))]
t2d_files <- ref_gwas[grep("^T2D_", names(ref_gwas))]
study_etal <- function(gwas_named) {
  s <- sub(".*?([A-Za-z]+)([0-9]{4}).*", "\\1 et al. \\2", names(gwas_named))
  sort(unique(s))
}
cad_src <- paste(study_etal(cad_files), collapse = ", ")
t2d_src <- paste(study_etal(t2d_files), collapse = ", ")

# =============================================================================
# Helper: write a labelled supplementary table (title row 1, headers row 3,
# data row 4+, legend below) matching the existing supplement format.
# =============================================================================
write_supp <- function(out, title, sheet, legend_lines, col_widths, num_cols, file) {
  wb <- createWorkbook()
  addWorksheet(wb, sheet)

  title_style <- createStyle(textDecoration = "bold", fontSize = 12)
  writeData(wb, 1, x = title, startRow = 1, startCol = 1)
  addStyle(wb, 1, title_style, rows = 1, cols = 1)
  mergeCells(wb, 1, cols = 1:ncol(out), rows = 1)

  hdr_style <- createStyle(textDecoration = "bold", border = "Bottom",
                           halign = "center", wrapText = TRUE)
  writeData(wb, 1, x = out, startRow = 3, headerStyle = hdr_style)

  if (length(num_cols) > 0) {
    num_style <- createStyle(numFmt = "0.0000")
    addStyle(wb, 1, num_style, rows = 4:(nrow(out) + 3), cols = num_cols,
             gridExpand = TRUE)
  }

  legend_start <- nrow(out) + 5
  legend_style <- createStyle(fontSize = 10, wrapText = TRUE)
  bold_legend  <- createStyle(textDecoration = "bold", fontSize = 10)
  for (i in seq_along(legend_lines)) {
    writeData(wb, 1, x = legend_lines[i], startRow = legend_start + i - 1, startCol = 1)
    sty <- if (i == 1) bold_legend else legend_style
    addStyle(wb, 1, sty, rows = legend_start + i - 1, cols = 1)
    mergeCells(wb, 1, cols = 1:ncol(out), rows = legend_start + i - 1)
  }

  setColWidths(wb, 1, cols = seq_along(col_widths), widths = col_widths)
  saveWorkbook(wb, file, overwrite = TRUE)
  cat(sprintf("Saved: %s  (%d rows x %d cols)\n", file, nrow(out), ncol(out)))
}

# =============================================================================
# TABLE 1: cluster_weights_meta — W-matrix weights per variant
# =============================================================================
cat("=== Building cluster_weights_meta ===\n")
W <- fread(file.path(results_dir, "META", "W_matrix_META.tsv"))
kcols <- setdiff(names(W), "VAR_ID")             # K1..K10 in matrix order
stopifnot(all(kcols %in% names(meta_labels)))

# Annotation (rsID, Locus) reused from the cad_t2d_input_variants supplement
iv_path <- file.path(out_dir, "cad_t2d_input_variants.xlsx")
if (!file.exists(iv_path)) {
  stop("cad_t2d_input_variants.xlsx not found — run make_supptable_input_variants.R first")
}
# Column headers are on row 3 (title row 1, blank row 2), so skip the first 2 rows
iv <- as.data.table(read_excel(iv_path, sheet = "Input Variants", skip = 2))
iv <- iv[!is.na(`Variant ID (hg19)`)]
annot <- iv[, .(VAR_ID = `Variant ID (hg19)`, rsID, Locus)]

w_tab <- merge(W, annot, by = "VAR_ID", all.x = TRUE, sort = FALSE)
# Order to match the input-variants table (by chr, pos)
w_tab[, c("chr", "pos") := .(as.integer(tstrsplit(VAR_ID, "_")[[1]]),
                             as.integer(tstrsplit(VAR_ID, "_")[[2]]))]
setorder(w_tab, chr, pos)

# Round weights and assemble in label order
for (k in kcols) w_tab[[k]] <- round(w_tab[[k]], 4)
weights_out <- w_tab[, c(list(`Variant ID (hg19)` = VAR_ID, rsID = rsID, Locus = Locus),
                         setNames(lapply(kcols, function(k) w_tab[[k]]), meta_labels[kcols]))]
weights_out <- as.data.table(weights_out)

weights_legend <- c(
  "Legend:",
  "Variant ID (hg19) — Genomic coordinates in GRCh37/hg19 format: chr_position_ref_alt.",
  "rsID — dbSNP reference SNP identifier; left blank where no rsID is available.",
  "Locus — Nearest protein-coding gene (Ensembl GRCh37.75 GTF annotation).",
  "Cluster columns — bNMF W-matrix cluster weights representing each variant's contribution to each cluster.",
  cluster_name_line,
  "Variant universe derived from genome-wide significant loci in multi-ancestry CAD and T2D GWAS.",
  sprintf("CAD GWAS sources: %s.", cad_src),
  sprintf("T2D GWAS sources: %s.", t2d_src),
  "bNMF — Bayesian Non-negative Matrix Factorization."
)
write_supp(
  out = weights_out,
  title = "Supplementary Table X: Weights for CAD and T2D loci and associated traits in the multi-ancestry clusters",
  sheet = "Cluster Weights",
  legend_lines = weights_legend,
  col_widths = c(22, 14, 14, rep(16, length(kcols))),
  num_cols = 4:(3 + length(kcols)),
  file = file.path(out_dir, "cluster_weights_meta.xlsx")
)

# =============================================================================
# TABLE 2: cluster_trait_loadings_meta — H-matrix net loadings per trait
# =============================================================================
cat("\n=== Building cluster_trait_loadings_meta ===\n")
h_file <- file.path(results_dir, "META", "H_matrix_META.tsv")
hdt <- process_h_matrix(h_file, strip_source = TRUE)   # long: cluster, trait, loading
# Pivot to traits x clusters (net loading pos - neg, per-cluster source dedup)
hw <- dcast(hdt, trait ~ cluster, value.var = "loading")
kcols_h <- intersect(names(meta_labels), names(hw))     # K1..K10 present, in order
setorder(hw, trait)
for (k in kcols_h) hw[[k]] <- round(hw[[k]], 4)

loadings_out <- hw[, c(list(Trait = trait),
                       setNames(lapply(kcols_h, function(k) hw[[k]]), meta_labels[kcols_h]))]
loadings_out <- as.data.table(loadings_out)

loadings_legend <- c(
  "Legend:",
  "Trait — Base trait name after removing GWAS source suffix. When multiple GWAS sources were available for a trait, the source with the largest absolute loading per cluster was retained.",
  "Cluster columns — bNMF H-matrix trait loadings representing each trait's contribution to each cluster. Values are net loadings computed as the positive component minus the negative component (pos - neg).",
  "Positive values indicate the trait loads in the risk-increasing direction for that cluster; negative values indicate the risk-decreasing direction.",
  cluster_name_line,
  "bNMF — Bayesian Non-negative Matrix Factorization.",
  "Abbreviations: HDL, high-density lipoprotein; ApoA, apolipoprotein A; ApoB, apolipoprotein B; Lpa, lipoprotein(a); BMI, body mass index; BMR, basal metabolic rate; VAT, visceral adipose tissue; ASAT, abdominal subcutaneous adipose tissue; GFAT, gluteofemoral adipose tissue; SBP, systolic blood pressure; DBP, diastolic blood pressure; HR, heart rate; HbA1c, glycated hemoglobin; ALT, alanine aminotransferase; AST, aspartate aminotransferase; GGT, gamma-glutamyl transferase; ALP, alkaline phosphatase; CRP, C-reactive protein; uACR, urinary albumin-to-creatinine ratio; IGF1, insulin-like growth factor 1; SHBG, sex hormone-binding globulin; VitD, vitamin D; RBC, red blood cell count; MCV, mean corpuscular volume; PltCount, platelet count; MeanPltVol, mean platelet volume; PltDistWidth, platelet distribution width; RBCDistWidth, red blood cell distribution width."
)
write_supp(
  out = loadings_out,
  title = "Supplementary Table X: Net trait loadings for the multi-ancestry clusters",
  sheet = "Trait Loadings",
  legend_lines = loadings_legend,
  col_widths = c(20, rep(16, length(kcols_h))),
  num_cols = 2:(1 + length(kcols_h)),
  file = file.path(out_dir, "cluster_trait_loadings_meta.xlsx")
)

cat("\nDone.\n")
