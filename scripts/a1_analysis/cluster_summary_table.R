#!/usr/bin/env Rscript
# cluster_summary_table.R
# Per-cluster review table for the META bNMF factorization:
#   - n variants (dominant = max-W cluster; and W > 0.05 non-trivial loading)
#   - top traits (signed H-matrix loadings, ranked by |loading|)
#   - top genes by max W weight
#   - top genes by cluster-specificity score (mean W[i,k]/sum_k W[i,] over SNPs with W[i,k]>0)
# alongside the current canonical labels.
#
# Reuses canonical helpers (no reimplementation):
#   - process_h_matrix()              from 06_ancestry_trait_barplots.R
#   - parse_gtf_genes(), map_snps_to_genes() from figure_utils.R
#   - specificity formula matches plot_gene_bars()/lollipop in 07_gene_bars_metric_comparison.R
#
# Usage:
#   Rscript scripts/a1_analysis/cluster_summary_table.R [--ancestry META] [--top-traits 5] [--top-genes 3]

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

# --- Args ---
args <- commandArgs(trailingOnly = TRUE)
ancestry    <- if ("--ancestry" %in% args) args[which(args == "--ancestry") + 1] else "META"
top_traits_n <- if ("--top-traits" %in% args) as.integer(args[which(args == "--top-traits") + 1]) else 5L
top_genes_n  <- if ("--top-genes" %in% args) as.integer(args[which(args == "--top-genes") + 1]) else 3L
W_THRESH    <- 0.05  # "not close to zero" loading threshold for the soft variant count

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

cfg <- read_yaml("config/a1_config.yaml")
results_dir <- cfg$results_dir
gtf_path    <- cfg$gtf_path
anc_dir     <- file.path(results_dir, ancestry)
h_file <- file.path(anc_dir, sprintf("H_matrix_%s.tsv", ancestry))
w_file <- file.path(anc_dir, sprintf("W_matrix_%s.tsv", ancestry))

# --- Canonical META cluster labels (from 06/07 plotting scripts) ---
meta_labels <- c(
  K1 = "Lpa",                    K2 = "Adiponectin",
  K3 = "Platelet",               K4 = "SHBG",
  K5 = "Blood Pressure-Stature", K6 = "Metabolic",
  K7 = "Triglycerides-HDL",      K8 = "ALP-LDL",
  K9 = "Obesity",                K10 = "Glycemic"
)

# --- Reuse process_h_matrix() from script 06 (eval only its function defs) ---
exprs_06 <- parse(file.path("scripts", "a1_analysis", "06_ancestry_trait_barplots.R"))
for (e in exprs_06) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, envir = globalenv())
  }
}
# --- Reuse gene-mapping helpers from figure_utils.R (function library) ---
source("scripts/a1_analysis/figure_utils.R")

# =============================================================================
# Top traits per cluster (signed loading = pos - neg, ranked by |loading|)
# =============================================================================
hdt <- process_h_matrix(h_file, strip_source = TRUE)  # cluster, trait, loading, abs_loading
hdt <- as.data.table(hdt)
top_traits <- hdt[, .SD[order(-abs_loading)][seq_len(min(.N, top_traits_n))], by = cluster]
top_traits[, dir := fifelse(loading > 0, "+", "-")]
top_traits[, label := sprintf("%s(%s)", trait, dir)]
traits_by_k <- top_traits[, .(top_traits = paste(label, collapse = ", ")), by = cluster]

# =============================================================================
# Variant composition + gene metrics from the W matrix
# =============================================================================
w <- fread(w_file)
kcols <- setdiff(names(w), "VAR_ID")
parts <- tstrsplit(w$VAR_ID, "_", fixed = TRUE)
w[, chr := as.character(parts[[1]])]
w[, pos := as.numeric(parts[[2]])]
# Per-variant total loading across clusters (denominator for specificity)
w[, .row_total := rowSums(as.matrix(w[, ..kcols]))]

# --- Variant counts (all 885 variants): dominant (max-W) and W > threshold ---
wm <- as.matrix(w[, ..kcols])
dom <- kcols[max.col(wm, ties.method = "first")]
vc_by_k <- data.table(
  cluster        = kcols,
  n_dominant     = as.integer(table(factor(dom, levels = kcols))),
  n_w_gt_thresh  = as.integer(colSums(wm > W_THRESH))
)

# --- Map SNPs to nearest protein-coding genes (carries W cols + .row_total) ---
gene_df <- parse_gtf_genes(gtf_path)
mapped  <- as.data.table(map_snps_to_genes(as.data.frame(w), gene_df))

# Top genes by MAX W weight, and by SPECIFICITY (mean W[k]/row_total over W[k]>0)
genes_by_k <- rbindlist(lapply(kcols, function(k) {
  gmax <- mapped[, .(score = max(get(k), na.rm = TRUE)), by = gene_name][
    score > 0][order(-score)][seq_len(min(.N, top_genes_n))]

  sub <- mapped[get(k) > 0 & .row_total > 0]
  sub[, .spec := get(k) / .row_total]
  gspec <- sub[, .(score = mean(.spec, na.rm = TRUE)), by = gene_name][
    order(-score)][seq_len(min(.N, top_genes_n))]

  data.table(
    cluster                = k,
    top_genes_maxweight    = paste(gmax$gene_name, collapse = ", "),
    top_genes_specificity  = paste(gspec$gene_name, collapse = ", ")
  )
}))

# =============================================================================
# Assemble + print
# =============================================================================
k_order <- paste0("K", seq_along(meta_labels))
summary_dt <- data.table(cluster = k_order, current_label = unname(meta_labels[k_order]))
summary_dt <- merge(summary_dt, vc_by_k,     by = "cluster", all.x = TRUE, sort = FALSE)
summary_dt <- merge(summary_dt, traits_by_k, by = "cluster", all.x = TRUE, sort = FALSE)
summary_dt <- merge(summary_dt, genes_by_k,  by = "cluster", all.x = TRUE, sort = FALSE)
summary_dt <- summary_dt[match(k_order, cluster)]
setnames(summary_dt, "n_w_gt_thresh", sprintf("n_w_gt_%.2f", W_THRESH))

out_tsv <- file.path(anc_dir, sprintf("cluster_summary_%s.tsv", ancestry))
fwrite(summary_dt, out_tsv, sep = "\t")

cat(sprintf("\n=== %s bNMF cluster summary (top %d traits / %d genes; trait sign = risk direction) ===\n",
            ancestry, top_traits_n, top_genes_n))
cat(sprintf("    n variants: dominant = max-W cluster; n_w_gt_%.2f = loading > %.2f\n\n", W_THRESH, W_THRESH))
for (i in seq_len(nrow(summary_dt))) {
  r <- summary_dt[i]
  cat(sprintf("%s  [%s]  (dominant=%d, W>%.2f=%d)\n",
              r$cluster, r$current_label, r$n_dominant, W_THRESH, r[[sprintf("n_w_gt_%.2f", W_THRESH)]]))
  cat(sprintf("   traits          : %s\n", r$top_traits))
  cat(sprintf("   genes (maxweight): %s\n", r$top_genes_maxweight))
  cat(sprintf("   genes (specific) : %s\n\n", r$top_genes_specificity))
}
cat(sprintf("Wrote: %s\n", out_tsv))
