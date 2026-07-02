#!/usr/bin/env Rscript
# report_cluster_relabel.R
# Read-only reporting helper: prints, per cluster, the top traits + top genes and
# the current suggested label so the labels can be reviewed/relabeled.
#
# Reuses existing logic (no reimplementation):
#   - process_h_matrix()  from 06_ancestry_trait_barplots.R (trait loadings)
#   - parse_gtf_genes(), map_snps_to_genes() from figure_utils.R (gene mapping)
#
# Usage: Rscript scripts/a1_analysis/report_cluster_relabel.R

suppressMessages({
  library(tidyverse)
  library(data.table)
  library(yaml)
})

project_root <- getwd()
cfg <- read_yaml("config/a1_config.yaml")
gtf_path    <- cfg$gtf_path
results_dir <- cfg$results_dir

# --- Source helpers (gene mapping) ---
source(file.path(project_root, "scripts", "a1_analysis", "figure_utils.R"))

# --- Pull function definitions ONLY from script 06 (skip its main loop) ---
exprs_06 <- parse(file.path(project_root, "scripts", "a1_analysis",
                            "06_ancestry_trait_barplots.R"))
for (e in exprs_06) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, envir = globalenv())
  }
}
# Recover the current cluster_names mapping from script 06 as well
for (e in exprs_06) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("cluster_names"))) {
    eval(e, envir = globalenv())
  }
}

TOP_TRAITS <- 10
TOP_GENES  <- 5

gene_df <- parse_gtf_genes(gtf_path)

report_ancestry <- function(anc) {
  cat(sprintf("\n\n================  %s  ================\n", anc))

  # --- Top traits per cluster (process_h_matrix) ---
  h_file <- sprintf("%s/%s/H_matrix_%s.tsv", results_dir, anc, anc)
  hdt <- process_h_matrix(h_file, strip_source = TRUE)
  hdt[, dir := fifelse(loading > 0, "+", "-")]

  # --- Top genes per cluster (max W-loading, gene-body overlap) ---
  w_path <- sprintf("%s/%s/W_matrix_%s.tsv", results_dir, anc, anc)
  w_df <- read_tsv(w_path, show_col_types = FALSE) %>%
    mutate(chr = str_replace(VAR_ID, "^([^_]+)_.*", "\\1"),
           pos = as.numeric(str_replace(VAR_ID, "^[^_]+_([^_]+)_.*", "\\1")))
  cluster_cols <- setdiff(colnames(w_df), c("VAR_ID", "chr", "pos"))
  mapped <- map_snps_to_genes(w_df, gene_df)

  cn <- cluster_names[[anc]]

  for (kcol in cluster_cols) {
    facet_key <- paste(anc, kcol)
    cur_label <- if (!is.null(cn) && facet_key %in% names(cn)) cn[[facet_key]] else "(none)"

    traits <- hdt[cluster == kcol][order(-abs_loading)][seq_len(min(.N, TOP_TRAITS))]
    trait_str <- paste0(traits$trait, "(", traits$dir, ")", collapse = ", ")

    genes <- mapped %>%
      group_by(gene_name) %>%
      summarise(score = max(.data[[kcol]], na.rm = TRUE), .groups = "drop") %>%
      filter(score > 0) %>%
      slice_max(score, n = TOP_GENES, with_ties = FALSE)
    gene_str <- paste(genes$gene_name, collapse = ", ")

    cat(sprintf("\n%s  [current label: %s]\n", facet_key, cur_label))
    cat(sprintf("  Top traits: %s\n", trait_str))
    cat(sprintf("  Top genes : %s\n", gene_str))
  }
}

for (anc in c("META", "EUR")) report_ancestry(anc)
cat("\n")
