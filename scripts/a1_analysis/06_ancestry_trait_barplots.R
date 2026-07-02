#!/usr/bin/env Rscript
# 06_ancestry_trait_barplots.R
# Barplots of top 10 H-matrix trait loadings per cluster for each ancestry.
#
# Usage:
#   Rscript scripts/a1_analysis/06_ancestry_trait_barplots.R

library(data.table)
library(ggplot2)

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

# Shared signed-loading helper (build_trait_loadings) lives in figure_utils.R
source(file.path(base_dir, "scripts", "a1_analysis", "figure_utils.R"))

out_dir <- "results/a1_analysis/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helper: create barplot ---
make_barplot <- function(dt, top_n, arm_prefix, nrow_facet, ncol_facet, title,
                         cluster_names = NULL) {
  dt[, cluster_label := paste0(arm_prefix, " ", cluster)]
  k_levels <- paste0(arm_prefix, " K", seq_len(length(unique(dt$cluster))))

  if (!is.null(cluster_names)) {
    dt[, cluster_label := ifelse(cluster_label %in% names(cluster_names),
                                 cluster_names[cluster_label], cluster_label)]
    k_levels <- cluster_names[k_levels]
  }
  dt[, cluster_label := factor(cluster_label, levels = k_levels)]

  top <- dt[, {
    .SD[order(-abs_loading)][seq_len(min(.N, top_n))]
  }, by = cluster_label]

  top[, direction := fifelse(loading > 0, "Increased", "Decreased")]

  top[, trait_facet := paste(trait, cluster_label, sep = "__")]
  top[, trait_ordered := reorder(trait_facet, loading), by = cluster_label]

  ggplot(top, aes(x = trait_ordered, y = loading, fill = direction)) +
    geom_col(width = 0.7) +
    scale_x_discrete(labels = function(x) sub("__.*$", "", x)) +
    scale_fill_manual(values = c("Increased" = "#F4845F", "Decreased" = "#4682B4"),
                      name = "Cluster Weight") +
    coord_flip() +
    facet_wrap(~cluster_label, nrow = nrow_facet, ncol = ncol_facet,
               scales = "free") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 24),
      strip.text = element_text(face = "bold", size = 18),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 15),
      legend.position = "bottom",
      legend.title = element_text(size = 15),
      legend.text  = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# --- Cluster name mappings (NULL = default K labels) ---
cluster_names <- list(
  META = c(
    "META K1"  = "Lpa",
    "META K2"  = "Adiponectin",
    "META K3"  = "Platelet",
    "META K4"  = "SHBG",
    "META K5"  = "Blood Pressure-Stature",
    "META K6"  = "Metabolic",
    "META K7"  = "Triglycerides-HDL",
    "META K8"  = "ALP-LDL",
    "META K9"  = "Obesity",
    "META K10" = "Glycemic"
  ),
  EUR  = c(
    "EUR K1" = "Insulin Sensitivity",
    "EUR K2" = "Lpa",
    "EUR K3" = "ApoB Triglycerides",
    "EUR K4" = "Body Composition-",
    "EUR K5" = "CRP ApoB",
    "EUR K6" = "ALP",
    "EUR K7" = "Hematologic",
    "EUR K8" = "Metabolic Syndrome",
    "EUR K9" = "Triglycerides-"
  ),
  AFR  = NULL,
  EAS  = NULL,
  SAS  = NULL
)

# --- Ancestry-specific barplots ---
ancestries <- list(
  META = list(nrow = 5, ncol = 2, width = 14, height = 22),
  EUR  = list(nrow = 3, ncol = 3, width = 18, height = 14),
  AFR  = list(nrow = 2, ncol = 2, width = 14, height = 10),
  EAS  = list(nrow = 2, ncol = 3, width = 18, height = 10),
  SAS  = list(nrow = 2, ncol = 2, width = 14, height = 10),
  # CAD-seeded run (optimal K=8); generic K labels until clusters are reviewed
  META_CAD = list(nrow = 4, ncol = 2, width = 14, height = 18)
)

for (anc in names(ancestries)) {
  cat(sprintf("Processing %s...\n", anc))
  params <- ancestries[[anc]]

  h_file <- sprintf("results/a1_analysis/%s/H_matrix_%s.tsv", anc, anc)
  dt <- build_trait_loadings(h_file, strip_source = TRUE)

  k <- length(unique(dt$cluster))
  title <- "Contributions of Traits to Clusters"

  p <- make_barplot(dt, top_n = 10, arm_prefix = anc,
                    nrow_facet = params$nrow, ncol_facet = params$ncol,
                    title = title,
                    cluster_names = cluster_names[[anc]])

  out_file <- file.path(out_dir, sprintf("barplot_%s_trait_loadings.png", anc))
  ggsave(out_file, p, width = params$width, height = params$height, dpi = 300)
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("Done.\n")
