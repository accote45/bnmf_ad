#!/usr/bin/env Rscript
# 08_combined_trait_gene_facets.R
# Combined per-cluster figure: each cluster appears ONCE as a composite holding
# its top-trait barplot (left, from the H matrix) and its gene scatter
# (right, max-weight vs specificity, from the W matrix). The 10 cluster
# composites are tiled in a 2-cluster-column x 5-row grid, replacing the two
# separate 10-facet panels (barplot_<anc>_trait_loadings + barplot_<anc>_gene_scatter)
# and their duplicated cluster headers.
#
# Usage:
#   Rscript scripts/a1_analysis/08_combined_trait_gene_facets.R \
#     --config config/a1_config.yaml --ancestry META

library(tidyverse)
library(data.table)
library(cowplot)
library(ggrepel)
library(yaml)

# --- Parse CLI args ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/a1_config.yaml"
ancestry    <- "META"
if ("--config" %in% args)   config_path <- args[which(args == "--config") + 1]
if ("--ancestry" %in% args) ancestry    <- args[which(args == "--ancestry") + 1]

cfg <- read_yaml(config_path)
project_root <- getwd()
results_dir  <- file.path(project_root, cfg$results_dir)
gtf_path     <- if (startsWith(cfg$gtf_path, "/")) cfg$gtf_path else file.path(project_root, cfg$gtf_path)
figures_dir  <- file.path(results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# --- Shared data-prep helpers (build_trait_loadings, build_gene_scatter_df) ---
source(file.path(project_root, "scripts", "a1_analysis", "figure_utils.R"))

# --- Cluster labels (per-script copy; canonical labels from 06/07) ---
cluster_labels <- list(
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
  EUR = c(
    "EUR K1" = "Insulin Sensitivity",
    "EUR K2" = "Lpa",
    "EUR K3" = "ApoB Triglycerides",
    "EUR K4" = "Body Composition-",
    "EUR K5" = "CRP ApoB",
    "EUR K6" = "ALP",
    "EUR K7" = "Hematologic",
    "EUR K8" = "Metabolic Syndrome",
    "EUR K9" = "Triglycerides-"
  )
)

# --- Layout: 2 cluster columns; figure size per ancestry ---
layout_config <- list(
  META = list(ncol = 2, width = 26, height = 24)
)
default_layout <- list(ncol = 2, width = 24, height = 22)

# Shared colors (identical to standalone panels b/c)
DIRECTION_COLS <- c("Increased" = "#F4845F", "Decreased" = "#8b769a")
GENE_TOP_COL   <- "#50B878"
GENE_LABEL_COL <- "#2E7D4F"

# --- Single-cluster trait barplot (top-N signed loadings) ---
# font_scale multiplies every text size (default 1 = unchanged); >1 enlarges fonts
# for figures rendered with many small panels (e.g. the split layout in 08b).
make_trait_panel <- function(dt_k, top_n = 10, font_scale = 1) {
  dt_k <- as.data.table(copy(dt_k))
  setorder(dt_k, -abs_loading)
  top <- dt_k[seq_len(min(.N, top_n))]
  top[, direction := fifelse(loading > 0, "Increased", "Decreased")]
  # Order bars by signed loading (smallest at the bottom after coord_flip)
  top[, trait := factor(trait, levels = top[order(loading), trait])]

  ggplot(top, aes(x = trait, y = loading, fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = DIRECTION_COLS, name = "Cluster Weight") +
    coord_flip() +
    labs(x = NULL, y = "Cluster Weight") +
    theme_minimal(base_size = 11 * font_scale) +
    theme(
      axis.text.y        = element_text(size = 10 * font_scale),
      axis.text.x        = element_text(size = 10 * font_scale),
      axis.title.x       = element_text(size = 11 * font_scale),
      legend.position    = "none",
      panel.grid.major.y = element_blank(),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      # extra left margin so the (longest) trait labels aren't clipped
      plot.margin        = margin(t = 5, r = 5, b = 5, l = 12)
    )
}

# --- Single-cluster gene scatter (max weight vs specificity, top-N labelled) ---
# font_scale multiplies every text size (default 1 = unchanged); see make_trait_panel.
make_gene_panel <- function(scatter_k, top_n = 5, font_scale = 1) {
  ggplot(scatter_k, aes(x = max_loading, y = specificity)) +
    geom_point(data = ~ filter_out(.x, is_top),
               color = "grey70", alpha = 0.3, size = 1.4) +
    geom_point(data = ~ filter(.x, is_top),
               color = GENE_TOP_COL, size = 2.4) +
    ggrepel::geom_text_repel(
      data = ~ filter(.x, is_top),
      aes(label = gene_name), size = 3.6 * font_scale, fontface = "bold",
      color = GENE_LABEL_COL, min.segment.length = 0, max.overlaps = Inf,
      box.padding = 0.5, seed = 42) +
    labs(x = "Max Weight", y = "Gene Specificity") +
    theme_bw(base_size = 11 * font_scale) +
    theme(
      axis.text        = element_text(size = 10 * font_scale),
      axis.title       = element_text(size = 11 * font_scale),
      legend.position  = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

# --- Build the combined figure for one ancestry ---
build_combined <- function(results_dir, anc, gtf_path, cluster_labels,
                           trait_top_n = 10, gene_top_n = 5, ncol_clusters = 2) {
  # Trait loadings (long) and gene scatter (long) — both per cluster
  h_file <- file.path(results_dir, anc, sprintf("H_matrix_%s.tsv", anc))
  trait_dt   <- build_trait_loadings(h_file, strip_source = TRUE)
  scatter_df <- build_gene_scatter_df(results_dir, anc, gtf_path,
                                      cluster_labels = cluster_labels,
                                      top_n = gene_top_n)

  # Cluster display order: desired META order (Glycemic -> Lpa) for META,
  # otherwise natural K order (K1, K2, ... K10).
  kids <- unique(trait_dt$cluster)
  num_order <- kids[order(as.integer(str_replace(kids, "^K", "")))]
  if (anc == "META") {
    desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
    kids <- c(desired_k_order[desired_k_order %in% kids],
              setdiff(num_order, desired_k_order))
  } else {
    kids <- num_order
  }

  # One composite (trait | gene) under a single cluster header per cluster
  composites <- map(kids, function(kid) {
    facet_key <- paste(anc, kid)
    label <- if (!is.null(cluster_labels[[anc]])) {
      cluster_labels[[anc]][[facet_key]] %||% facet_key
    } else facet_key

    trait_p <- make_trait_panel(trait_dt[cluster == kid], top_n = trait_top_n)
    gene_p  <- make_gene_panel(scatter_df %>% filter(cluster_label == label),
                               top_n = gene_top_n)

    # trait panel slightly narrower than the gene panel so its y-axis labels fit
    body <- plot_grid(trait_p, gene_p, ncol = 2, rel_widths = c(0.9, 1.1))
    header <- ggdraw() + draw_label(label, fontface = "bold", size = 15)
    plot_grid(header, body, ncol = 1, rel_heights = c(0.12, 1))
  })

  grid <- plot_grid(plotlist = composites, ncol = ncol_clusters)

  # Shared trait-direction legend (extracted once, placed at the bottom)
  legend_source <- ggplot(
      data.frame(direction = factor(c("Increased", "Decreased"),
                                    levels = c("Increased", "Decreased")),
                 x = c(1, 2), y = c(1, 1)),
      aes(x, y, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = DIRECTION_COLS, name = "Cluster Weight") +
    # title centred above the swatches so it can't be clipped next to the keys
    guides(fill = guide_legend(keywidth = unit(1.1, "cm"),
                               keyheight = unit(1.1, "cm"),
                               title.position = "top", title.hjust = 0.5)) +
    theme_void() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title    = element_text(size = 18, face = "bold", hjust = 0.5,
                                         margin = margin(b = 6)),
          # right margin on each label separates the two swatches/labels
          legend.text     = element_text(size = 16, margin = margin(r = 24, l = 6)),
          legend.spacing.x = unit(0.5, "cm"),
          legend.margin   = margin(t = 6, r = 12, b = 6, l = 12))
  shared_legend <- cowplot::get_plot_component(legend_source, "guide-box-bottom")

  overall_title <- ggdraw() +
    draw_label("Trait and Gene Contributions to Clusters",
               fontface = "bold", size = 22)

  plot_grid(overall_title, grid, shared_legend, ncol = 1,
            rel_heights = c(0.03, 1, 0.07))
}

# --- Generate ---
cat(sprintf("\n=== Generating combined trait+gene figure: ancestry=%s ===\n", ancestry))
lay <- layout_config[[ancestry]] %||% default_layout

p_combined <- build_combined(results_dir, ancestry, gtf_path,
                             cluster_labels = cluster_labels,
                             trait_top_n = 10, gene_top_n = 5,
                             ncol_clusters = lay$ncol)

fname <- sprintf("barplot_%s_combined_trait_gene.png", ancestry)
# bg = "white": cowplot's canvas + the ggdraw header/title strips are transparent;
# without this they render black and hide the (dark) header text.
ggsave(file.path(figures_dir, fname), p_combined,
       width = lay$width, height = lay$height, dpi = 300, bg = "white")
cat(sprintf("  Saved: %s\n", file.path(figures_dir, fname)))

cat("\nDone.\n")
