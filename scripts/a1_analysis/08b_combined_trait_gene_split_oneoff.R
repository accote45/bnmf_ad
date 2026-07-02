#!/usr/bin/env Rscript
# 08b_combined_trait_gene_split_oneoff.R
# ONE-OFF / EXPLORATORY: alternative layout of barplot_<anc>_combined_trait_gene
# for visual comparison against the current figure (08_combined_trait_gene_facets.R).
#
# Difference vs 08: instead of pairing each cluster's trait barplot and gene scatter
# side-by-side, this pools ALL barplots into the LEFT half and ALL scatters into the
# RIGHT half. Each cluster header is duplicated so it appears above its barplot (left)
# AND above its scatter (right).
#
# It REUSES 08's panel builders / data-prep / labels / colors by sourcing 08 (which
# also regenerates the current figure, handy for comparison). No logic is reimplemented;
# only the top-level cowplot assembly differs. Writes a separate *_split.png so the
# current figure is left untouched.
#
# Usage:
#   Rscript scripts/a1_analysis/08b_combined_trait_gene_split_oneoff.R \
#     --config config/a1_config.yaml --ancestry META

# --- Source 08: defines make_trait_panel/make_gene_panel, build_trait_loadings/
#     build_gene_scatter_df (via figure_utils.R), cluster_labels, DIRECTION_COLS,
#     layout_config/default_layout, results_dir, gtf_path, figures_dir.
#     08 reads commandArgs(), so the same --config/--ancestry flags flow through. ---
source(file.path(getwd(), "scripts", "a1_analysis", "08_combined_trait_gene_facets.R"))

# --- Re-parse CLI (08 already did, but keep this script self-describing) ---
args <- commandArgs(trailingOnly = TRUE)
anc  <- "META"
if ("--ancestry" %in% args) anc <- args[which(args == "--ancestry") + 1]

# --- Split-layout assembler: reuses 08's make_trait_panel / make_gene_panel ---
build_combined_split <- function(results_dir, anc, gtf_path, cluster_labels,
                                 trait_top_n = 10, gene_top_n = 5, ncol_clusters = 2,
                                 font_scale = 1.5) {
  # Same data prep as 08's build_combined
  h_file <- file.path(results_dir, anc, sprintf("H_matrix_%s.tsv", anc))
  trait_dt   <- build_trait_loadings(h_file, strip_source = TRUE)
  scatter_df <- build_gene_scatter_df(results_dir, anc, gtf_path,
                                      cluster_labels = cluster_labels,
                                      top_n = gene_top_n)

  # Same cluster display order as 08 (desired META order; else natural K order)
  kids <- unique(trait_dt$cluster)
  num_order <- kids[order(as.integer(str_replace(kids, "^K", "")))]
  if (anc == "META") {
    desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
    kids <- c(desired_k_order[desired_k_order %in% kids],
              setdiff(num_order, desired_k_order))
  } else {
    kids <- num_order
  }

  # Resolve the display label for a cluster id (same lookup 08 uses)
  label_for <- function(kid) {
    facet_key <- paste(anc, kid)
    if (!is.null(cluster_labels[[anc]])) {
      cluster_labels[[anc]][[facet_key]] %||% facet_key
    } else facet_key
  }
  # header strip identical to 08 (bold, size 15) — built per side so labels duplicate
  header_for <- function(label) ggdraw() + draw_label(label, fontface = "bold", size = 15)

  # LEFT: header + trait barplot, stacked per cluster, tiled 2 cluster-cols x 5 rows
  left_cells <- map(kids, function(kid) {
    p <- make_trait_panel(trait_dt[cluster == kid], top_n = trait_top_n,
                          font_scale = font_scale)
    plot_grid(header_for(label_for(kid)), p, ncol = 1, rel_heights = c(0.12, 1))
  })
  # RIGHT: header + gene scatter, same tiling; header duplicated from the left side
  right_cells <- map(kids, function(kid) {
    p <- make_gene_panel(scatter_df %>% filter(cluster_label == label_for(kid)),
                         top_n = gene_top_n, font_scale = font_scale)
    plot_grid(header_for(label_for(kid)), p, ncol = 1, rel_heights = c(0.12, 1))
  })

  left_grid  <- plot_grid(plotlist = left_cells,  ncol = ncol_clusters)
  right_grid <- plot_grid(plotlist = right_cells, ncol = ncol_clusters)

  # Combine halves; 0.9/1.1 mirrors 08's per-composite trait/gene width ratio
  body <- plot_grid(left_grid, right_grid, ncol = 2, rel_widths = c(0.9, 1.1))

  # --- Shared legend + overall title (presentation glue, mirrors 08:167-194) ---
  legend_source <- ggplot(
      data.frame(direction = factor(c("Increased", "Decreased"),
                                    levels = c("Increased", "Decreased")),
                 x = c(1, 2), y = c(1, 1)),
      aes(x, y, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = DIRECTION_COLS, name = "Cluster Weight") +
    guides(fill = guide_legend(keywidth = unit(1.1, "cm"),
                               keyheight = unit(1.1, "cm"),
                               title.position = "top", title.hjust = 0.5)) +
    theme_void() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title    = element_text(size = 18, face = "bold", hjust = 0.5,
                                         margin = margin(b = 6)),
          legend.text     = element_text(size = 16, margin = margin(r = 24, l = 6)),
          legend.spacing.x = unit(0.5, "cm"),
          legend.margin   = margin(t = 6, r = 12, b = 6, l = 12))
  shared_legend <- cowplot::get_plot_component(legend_source, "guide-box-bottom")

  overall_title <- ggdraw() +
    draw_label("Trait and Gene Contributions to Clusters",
               fontface = "bold", size = 22)

  plot_grid(overall_title, body, shared_legend, ncol = 1,
            rel_heights = c(0.03, 1, 0.07))
}

# --- Generate split variant ---
cat(sprintf("\n=== Generating SPLIT trait+gene figure (one-off): ancestry=%s ===\n", anc))
lay <- layout_config[[anc]] %||% default_layout

p_split <- build_combined_split(results_dir, anc, gtf_path,
                                cluster_labels = cluster_labels,
                                trait_top_n = 10, gene_top_n = 5,
                                ncol_clusters = lay$ncol)

# narrower than the old 44 in to remove the horizontal stretch; height matches the
# original per-cluster top-10 layout. Tunable.
split_width  <- 26
split_height <- 24
fname <- sprintf("barplot_%s_combined_trait_gene_split.png", anc)
ggsave(file.path(figures_dir, fname), p_split,
       width = split_width, height = split_height, dpi = 300, bg = "white")
cat(sprintf("  Saved: %s\n", file.path(figures_dir, fname)))

cat("\nDone.\n")
