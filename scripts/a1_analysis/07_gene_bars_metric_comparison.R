#!/usr/bin/env Rscript
# 07_gene_bars_metric_comparison.R
# Generate gene contribution barplots using max-loading and cluster-specificity
# metrics for comparison. Uses the same aesthetics as the existing panel_C plots.
#
# Usage:
#   Rscript scripts/a1_analysis/07_gene_bars_metric_comparison.R \
#     --config config/a1_config.yaml --ancestry META

library(tidyverse)
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

# --- Source helpers (theme, gene mapping) ---
source(file.path(project_root, "scripts", "a1_analysis", "figure_utils.R"))

# --- Cluster labels ---
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

# --- Core plotting function ---
# Supports metric = "max" (peak W-loading per gene) or "specificity"
# (mean cluster-specificity score: W[i,k] / sum(W[i,]) for SNPs with W[i,k]>0)
plot_gene_bars <- function(results_dir, anc, gtf_path, metric, cluster_labels,
                           top_n = 5, ncol_grid = 3) {
  metric_label <- switch(metric,
    mean        = "Mean Loading",
    max         = "Max Loading",
    specificity = "Cluster Specificity"
  )
  metric_title <- switch(metric,
    mean        = "Mean Gene Contributions per Cluster",
    max         = "Maximum Gene Contributions per Cluster",
    specificity = "Specificity-scored Gene Contributions per Cluster"
  )

  gene_df <- parse_gtf_genes(gtf_path)
  w_path <- file.path(results_dir, anc, sprintf("W_matrix_%s.tsv", anc))
  if (!file.exists(w_path)) stop(sprintf("W matrix not found: %s", w_path))

  w_df <- read_tsv(w_path, show_col_types = FALSE) %>%
    mutate(chr = str_replace(VAR_ID, "^([^_]+)_.*", "\\1"),
           pos = as.numeric(str_replace(VAR_ID, "^[^_]+_([^_]+)_.*", "\\1")))
  cluster_cols <- setdiff(colnames(w_df), c("VAR_ID", "chr", "pos"))

  # Pre-compute row sums for specificity
  if (metric == "specificity") {
    w_df <- w_df %>%
      mutate(.row_total = rowSums(across(all_of(cluster_cols)), na.rm = TRUE))
  }

  mapped <- map_snps_to_genes(w_df, gene_df)
  cat(sprintf("  %s: %d SNPs mapped to genes (of %d total)\n",
              anc, n_distinct(mapped$VAR_ID), nrow(w_df)))

  track_data_list <- list()
  for (kcol in cluster_cols) {
    gene_weights <- switch(metric,
      mean = mapped %>%
        group_by(gene_name) %>%
        summarise(score = mean(.data[[kcol]], na.rm = TRUE), .groups = "drop"),
      max = mapped %>%
        group_by(gene_name) %>%
        summarise(score = max(.data[[kcol]], na.rm = TRUE), .groups = "drop"),
      specificity = mapped %>%
        filter(.data[[kcol]] > 0, .row_total > 0) %>%
        mutate(.spec = .data[[kcol]] / .row_total) %>%
        group_by(gene_name) %>%
        summarise(score = mean(.spec, na.rm = TRUE), .groups = "drop")
    )
    gene_weights <- gene_weights %>%
      filter(score > 0) %>%
      slice_max(score, n = top_n, with_ties = FALSE) %>%
      mutate(norm_weight = score / max(score))

    facet_key <- paste(anc, kcol)
    display_label <- if (!is.null(cluster_labels[[anc]])) {
      cluster_labels[[anc]][[facet_key]] %||% facet_key
    } else {
      facet_key
    }
    track_data_list[[facet_key]] <- list(data = gene_weights, label = display_label)
  }

  # Build individual bar plots
  bar_plots <- map(track_data_list, function(ti) {
    td <- ti$data %>% arrange(norm_weight) %>% mutate(gene_name = fct_inorder(gene_name))
    ggplot(td, aes(x = gene_name, y = norm_weight)) +
      geom_col(fill = "#50B878", width = 0.7) + coord_flip() +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 11),
            axis.title  = element_text(size = 10),
            panel.background = element_rect(fill = "white"),
            plot.background  = element_rect(fill = "white", color = NA)) +
      ggtitle(ti$label) + labs(x = NULL, y = metric_label)
  })

  n_clusters <- length(bar_plots)
  ncol_grid <- min(ncol_grid, n_clusters)
  title_gg <- ggdraw() +
    draw_label(metric_title, size = 16, fontface = "bold")
  grid <- plot_grid(plotlist = bar_plots, ncol = ncol_grid)
  plot_grid(title_gg, grid, ncol = 1, rel_heights = c(0.04, 1))
}

# --- Lollipop plot: max loading (x-axis) + specificity (dot size) ---
plot_gene_lollipop <- function(results_dir, anc, gtf_path, cluster_labels,
                               top_n = 5, ncol_grid = 3) {
  gene_df <- parse_gtf_genes(gtf_path)
  w_path  <- file.path(results_dir, anc, sprintf("W_matrix_%s.tsv", anc))
  if (!file.exists(w_path)) stop(sprintf("W matrix not found: %s", w_path))

  w_df <- read_tsv(w_path, show_col_types = FALSE) %>%
    mutate(chr = str_replace(VAR_ID, "^([^_]+)_.*", "\\1"),
           pos = as.numeric(str_replace(VAR_ID, "^[^_]+_([^_]+)_.*", "\\1")))
  cluster_cols <- setdiff(colnames(w_df), c("VAR_ID", "chr", "pos"))

  w_df <- w_df %>%
    mutate(.row_total = rowSums(across(all_of(cluster_cols)), na.rm = TRUE))

  mapped <- map_snps_to_genes(w_df, gene_df)
  cat(sprintf("  %s: %d SNPs mapped to genes (of %d total)\n",
              anc, n_distinct(mapped$VAR_ID), nrow(w_df)))

  track_data_list <- list()
  for (kcol in cluster_cols) {
    gene_max <- mapped %>%
      group_by(gene_name) %>%
      summarise(max_loading = max(.data[[kcol]], na.rm = TRUE), .groups = "drop") %>%
      filter(max_loading > 0)

    gene_spec <- mapped %>%
      filter(.data[[kcol]] > 0, .row_total > 0) %>%
      mutate(.spec = .data[[kcol]] / .row_total) %>%
      group_by(gene_name) %>%
      summarise(specificity = mean(.spec, na.rm = TRUE), .groups = "drop")

    gene_combined <- gene_max %>%
      slice_max(max_loading, n = top_n, with_ties = FALSE) %>%
      left_join(gene_spec, by = "gene_name") %>%
      mutate(specificity = replace_na(specificity, 0))

    facet_key <- paste(anc, kcol)
    display_label <- if (!is.null(cluster_labels[[anc]])) {
      cluster_labels[[anc]][[facet_key]] %||% facet_key
    } else {
      facet_key
    }
    track_data_list[[facet_key]] <- list(data = gene_combined, label = display_label)
  }

  global_max_spec <- max(map_dbl(track_data_list, ~ max(.x$data$specificity, na.rm = TRUE)), na.rm = TRUE)

  lollipop_plots <- map(track_data_list, function(ti) {
    td <- ti$data %>%
      arrange(max_loading) %>%
      mutate(gene_name = fct_inorder(gene_name))

    ggplot(td, aes(x = gene_name, y = max_loading)) +
      geom_segment(aes(xend = gene_name, y = 0, yend = max_loading),
                   color = "grey60", linewidth = 0.8) +
      geom_point(aes(size = specificity), color = "#50B878") +
      coord_flip() +
      scale_size_continuous(name = "Specificity", range = c(2, 7),
                            limits = c(0, global_max_spec)) +
      theme_bw(base_size = 11) +
      theme(
        plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y      = element_text(size = 14),
        axis.text.x      = element_text(size = 13),
        axis.title       = element_text(size = 12),
        legend.position  = "none",
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white", color = NA)
      ) +
      ggtitle(ti$label) +
      labs(x = NULL, y = "Max Loading")
  })

  legend_source <- ggplot(data.frame(x = 1, y = 1, specificity = 0.5),
                          aes(x, y, size = specificity)) +
    geom_point(color = "#50B878") +
    scale_size_continuous(name = "Specificity", range = c(2, 7),
                          limits = c(0, global_max_spec)) +
    theme_void() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 13))
  shared_legend <- cowplot::get_plot_component(legend_source, "guide-box-bottom")

  n_clusters <- length(lollipop_plots)
  ncol_grid  <- min(ncol_grid, n_clusters)
  title_gg <- ggdraw() +
    draw_label("Top Gene Contributions per Cluster", size = 20, fontface = "bold")
  grid <- plot_grid(plotlist = lollipop_plots, ncol = ncol_grid)
  plot_grid(title_gg, grid, shared_legend, ncol = 1, rel_heights = c(0.04, 1, 0.06))
}

# --- Scatter plot: max loading (x) vs specificity (y), faceted by cluster ---
# Every mapped gene is plotted (grey, translucent); the top-N genes by max
# loading per cluster are highlighted green and labelled with ggrepel.
plot_gene_scatter <- function(results_dir, anc, gtf_path, cluster_labels,
                              top_n = 5, ncol_grid = 3) {
  # Per-cluster (gene, max_loading, specificity, is_top) table — shared helper
  scatter_df <- build_gene_scatter_df(results_dir, anc, gtf_path,
                                      cluster_labels = cluster_labels, top_n = top_n)

  ncol_grid <- min(ncol_grid, nlevels(scatter_df$cluster_label))

  ggplot(scatter_df, aes(x = max_loading, y = specificity)) +
    geom_point(data = ~ filter_out(.x, is_top),
               color = "grey70", alpha = 0.3, size = 1.6) +
    geom_point(data = ~ filter(.x, is_top),
               color = "#50B878", size = 2.8) +
    ggrepel::geom_text_repel(
      data = ~ filter(.x, is_top),
      aes(label = gene_name), size = 4.5, fontface = "bold",
      color = "#2E7D4F", min.segment.length = 0, max.overlaps = Inf,
      box.padding = 0.5, seed = 42) +
    facet_wrap(~cluster_label, ncol = ncol_grid, scales = "free") +
    labs(title = "Contributions of Genes to Clusters",
         x = "Max Weight", y = "Gene Specificity to Cluster Score") +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 24),
      strip.text       = element_text(face = "bold", size = 18),
      axis.text.y      = element_text(size = 15),
      axis.text.x      = element_text(size = 15),
      axis.title       = element_text(size = 15),
      legend.position  = "none",
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

# --- Layout config per ancestry ---
layout_config <- list(
  META = list(ncol = 2, width = 14, height = 22)
)
default_layout <- list(ncol = 3, width = 18, height = 18)

# --- Generate plots ---
metrics <- c("mean", "max", "specificity")
suffixes <- c(mean = "mean", max = "max", specificity = "specificity")

for (m in metrics) {
  cat(sprintf("\n=== Generating gene bars: metric=%s, ancestry=%s ===\n", m, ancestry))
  lay <- layout_config[[ancestry]] %||% default_layout
  p <- plot_gene_bars(results_dir, ancestry, gtf_path, metric = m,
                      cluster_labels = cluster_labels, ncol_grid = lay$ncol)
  fname <- sprintf("barplot_%s_gene_%s.png", ancestry, suffixes[m])
  ggsave(file.path(figures_dir, fname), p, width = lay$width, height = lay$height, dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

# --- Generate lollipop plot ---
cat(sprintf("\n=== Generating gene lollipop: ancestry=%s ===\n", ancestry))
lay <- layout_config[[ancestry]] %||% default_layout
p_lollipop <- plot_gene_lollipop(results_dir, ancestry, gtf_path,
                                  cluster_labels = cluster_labels,
                                  ncol_grid = lay$ncol)
fname_lollipop <- sprintf("barplot_%s_gene_lollipop.png", ancestry)
ggsave(file.path(figures_dir, fname_lollipop), p_lollipop,
       width = lay$width, height = lay$height, dpi = 300)
cat(sprintf("  Saved: %s\n", fname_lollipop))

# --- Generate scatter plot (max loading vs specificity) ---
cat(sprintf("\n=== Generating gene scatter: ancestry=%s ===\n", ancestry))
lay <- layout_config[[ancestry]] %||% default_layout
p_scatter <- plot_gene_scatter(results_dir, ancestry, gtf_path,
                               cluster_labels = cluster_labels,
                               ncol_grid = lay$ncol)
fname_scatter <- sprintf("barplot_%s_gene_scatter.png", ancestry)
ggsave(file.path(figures_dir, fname_scatter), p_scatter,
       width = lay$width, height = lay$height, dpi = 300)
cat(sprintf("  Saved: %s\n", fname_scatter))

cat("\nDone.\n")
