#!/usr/bin/env Rscript
# 02_visualize_results.R
# Generates all visualizations for A1 bNMF analysis:
#   Panel A: Heatmap grid of cluster-trait weights
#   Panel B: Spearman rank correlation dot plot
#   Panel C: Bar plots of pos/neg weights per cluster
#   Panel D: Cluster frequency across ancestries
#   Panel E: Pearson correlation dot plot
#   Panel F: Variant selection frequency bar plot
#   Panel G: Gene-level contribution bar plots
#   Cross-ancestry: Trait-correlation scatter per matched cluster pair
#   Supplementary: bNMF convergence scatter
#
# Usage:
#   Rscript scripts/a1_analysis/02_visualize_results.R \
#     --config config/a1_config.yaml

library(tidyverse)
library(data.table)
library(patchwork)
library(yaml)

# --- Parse CLI args ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/a1_config.yaml"
if ("--config" %in% args) config_path <- args[which(args == "--config") + 1]

cfg <- read_yaml(config_path)

project_root <- getwd()
script_dir <- "scripts/a1_analysis"
source(file.path(project_root, script_dir, "figure_utils.R"))

results_dir <- file.path(project_root, cfg$results_dir)
figures_dir <- file.path(project_root, cfg$figures_dir)
ancestries  <- cfg$ancestries
gtf_path    <- cfg$gtf_path

# Resolve trait names from config (union across ancestries)
gwas_traits <- unique(unlist(lapply(cfg$trait_gwas, names)))

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# Load data
# ============================================================================

cat("Loading H matrices...\n")
h_list <- load_h_matrices(results_dir, ancestries)
if (length(h_list) == 0) stop("No H matrices found. Check results_dir.")

# Filter out ancestries with K=0
has_active <- function(h) {
  cols <- setdiff(colnames(h), "Cluster")
  length(cols) > 0 && any(as.matrix(h[, cols]) != 0)
}
h_list <- Filter(has_active, h_list)
if (length(h_list) == 0) stop("All ancestries have K=0.")
active_ancestries <- names(h_list)
cat(sprintf("Active ancestries (K>0): %s\n", paste(active_ancestries, collapse = ", ")))

# Cluster matching
cat("Matching clusters across ancestries...\n")
cluster_matches <- match_all_clusters(h_list, sim_threshold = 0.3)


# ============================================================================
# Panel A: Heatmap grid
# ============================================================================

plot_panel_a <- function(h_list) {
  h_all <- map_dfr(names(h_list), function(anc) {
    h_df <- h_list[[anc]]
    n_clusters <- nrow(h_df)
    trait_cols <- setdiff(colnames(h_df), "Cluster")
    # Extract unique trait names (strip _pos/_neg suffix)
    traits <- unique(str_replace(trait_cols, "_(pos|neg)$", ""))

    # Compute net weight per trait: pos - neg (keep raw if only one direction exists)
    net_df <- map_dfc(traits, function(tr) {
      pos_col <- paste0(tr, "_pos")
      neg_col <- paste0(tr, "_neg")
      has_pos <- pos_col %in% trait_cols
      has_neg <- neg_col %in% trait_cols
      vals <- if (has_pos && has_neg) {
        h_df[[pos_col]] - h_df[[neg_col]]
      } else if (has_pos) {
        h_df[[pos_col]]
      } else {
        -h_df[[neg_col]]
      }
      tibble(!!tr := vals)
    })

    bind_cols(tibble(Cluster = h_df$Cluster), net_df) %>%
      pivot_longer(-Cluster, names_to = "trait", values_to = "net_weight") %>%
      mutate(ancestry = anc,
             ancestry_label = sprintf("%s (K=%d)", anc, n_clusters))
  }) %>%
    mutate(Cluster = factor(Cluster, levels = sort(unique(Cluster))))

  ggplot(h_all, aes(x = trait, y = Cluster, fill = net_weight)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3",
                         midpoint = 0, name = "Net Weight\n(pos \u2212 neg)") +
    facet_wrap(~ ancestry_label, ncol = 2, scales = "free_y") +
    theme_big_text() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          strip.text = element_text(face = "bold")) +
    labs(x = "Trait", y = "Cluster")
}


# ============================================================================
# Panel B: Spearman rank correlation
# ============================================================================

plot_panel_b <- function(h_list, cluster_matches) {
  valid <- cluster_matches %>% filter_out(when_any(is.na(cluster_1), is.na(cluster_2)))
  if (nrow(valid) == 0) {
    return(ggplot() + annotate("text", x = .5, y = .5, label = "No matched clusters", size = 6) +
             theme_big_text() + labs(title = "Spearman Rank Correlation"))
  }
  spearman_df <- map_dfr(seq_len(nrow(valid)), function(i) {
    h1 <- h_list[[valid$anc1[i]]]; h2 <- h_list[[valid$anc2[i]]]
    common_cols <- intersect(setdiff(colnames(h1), "Cluster"), setdiff(colnames(h2), "Cluster"))
    vec1 <- h1 %>% filter(Cluster == valid$cluster_1[i]) %>% select(all_of(common_cols)) %>% unlist(use.names = FALSE)
    vec2 <- h2 %>% filter(Cluster == valid$cluster_2[i]) %>% select(all_of(common_cols)) %>% unlist(use.names = FALSE)
    tibble(pair_label = paste(valid$anc1[i], "-", valid$anc2[i]),
           match_label = paste(valid$cluster_1[i], "-", valid$cluster_2[i]),
           spearman_rho = cor(rank(-vec1), rank(-vec2), method = "spearman"))
  })
  ggplot(spearman_df, aes(x = pair_label, y = match_label,
                           size = abs(spearman_rho), color = spearman_rho)) +
    geom_point() +
    scale_color_gradient2(low = "#D94A4A", mid = "white", high = "#4A90D9",
                          midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
    scale_size_continuous(range = c(4, 14), limits = c(0, 1), name = "|rho|") +
    theme_big_text() + theme(panel.grid.major = element_blank()) +
    labs(x = "Ancestry Pair", y = "Matched Clusters")
}


# ============================================================================
# Panel C: Pos/neg bar plots
# ============================================================================

plot_panel_c <- function(h_list) {
  bar_df <- map_dfr(names(h_list), function(anc) {
    h_list[[anc]] %>%
      pivot_longer(-Cluster, names_to = "trait_col", values_to = "weight") %>%
      mutate(ancestry = anc,
             trait = str_replace(trait_col, "_(pos|neg)$", ""),
             direction = str_extract(trait_col, "(pos|neg)$"),
             weight_display = if_else(direction == "neg", -weight, weight)) %>%
      rename(cluster = Cluster) %>%
      select(ancestry, cluster, trait, direction, weight_display)
  }) %>%
    mutate(facet_label = paste(ancestry, "/", cluster))

  ggplot(bar_df, aes(x = trait, y = weight_display, fill = direction)) +
    geom_col(position = "identity", width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c("pos" = "#4A90D9", "neg" = "#D94A4A"),
                      labels = c("pos" = "Positive", "neg" = "Negative"), name = "Direction") +
    facet_wrap(~ facet_label, scales = "free_y") +
    theme_big_text() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(x = "Trait", y = "Weight")
}


# ============================================================================
# Panel D: Cluster frequency across ancestries
# ============================================================================

plot_panel_d <- function(h_list, cluster_matches) {
  all_nodes <- map_dfr(names(h_list), function(anc) {
    tibble(ancestry = anc, cluster = h_list[[anc]]$Cluster)
  }) %>% mutate(node_id = paste(ancestry, cluster, sep = "_"))

  valid_matches <- cluster_matches %>% filter_out(when_any(is.na(cluster_1), is.na(cluster_2)))

  parent <- set_names(all_nodes$node_id, all_nodes$node_id)
  find_root <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
  union_nodes <- function(a, b) { ra <- find_root(a); rb <- find_root(b); if (ra != rb) parent[ra] <<- rb }

  for (i in seq_len(nrow(valid_matches))) {
    a <- paste(valid_matches$anc1[i], valid_matches$cluster_1[i], sep = "_")
    b <- paste(valid_matches$anc2[i], valid_matches$cluster_2[i], sep = "_")
    if (a %in% names(parent) && b %in% names(parent)) union_nodes(a, b)
  }
  all_nodes <- all_nodes %>% mutate(canonical = map_chr(node_id, find_root))
  lbl_map <- set_names(paste("Cluster", seq_along(unique(all_nodes$canonical))), unique(all_nodes$canonical))
  all_nodes <- all_nodes %>% mutate(canonical_label = lbl_map[canonical])
  freq_df <- all_nodes %>% group_by(canonical_label) %>% summarise(count = n_distinct(ancestry), .groups = "drop")

  ggplot(freq_df, aes(x = count, y = reorder(canonical_label, count))) +
    geom_col(fill = "#5B8FA8", width = 0.6) +
    scale_x_continuous(breaks = seq(0, length(h_list)), limits = c(0, length(h_list) + 0.5)) +
    coord_flip() + theme_big_text() +
    labs(x = "Number of Ancestries", y = "Cluster")
}


# ============================================================================
# Panel E: Pearson correlation dot plot
# ============================================================================

plot_panel_e <- function(h_list, cluster_matches) {
  valid <- cluster_matches %>% filter_out(when_any(is.na(cluster_1), is.na(cluster_2)))
  if (nrow(valid) == 0) {
    return(ggplot() + annotate("text", x = .5, y = .5, label = "No matched clusters", size = 6) +
             theme_big_text() + labs(title = "Pearson Correlation"))
  }
  pearson_df <- map_dfr(seq_len(nrow(valid)), function(i) {
    h1 <- h_list[[valid$anc1[i]]]; h2 <- h_list[[valid$anc2[i]]]
    common_cols <- intersect(setdiff(colnames(h1), "Cluster"), setdiff(colnames(h2), "Cluster"))
    vec1 <- h1 %>% filter(Cluster == valid$cluster_1[i]) %>% select(all_of(common_cols)) %>% unlist(use.names = FALSE)
    vec2 <- h2 %>% filter(Cluster == valid$cluster_2[i]) %>% select(all_of(common_cols)) %>% unlist(use.names = FALSE)
    tibble(pair_label = paste(valid$anc1[i], "-", valid$anc2[i]),
           match_label = paste(valid$cluster_1[i], "-", valid$cluster_2[i]),
           pearson_r = cor(vec1, vec2, method = "pearson"))
  })
  ggplot(pearson_df, aes(x = pair_label, y = match_label, size = abs(pearson_r), color = pearson_r)) +
    geom_point() +
    scale_color_gradient2(low = "#D94A4A", mid = "white", high = "#4A90D9",
                          midpoint = 0, limits = c(-1, 1), name = "Pearson\nr") +
    scale_size_continuous(range = c(4, 14), limits = c(0, 1), name = "|r|") +
    theme_big_text() + theme(panel.grid.major = element_blank()) +
    labs(x = "Ancestry Pair", y = "Matched Clusters")
}


# ============================================================================
# Panel F: Variant selection frequency
# ============================================================================

plot_panel_f <- function(results_dir, ancestries) {
  top_n <- 10
  top_entries <- map_dfr(ancestries, function(anc) {
    w_path <- file.path(results_dir, anc, sprintf("W_matrix_%s.tsv", anc))
    if (!file.exists(w_path)) return(NULL)
    w_df <- read_tsv(w_path, show_col_types = FALSE)
    cluster_cols <- setdiff(colnames(w_df), "VAR_ID")
    map_dfr(cluster_cols, function(kcol) {
      w_df %>% arrange(desc(.data[[kcol]])) %>%
        slice_head(n = min(top_n, nrow(w_df))) %>%
        select(VAR_ID) %>% mutate(ancestry = anc, cluster = kcol)
    })
  })
  if (is.null(top_entries) || nrow(top_entries) == 0) {
    return(ggplot() + annotate("text", x = .5, y = .5, label = "No W matrices found", size = 6) +
             theme_big_text() + labs(title = "Variant Selection Frequency"))
  }
  var_counts <- top_entries %>% count(VAR_ID, name = "times_selected") %>%
    arrange(desc(times_selected), VAR_ID) %>% slice_head(n = top_n) %>%
    mutate(VAR_ID = fct_inorder(VAR_ID))

  ggplot(var_counts, aes(x = VAR_ID, y = times_selected)) +
    geom_col(fill = "gray50", width = 0.7) +
    scale_y_continuous(breaks = seq(0, max(var_counts$times_selected))) +
    theme_big_text() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(x = "Variant", y = "Times Selected (Top 10)")
}


# ============================================================================
# Panel G: Gene-level contribution bar plots
# ============================================================================

plot_panel_g <- function(results_dir, ancestries, gtf_path,
                         metric = c("mean", "max", "specificity"),
                         cluster_labels = NULL, top_n = 10) {
  metric <- match.arg(metric)
  metric_label <- switch(metric,
    mean        = "Normalized Weight",
    max         = "Max Loading",
    specificity = "Cluster Specificity"
  )

  cat("  Parsing GTF for protein-coding genes...\n")
  gene_df <- parse_gtf_genes(gtf_path)
  cat(sprintf("  Found %d protein-coding genes\n", nrow(gene_df)))

  track_data_list <- list()
  for (anc in ancestries) {
    w_path <- file.path(results_dir, anc, sprintf("W_matrix_%s.tsv", anc))
    if (!file.exists(w_path)) next
    w_df <- read_tsv(w_path, show_col_types = FALSE) %>%
      mutate(chr = str_replace(VAR_ID, "^([^_]+)_.*", "\\1"),
             pos = as.numeric(str_replace(VAR_ID, "^[^_]+_([^_]+)_.*", "\\1")))
    cluster_cols <- setdiff(colnames(w_df), c("VAR_ID", "chr", "pos"))

    # Pre-compute row sums across all clusters (used for specificity metric)
    if (metric == "specificity") {
      w_df <- w_df %>%
        mutate(.row_total = rowSums(across(all_of(cluster_cols)), na.rm = TRUE))
    }

    mapped <- map_snps_to_genes(w_df, gene_df)
    cat(sprintf("  %s: %d SNPs mapped to genes (of %d total)\n",
                anc, n_distinct(mapped$VAR_ID), nrow(w_df)))
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
  }

  if (length(track_data_list) == 0) {
    return(list())
  }

  # Group bar plots by ancestry and return one combined plot per ancestry
  anc_plots <- list()
  for (anc in ancestries) {
    anc_keys <- grep(paste0("^", anc, " "), names(track_data_list), value = TRUE)
    anc_entries <- track_data_list[anc_keys]
    if (length(anc_entries) == 0) next
    bar_plots <- map(anc_entries, function(ti) {
      td <- ti$data %>% arrange(norm_weight) %>% mutate(gene_name = fct_inorder(gene_name))
      ggplot(td, aes(x = gene_name, y = norm_weight)) +
        geom_col(fill = "#50B878", width = 0.7) + coord_flip() +
        theme_big_text() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(ti$label) + labs(x = NULL, y = metric_label)
    })
    anc_plots[[anc]] <- patchwork::wrap_plots(bar_plots, ncol = 2) +
      patchwork::plot_annotation(title = sprintf("Gene-Level Contributions: %s", anc),
                                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
  }
  anc_plots
}


# ============================================================================
# Cross-ancestry: Trait-correlation scatter per matched pair
# ============================================================================

plot_cross_ancestry_scatter <- function(h_list, gwas_traits) {
  use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
  anc_names <- names(h_list)
  if (length(anc_names) < 2) return(NULL)

  pairs <- combn(anc_names, 2, simplify = FALSE)
  scatter_plots <- list()

  for (pair in pairs) {
    anc1 <- pair[1]; anc2 <- pair[2]
    h1 <- h_list[[anc1]]; h2 <- h_list[[anc2]]

    match_result <- match_clusters_by_trait_correlation(h1, h2)
    assignment <- match_result$assignment
    common_traits <- match_result$common_traits

    for (row_i in seq_len(nrow(assignment))) {
      cl1 <- assignment$cluster_1[row_i]
      cl2 <- assignment$cluster_2[row_i]
      r_val <- assignment$pearson_r[row_i]

      net1 <- match_result$net_vecs_1[[cl1]]
      net2 <- match_result$net_vecs_2[[cl2]]

      scatter_df <- tibble(
        trait = common_traits,
        weight_1 = net1[common_traits],
        weight_2 = net2[common_traits],
        point_type = factor(ifelse(common_traits %in% gwas_traits, "GWAS trait", "Reference"),
                            levels = c("GWAS trait", "Reference"))
      )

      p <- ggplot(scatter_df, aes(x = weight_1, y = weight_2, color = point_type)) +
        geom_point(size = 3) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        scale_color_manual(values = c("GWAS trait" = "#E74C3C", "Reference" = "#3498DB"),
                           name = "Type") +
        theme_big_text() +
        labs(x = sprintf("%s %s net weight", anc1, cl1),
             y = sprintf("%s %s net weight", anc2, cl2),
             title = sprintf("%s %s vs %s %s (r = %.2f)", anc1, cl1, anc2, cl2, r_val))

      if (use_ggrepel) {
        p <- p + ggrepel::geom_text_repel(aes(label = trait), size = 3,
                                           max.overlaps = Inf, show.legend = FALSE)
      }
      scatter_plots[[paste(anc1, cl1, anc2, cl2)]] <- p
    }
  }
  scatter_plots
}


# ============================================================================
# Supplementary: bNMF convergence scatter
# ============================================================================

plot_convergence <- function(results_dir, ancestries) {
  dt <- rbindlist(lapply(ancestries, function(anc) {
    path <- file.path(results_dir, anc, sprintf("run_summary_%s.csv", anc))
    if (!file.exists(path)) return(NULL)
    d <- fread(path)
    d[, ancestry := anc]
  }), fill = TRUE)

  if (nrow(dt) == 0) return(NULL)
  dt[, ancestry := factor(ancestry, levels = ancestries)]
  dt[, is_optimal_k := as.logical(is_optimal_k)]

  ggplot(dt, aes(x = iterations, y = K_converged, color = ancestry,
                  shape = is_optimal_k)) +
    geom_jitter(size = 3, stroke = 1, width = 0, height = 0.1) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "Yes", "FALSE" = "No"), name = "Optimal K") +
    scale_y_continuous(breaks = seq(0, max(dt$K_converged))) +
    labs(x = "Iterations", y = expression(K[converged]), color = "Ancestry") +
    theme_big_text() + theme(legend.position = "right")
}


# ============================================================================
# Generate all figures
# ============================================================================

cat("Generating panels...\n")

cat("  Panel A: Heatmap...\n")
pa <- plot_panel_a(h_list)

cat("  Panel B: Bar plots...\n")
pc <- plot_panel_c(h_list)

cat("  Panel C: Gene bars...\n")
pg_list <- plot_panel_g(results_dir, active_ancestries, gtf_path)


# --- Save individual panels ---
cat("Saving individual panel PNGs...\n")
ggsave(file.path(figures_dir, "panel_A_heatmap.png"), pa,
       width = 14, height = 8, dpi = 300)
cat("  Saved: panel_A_heatmap.png\n")
ggsave(file.path(figures_dir, "panel_B_posneg_bars.png"), pc,
       width = 14, height = 10, dpi = 300)
cat("  Saved: panel_B_posneg_bars.png\n")

# Save one gene bars plot per ancestry
for (anc_name in names(pg_list)) {
  fname <- sprintf("panel_C_gene_bars_%s.png", anc_name)
  ggsave(file.path(figures_dir, fname), pg_list[[anc_name]],
         width = 14, height = 10, dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}


# --- Main figure assembly (use first ancestry gene bars as representative) ---
cat("Assembling main figure...\n")
pg_representative <- if (length(pg_list) > 0) pg_list[[1]] else ggplot() + theme_void()
main_fig <- wrap_plots(list(pa, pc, pg_representative), ncol = 1) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 22, face = "bold")))

ggsave(file.path(figures_dir, "figure_main.png"), main_fig,
       width = 16, height = 28, dpi = 300)
cat("Saved: figure_main.png\n")


# --- Cross-ancestry scatter plots ---
cat("  Cross-ancestry scatter...\n")
scatter_plots <- plot_cross_ancestry_scatter(h_list, gwas_traits)
if (!is.null(scatter_plots) && length(scatter_plots) > 0) {
  cross_fig <- wrap_plots(scatter_plots, ncol = 2)
  scatter_height <- min(7 * ceiling(length(scatter_plots) / 2), 48)
  ggsave(file.path(figures_dir, "cross_ancestry_scatter.png"), cross_fig,
         width = 14, height = scatter_height, dpi = 300, limitsize = FALSE)
  cat("Saved: cross_ancestry_scatter.png\n")

  # Save individual scatter panels
  for (panel_name in names(scatter_plots)) {
    safe_name <- gsub(" ", "_", panel_name)
    ggsave(file.path(figures_dir, sprintf("panel_scatter_%s.png", safe_name)),
           scatter_plots[[panel_name]], width = 8, height = 7, dpi = 300)
  }
  cat(sprintf("  Saved %d individual scatter panels\n", length(scatter_plots)))
}


# --- Supplementary: convergence ---
cat("  Supp: Convergence...\n")
conv_plot <- plot_convergence(results_dir, active_ancestries)
if (!is.null(conv_plot)) {
  ggsave(file.path(figures_dir, "supp_bnmf_convergence.png"), conv_plot,
         width = 8, height = 6, dpi = 300)
  cat("Saved: supp_bnmf_convergence.png\n")
}

cat("\n=== Visualization complete! ===\n")
