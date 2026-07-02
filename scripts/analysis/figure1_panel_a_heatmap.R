# figure1_panel_a_heatmap.R
# Panel A: Grid of cluster-trait heatmaps, one per ancestry.
# Recreates H-matrix heatmaps using ggplot2 (instead of pheatmap) for
# consistent styling across the multi-panel figure.

# Assumes figure_utils.R has been sourced.


#' Create Panel A: Heatmap grid of cluster-trait weights
#'
#' @param h_list Named list of H-matrix tibbles (from load_h_matrices)
#' @return ggplot object
plot_panel_a <- function(h_list) {

  # Melt each ancestry's H-matrix to long format and bind
  h_all <- map_dfr(names(h_list), function(anc) {
    h_df <- h_list[[anc]]
    n_clusters <- nrow(h_df)

    # Per-ancestry normalization: scale each ancestry's weights to [0,1]
    # so that ancestries with very small absolute values are still visible
    mat <- h_df %>% select(-Cluster) %>% as.matrix()
    row_max <- max(mat)
    if (row_max > 0) mat <- mat / row_max

    bind_cols(tibble(Cluster = h_df$Cluster), as_tibble(mat)) %>%
      pivot_longer(-Cluster, names_to = "trait_col", values_to = "weight") %>%
      mutate(
        ancestry = anc,
        ancestry_label = sprintf("%s (K=%d)", anc, n_clusters)
      )
  }) %>%
    mutate(Cluster = factor(Cluster, levels = sort(unique(Cluster))))

  ggplot(h_all, aes(x = trait_col, y = Cluster, fill = weight)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient(low = "white", high = "darkblue",
                        name = "Normalized\nWeight") +
    facet_wrap(~ ancestry_label, ncol = 2, scales = "free_y") +
    theme_big_text() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      strip.text  = element_text(face = "bold")
    ) +
    labs(x = "Trait", y = "Cluster")
}
