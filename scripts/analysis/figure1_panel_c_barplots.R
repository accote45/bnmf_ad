# figure1_panel_c_barplots.R
# Panel C: Grid of bar plots showing positive (blue, up) and negative (red, down)
# cluster-trait weights, one subplot per ancestry-cluster combination.

# Assumes figure_utils.R has been sourced.


#' Create Panel C: Bar plots of cluster weights per ancestry-cluster
#'
#' For each ancestry-cluster combination, shows separate bars for _pos (blue, up)
#' and _neg (red, pointing down from zero) weights for each trait.
#'
#' @param h_list Named list of H-matrix tibbles (from load_h_matrices)
#' @return ggplot object
plot_panel_c <- function(h_list) {

  bar_df <- map_dfr(names(h_list), function(anc) {
    h_list[[anc]] %>%
      pivot_longer(-Cluster, names_to = "trait_col", values_to = "weight") %>%
      mutate(
        ancestry  = anc,
        trait     = str_replace(trait_col, "_(pos|neg)$", ""),
        direction = str_extract(trait_col, "(pos|neg)$"),
        weight_display = if_else(direction == "neg", -weight, weight)
      ) %>%
      rename(cluster = Cluster) %>%
      select(ancestry, cluster, trait, direction, weight_display)
  }) %>%
    mutate(facet_label = paste(ancestry, "/", cluster))

  ggplot(bar_df, aes(x = trait, y = weight_display, fill = direction)) +
    geom_col(position = "identity", width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(
      values = c("pos" = "#4A90D9", "neg" = "#D94A4A"),
      labels = c("pos" = "Positive", "neg" = "Negative"),
      name   = "Direction"
    ) +
    facet_wrap(~ facet_label, scales = "free_y") +
    theme_big_text() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = "Trait", y = "Weight")
}
