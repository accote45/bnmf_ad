# figure1_panel_e_pearson.R
# Panel E: Dot plot of Pearson correlation coefficients of cluster weight
# vectors between matched clusters across ancestry pairs.

# Assumes figure_utils.R has been sourced.


#' Create Panel E: Pearson correlation dot plot
#'
#' For each matched cluster pair across ancestries, computes Pearson r
#' of the full H-matrix weight vectors.
#'
#' @param h_list Named list of H-matrix tibbles (from load_h_matrices)
#' @param cluster_matches tibble from match_all_clusters()
#' @return ggplot object
plot_panel_e <- function(h_list, cluster_matches) {

  valid <- cluster_matches %>% filter_out(when_any(is.na(cluster_1), is.na(cluster_2)))

  if (nrow(valid) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No matched clusters",
                 size = 6) +
        theme_big_text() +
        labs(title = "Pearson Correlation")
    )
  }

  pearson_df <- map_dfr(seq_len(nrow(valid)), function(i) {
    anc1 <- valid$anc1[i]
    anc2 <- valid$anc2[i]
    cl1  <- valid$cluster_1[i]
    cl2  <- valid$cluster_2[i]

    h1 <- h_list[[anc1]]
    h2 <- h_list[[anc2]]

    common_cols <- intersect(
      setdiff(colnames(h1), "Cluster"),
      setdiff(colnames(h2), "Cluster")
    )

    vec1 <- h1 %>% filter(Cluster == cl1) %>% select(all_of(common_cols)) %>% unlist(use.names = FALSE)
    vec2 <- h2 %>% filter(Cluster == cl2) %>% select(all_of(common_cols)) %>% unlist(use.names = FALSE)

    r <- cor(vec1, vec2, method = "pearson")

    tibble(
      pair_label  = paste(anc1, "-", anc2),
      match_label = paste(cl1, "-", cl2),
      pearson_r   = r
    )
  })

  ggplot(pearson_df, aes(x = pair_label, y = match_label,
                          size = abs(pearson_r), color = pearson_r)) +
    geom_point() +
    scale_color_gradient2(low = "#D94A4A", mid = "white", high = "#4A90D9",
                          midpoint = 0, limits = c(-1, 1),
                          name = "Pearson\nr") +
    scale_size_continuous(range = c(4, 14), limits = c(0, 1), name = "|r|") +
    theme_big_text() +
    theme(
      panel.grid.major = element_blank()
    ) +
    labs(x = "Ancestry Pair", y = "Matched Clusters")
}
