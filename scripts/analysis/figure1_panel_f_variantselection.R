# figure1_panel_f_variantselection.R
# Panel F: Barplot of the most recurrently top-weighted variants across
# all bNMF clusters and ancestries.

# Assumes figure_utils.R has been sourced.


#' Create Panel F: Variant selection frequency barplot
#'
#' For each ancestry, reads W_matrix.tsv and identifies the top 10 variants
#' (by weight) in each cluster. Counts how many cluster x ancestry combinations
#' each variant appears in, then plots the overall top 10 most recurrent.
#'
#' @param results_dir Path to results directory (default: "results")
#' @param ancestries Character vector of ancestry codes
#' @return ggplot object
plot_panel_f <- function(results_dir = "results",
                         ancestries = c("EUR", "AFR", "META")) {

  top_n <- 10

  # Collect top-10 variants per cluster per ancestry
  top_entries <- map_dfr(ancestries, function(anc) {
    w_path <- file.path(results_dir, anc, sprintf("W_matrix_%s.tsv", anc))
    if (!file.exists(w_path)) {
      warning(sprintf("W_matrix_%s.tsv not found for %s at %s", anc, anc, w_path))
      return(NULL)
    }
    w_df <- read_tsv(w_path, show_col_types = FALSE)
    cluster_cols <- setdiff(colnames(w_df), "VAR_ID")

    map_dfr(cluster_cols, function(kcol) {
      w_df %>%
        arrange(desc(.data[[kcol]])) %>%
        slice_head(n = min(top_n, nrow(w_df))) %>%
        select(VAR_ID) %>%
        mutate(ancestry = anc, cluster = kcol)
    })
  })

  if (is.null(top_entries) || nrow(top_entries) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No W matrices found",
                 size = 6) +
        theme_big_text() +
        labs(title = "Variant Selection Frequency")
    )
  }

  # Count appearances per variant and take overall top 10
  var_counts <- top_entries %>%
    count(VAR_ID, name = "times_selected") %>%
    arrange(desc(times_selected), VAR_ID) %>%
    slice_head(n = min(top_n, n())) %>%
    mutate(VAR_ID = fct_inorder(VAR_ID))

  ggplot(var_counts, aes(x = VAR_ID, y = times_selected)) +
    geom_col(fill = "gray50", width = 0.7) +
    scale_y_continuous(breaks = seq(0, max(var_counts$times_selected))) +
    theme_big_text() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = "Variant", y = "Times Selected (Top 10)")
}
