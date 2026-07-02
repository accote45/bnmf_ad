# figure2_panel_b_corr.R
# Panel B: Pearson correlation heatmap of all PRS columns (including GRS total).
#
# Assumes figure_utils.R and figure2_utils.R have been sourced.


#' Create Panel B: Pearson correlation heatmap
#'
#' @param dat              tibble with standardized PRS columns
#' @param prs_meta         tibble with PRS column metadata
#' @param ancestry_filter  Optional character vector of ancestry_pred values to
#'                         subset individuals (e.g., "eur", "afr"). NULL = all.
#' @return ggplot object
plot_panel_b_corr <- function(dat, prs_meta, ancestry_filter = NULL) {

  # Optionally filter to specific individual ancestries
  if (!is.null(ancestry_filter)) {
    dat <- dat %>% filter(ancestry_pred %in% ancestry_filter)
    cat(sprintf("  Panel B: filtered to %d individuals (%s)\n",
      nrow(dat), paste(ancestry_filter, collapse = ", ")))
  }

  # Order: cluster PRS first (by ancestry, then K), then GRS total
  meta_ordered <- bind_rows(
    prs_meta %>% filter_out(is_grs) %>% arrange(gwas_ancestry, prs_col),
    prs_meta %>% filter(is_grs) %>% arrange(gwas_ancestry)
  )

  prs_cols <- meta_ordered$prs_col
  labels   <- meta_ordered$label

  # Drop PRS columns with no data in this subset
  has_data <- map_lgl(prs_cols, ~ !all(is.na(dat[[.x]])))
  prs_cols <- prs_cols[has_data]
  labels   <- labels[has_data]

  # Extract matrix and compute Pearson correlation
  mat <- dat %>% select(all_of(prs_cols)) %>% as.matrix()
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")

  # Assign display labels
  rownames(cor_mat) <- labels
  colnames(cor_mat) <- labels

  # Melt to long format for ggplot
  cor_long <- cor_mat %>%
    as_tibble(rownames = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
    mutate(
      var1 = factor(var1, levels = labels),
      var2 = factor(var2, levels = labels)
    )

  ggplot(cor_long, aes(x = var2, y = var1, fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "#D94A4A", mid = "white", high = "#4A90D9",
      midpoint = 0, limits = c(-1, 1),
      name = "Pearson r"
    ) +
    theme_big_text() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12)
    ) +
    labs(x = NULL, y = NULL) +
    coord_fixed()
}
