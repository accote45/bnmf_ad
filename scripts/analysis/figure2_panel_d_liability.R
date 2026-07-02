# figure2_panel_d_liability.R
# Panel D: Bar plot of R-squared for each PRS cluster,
# faceted by outcome (T2D, CAD, Confluence).
#
# Supports two metrics:
#   - "liability" (default): Lee et al. 2012 liability-scale R²
#   - "nagelkerke": incremental Nagelkerke pseudo-R²
#
# Assumes figure_utils.R and figure2_utils.R have been sourced.


#' Create Panel D: R-squared bar plot
#'
#' @param dat       tibble with standardized PRS columns and outcomes
#' @param prs_meta  tibble with PRS column metadata
#' @param outcomes  Named character vector (display name = column name)
#' @param metric    Which R² metric to plot: "liability" or "nagelkerke"
#' @return ggplot object
plot_panel_d_liability <- function(dat, prs_meta,
    outcomes = c("Type 2 Diabetes" = "t2d",
                 "CAD" = "cad",
                 "Confluence" = "confluence"),
    metric = c("liability", "nagelkerke")) {

  metric <- match.arg(metric)

  cat(sprintf("Computing R-squared (metric = %s)...\n", metric))
  r2_results <- run_all_liability_r2(dat, prs_meta, outcomes)

  # Drop clusters with no valid results (non-existent clusters)
  if (metric == "liability") {
    r2_results <- r2_results %>% filter_out(is.na(r2_liability))
  } else {
    r2_results <- r2_results %>% filter_out(is.na(r2_obs))
  }

  # Order labels: cluster PRS by ancestry then K, then GRS
  # Only include labels that have valid results
  label_order <- c(
    prs_meta %>% filter_out(is_grs) %>% arrange(gwas_ancestry, prs_col) %>% pull(label),
    prs_meta %>% filter(is_grs) %>% arrange(gwas_ancestry) %>% pull(label)
  )
  label_order <- intersect(label_order, unique(r2_results$label))
  r2_results <- r2_results %>%
    mutate(
      label = factor(label, levels = label_order),
      outcome_label = factor(outcome_label, levels = names(outcomes))
    )

  # Select y-axis column and label based on metric
  if (metric == "liability") {
    y_col <- "r2_liability"
    y_lab <- expression(Liability ~ R^2)
  } else {
    y_col <- "r2_obs"
    y_lab <- expression(Nagelkerke ~ R^2)
  }

  # Ancestry colors
  anc_colors <- c("EUR" = "#4A90D9", "AFR" = "#D94A4A", "META" = "#50B878")

  ggplot(r2_results, aes(x = label, y = .data[[y_col]], fill = gwas_ancestry)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~ outcome_label, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = anc_colors, name = "GWAS\nAncestry") +
    theme_big_text() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = "Cluster",
      y = y_lab
    )
}
