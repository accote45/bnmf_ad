# figure2_panel_c_forest.R
# Panel C: Forest plot of PRS associations (OR per SD) for T2D, CAD, and
# Confluence (CAD or T2D), arranged as 3 rows x 1 column.
#
# Covariates: age + age^2 + sex + PC1-PC16.
# Color by GWAS ancestry.
#
# Assumes figure_utils.R and figure2_utils.R have been sourced.

library(patchwork)


#' Create Panel C: Forest plot (3 outcomes stacked)
#'
#' @param dat              tibble with standardized PRS columns and outcomes
#' @param prs_meta         tibble with PRS column metadata
#' @param outcomes         Named character vector (display name = column name)
#' @param ancestry_filter  Optional character vector of ancestry_pred values to
#'                         subset individuals (e.g., "eur", "afr"). NULL = all.
#' @return ggplot object (wrapped patchwork)
plot_panel_c_forest <- function(dat, prs_meta,
    outcomes = c("Type 2 Diabetes" = "t2d",
                 "CAD" = "cad",
                 "Confluence" = "confluence"),
    ancestry_filter = NULL) {

  # Optionally filter to specific individual ancestries
  if (!is.null(ancestry_filter)) {
    dat <- dat %>% filter(ancestry_pred %in% ancestry_filter)
    cat(sprintf("  Panel C: filtered to %d individuals (%s)\n",
      nrow(dat), paste(ancestry_filter, collapse = ", ")))
  }

  cat("Running logistic regressions for forest plot...\n")
  results <- run_all_regressions(dat, prs_meta, outcomes)

  # Drop clusters with no valid results (non-existent clusters)
  results <- results %>% filter_out(is.na(OR))

  # Order labels: group by ancestry, within each ancestry clusters first then GRS
  # Only include labels that have valid results
  label_order <- prs_meta %>% arrange(gwas_ancestry, is_grs, prs_col) %>% pull(label)
  label_order <- intersect(label_order, unique(results$label))
  results <- results %>%
    mutate(label = factor(label, levels = rev(label_order)))

  # Ancestry colors
  anc_colors <- c("EUR" = "#4A90D9", "AFR" = "#D94A4A", "META" = "#50B878")

  # Outcome display order
  results <- results %>%
    mutate(outcome_label = factor(outcome_label, levels = names(outcomes)))

  # Build one forest plot per outcome
  plot_list <- map(names(outcomes), function(olabel) {
    sub <- results %>% filter(outcome_label == olabel)

    p <- ggplot(sub, aes(x = OR, y = label, color = gwas_ancestry)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25) +
      scale_color_manual(values = anc_colors, name = "GWAS\nAncestry") +
      theme_big_text() +
      labs(
        x = "OR per SD (95% CI)",
        y = NULL,
        title = olabel
      ) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = if (olabel == tail(names(outcomes), 1)) "right" else "none"
      )

    p
  })

  # Stack 3 forest plots vertically and wrap as a single element
  # so patchwork assigns one "C" tag to the group
  stacked <- wrap_plots(plot_list, ncol = 1)
  wrap_elements(stacked)
}
