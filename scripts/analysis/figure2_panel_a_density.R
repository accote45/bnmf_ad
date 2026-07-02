# figure2_panel_a_density.R
# Panel A: Density curves for each cluster PRS + GRS total,
# faceted by individual predicted ancestry (ancestry_pred).
#
# Each GWAS ancestry gets a distinct base color (blue for EUR, red for AFR,
# green for META). Clusters within an ancestry are rendered as shades of that
# base color with distinct linetypes so they remain easy to tell apart.
#
# Assumes figure_utils.R and figure2_utils.R have been sourced.


#' Create Panel A: PRS density curves faceted by individual ancestry
#'
#' @param dat             tibble with standardized PRS columns and ancestry_pred
#' @param prs_meta        tibble with PRS column metadata
#' @param ancestries      Character vector of ancestry_pred values to facet by
#'                        (individual ancestry, lowercase: "eur", "afr")
#' @param gwas_ancestries Character vector of GWAS ancestries whose PRS columns
#'                        to display (uppercase: "EUR", "AFR", "META")
#' @return ggplot object
plot_panel_a_density <- function(dat, prs_meta,
    ancestries = c("eur", "afr"),
    gwas_ancestries = c("EUR", "AFR")) {

  # Filter prs_meta to requested GWAS ancestries (keep clusters + GRS)
  plot_meta <- prs_meta %>%
    filter(gwas_ancestry %in% gwas_ancestries)

  if (nrow(plot_meta) == 0) {
    stop("No PRS columns match the requested gwas_ancestries: ",
         paste(gwas_ancestries, collapse = ", "))
  }

  # Filter to requested individual ancestries
  dat_sub <- dat %>% filter(ancestry_pred %in% ancestries)
  cat(sprintf("  Panel A: %d individuals after ancestry filter (%s)\n",
    nrow(dat_sub), paste(ancestries, collapse = ", ")))
  cat(sprintf("  Panel A: showing GWAS ancestries: %s\n",
    paste(gwas_ancestries, collapse = ", ")))

  # Re-standardize PRS within each individual ancestry group
  prs_cols <- plot_meta$prs_col
  dat_sub <- dat_sub %>%
    group_by(ancestry_pred) %>%
    mutate(across(all_of(prs_cols), ~ {
      if (all(is.na(.x)) || sd(.x, na.rm = TRUE) == 0) .x
      else (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)
    })) %>%
    ungroup()

  # Pivot to long format
  long <- dat_sub %>%
    select(person_id, ancestry_pred, all_of(prs_cols)) %>%
    pivot_longer(cols = all_of(prs_cols), names_to = "prs_col", values_to = "prs_value") %>%
    filter_out(is.na(prs_value)) %>%
    left_join(plot_meta %>% select(prs_col, label, gwas_ancestry, is_grs), by = "prs_col")

  # --- Build color + linetype palettes ---
  # Each GWAS ancestry gets shades of a base hue (dark → light).
  # Clusters are ordered first (K1, K2, …), then GRS total last.
  ancestry_shades <- list(
    EUR  = c("#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1"),
    AFR  = c("#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1"),
    META = c("#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0")
  )

  # Linetypes: clusters get solid/dashed/dotdash, GRS total always dotted
  cluster_linetypes <- c("solid", "dashed", "dotdash", "longdash", "twodash")

  color_map <- c()
  ltype_map <- c()
  label_order <- character(0)

  for (anc in gwas_ancestries) {
    anc_meta <- plot_meta %>%
      filter(gwas_ancestry == anc) %>%
      arrange(is_grs, prs_col)   # clusters first, GRS last

    shades <- ancestry_shades[[anc]]
    if (is.null(shades)) {
      # Fallback grey palette for unexpected ancestry codes
      shades <- c("#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9")
    }

    for (i in seq_len(nrow(anc_meta))) {
      lbl <- anc_meta$label[i]
      color_map[lbl] <- shades[min(i, length(shades))]

      if (anc_meta$is_grs[i]) {
        ltype_map[lbl] <- "dotted"
      } else {
        # K1 → solid, K2 → dashed, K3 → dotdash, …
        k_idx <- sum(!anc_meta$is_grs[seq_len(i)])
        ltype_map[lbl] <- cluster_linetypes[min(k_idx, length(cluster_linetypes))]
      }

      label_order <- c(label_order, lbl)
    }
  }

  # Drop labels with no data
  label_order <- intersect(label_order, unique(long$label))
  long <- long %>%
    filter(label %in% label_order) %>%
    mutate(label = factor(label, levels = label_order))

  # Nice facet labels
  anc_labels <- c(
    "eur" = "European", "afr" = "African", "amr" = "Admixed American",
    "eas" = "East Asian", "sas" = "South Asian"
  )
  long <- long %>%
    mutate(ancestry_label = factor(
      if_else(ancestry_pred %in% names(anc_labels),
              anc_labels[ancestry_pred],
              toupper(ancestry_pred)),
      levels = anc_labels[intersect(ancestries, names(anc_labels))]
    ))

  ggplot(long, aes(x = prs_value, color = label, linetype = label)) +
    geom_density(linewidth = 0.8) +
    facet_wrap(~ ancestry_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = color_map[label_order]) +
    scale_linetype_manual(values = ltype_map[label_order]) +
    guides(color = guide_legend(title = "Score"),
           linetype = guide_legend(title = "Score")) +
    theme_big_text() +
    labs(
      x        = "Standardized PRS",
      y        = "Density"
    )
}
