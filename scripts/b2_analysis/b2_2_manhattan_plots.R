#!/usr/bin/env Rscript
#
# B2.2 PheWAS Manhattan Plots
#
# Generates per-cluster Manhattan-style plots of PheWAS results:
#   x-axis = ICD10 phenotype (grouped by disease chapter)
#   y-axis = -log10(p-value)
#   color  = disease chapter
#   shape  = effect direction (up triangle = positive beta, down = negative)
#   Bonferroni threshold line + labels for significant associations
#
# Also produces a combined 2-column multi-panel figure (cowplot).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(cowplot)
  library(RSQLite)
  library(yaml)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
source(file.path(base_dir, "scripts/a1_analysis/figure_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
hard_assignment <- "--hard-assignment" %in% args

# ── Config ───────────────────────────────────────────────────────────────────

cfg       <- read_yaml(file.path(base_dir, "config/b2_2_config.yaml"))
label_map <- unlist(cfg$cluster_labels)

if (hard_assignment) {
  results_file <- file.path(base_dir, "results/b2_2_analysis/prs_hard_assignment/phewas_results.csv")
  out_dir      <- file.path(base_dir, "results/b2_2_analysis/prs_hard_assignment")
  message("  ** HARD ASSIGNMENT MODE **")
} else {
  results_file <- file.path(base_dir, "results/b2_2_analysis/phewas_results.csv")
  out_dir      <- file.path(base_dir, "results/b2_2_analysis")
}

chapter_labels <- c(
  E = "Endocrine",     F = "Mental",       G = "Nervous",
  H = "Eye/Ear",       I = "Circulatory",  J = "Respiratory",
  K = "Digestive",     L = "Skin",         M = "Musculoskeletal",
  N = "Genitourinary"
)

# 10-color qualitative palette (colorblind-friendly, from Okabe-Ito + extensions)
chapter_colors <- c(
  E = "#E69F00", F = "#56B4E9", G = "#009E73", H = "#F0E442",
  I = "#CC79A7", J = "#D55E00", K = "#0072B2", L = "#999999",
  M = "#882255", N = "#44AA99"
)

# ── Load ICD10 disease names from UKB database ──────────────────────────────

ukb_db <- file.path(base_dir, "../../../data/ukb/phenotype/ukb18177.db")
con    <- dbConnect(SQLite(), ukb_db)
icd10_names <- dbGetQuery(con, "SELECT value, meaning FROM code WHERE code_id = 19 AND length(value) = 3")
dbDisconnect(con)

# Strip the code prefix from the meaning (e.g., "E11 Non-insulin-dependent..." → "Non-insulin-dependent...")
icd10_names <- icd10_names %>%
  mutate(disease_name = str_trim(str_remove(meaning, "^[A-Z]\\d+\\s+"))) %>%
  select(icd10_code = value, disease_name)

# ── Load and prepare data ────────────────────────────────────────────────────

message("Loading PheWAS results...")
results <- read_csv(results_file, show_col_types = FALSE) %>%
  left_join(icd10_names, by = "icd10_code")

# Remap cluster labels from config (overrides old labels in CSV)
results$cluster_label <- label_map[results$cluster]

bonferroni_threshold <- 0.05 / nrow(results)
bonferroni_y         <- -log10(bonferroni_threshold)

message(sprintf("  %d associations, Bonferroni threshold = %.2e (-log10 = %.2f)",
                nrow(results), bonferroni_threshold, bonferroni_y))

# Compute cumulative x-positions: each ICD10 code gets a unique integer,
# ordered by chapter then code, so chapters cluster together
code_order <- results %>%
  distinct(icd10_code, chapter) %>%
  arrange(chapter, icd10_code) %>%
  mutate(x_pos = row_number())

# Chapter midpoints for x-axis labels
chapter_midpoints <- code_order %>%
  group_by(chapter) %>%
  summarise(center = mean(x_pos), .groups = "drop")

# Merge positions back
plot_data <- results %>%
  left_join(code_order, by = c("icd10_code", "chapter")) %>%
  mutate(
    log_p     = -log10(pmax(p_value, 1e-300)),
    direction = if_else(beta >= 0, "Positive", "Negative"),
    chapter_name = chapter_labels[chapter]
  )

# ── Plotting function ────────────────────────────────────────────────────────

make_manhattan <- function(df, cluster_id, cluster_label,
                           base_size = 14, show_legend = TRUE,
                           size_by_y = FALSE) {

  cluster_df <- df %>% filter(cluster == cluster_id)
  # Label only Bonferroni-significant associations (points above the dashed line)
  sig_df     <- cluster_df %>% filter(log_p > bonferroni_y)

  # Optionally scale triangle size with the y-axis value (-log10 P): stronger
  # associations are drawn as larger markers. Size legend is suppressed since
  # marker size is redundant with the y position.
  if (size_by_y) {
    point_layer <- geom_point(aes(size = log_p), alpha = 0.7)
    size_scale  <- scale_size_continuous(range = c(1, 6), guide = "none")
  } else {
    point_layer <- geom_point(size = 2, alpha = 0.7)
    size_scale  <- NULL
  }

  p <- ggplot(cluster_df, aes(x = x_pos, y = log_p,
                               color = chapter, shape = direction)) +
    point_layer +
    size_scale +
    geom_hline(yintercept = bonferroni_y, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    scale_color_manual(
      values = chapter_colors,
      labels = chapter_labels,
      name   = "Disease Category"
    ) +
    scale_shape_manual(
      values = c("Positive" = 24, "Negative" = 25),
      name   = "Effect Direction"
    ) +
    scale_x_continuous(
      breaks = chapter_midpoints$center,
      labels = chapter_labels[chapter_midpoints$chapter],
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    coord_cartesian(ylim = c(0, max(cluster_df$log_p, bonferroni_y) * 1.15)) +
    labs(
      x     = NULL,
      y     = expression(-log[10](italic(P))),
      title = cluster_label
    ) +
    theme_big_text(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = base_size * 0.7),
      plot.title          = element_text(size = base_size * 1.1, face = "bold")
    )

  # Label Bonferroni-significant points
  if (nrow(sig_df) > 0) {
    p <- p +
      geom_text_repel(
        data        = sig_df,
        aes(x = x_pos, y = log_p, label = disease_name, color = chapter),
        size        = base_size / 3.5,
        max.overlaps = 30,
        segment.size = 0.3,
        box.padding  = 0.4,
        show.legend  = FALSE,
        inherit.aes  = FALSE
      )
  }

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}

# ── Generate individual plots ────────────────────────────────────────────────

# Desired META cluster display order (Glycemic -> Lpa); clusters not in the list
# fall back to natural K order at the end.
desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
clusters <- results %>%
  distinct(cluster, cluster_label) %>%
  mutate(.ord = match(cluster, desired_k_order),
         .ord = ifelse(is.na(.ord),
                       100 + as.numeric(str_extract(cluster, "\\d+")), .ord)) %>%
  arrange(.ord) %>%
  select(-.ord)

message("Generating individual Manhattan plots...")

for (i in seq_len(nrow(clusters))) {
  cl  <- clusters$cluster[i]
  lab <- clusters$cluster_label[i]

  p <- make_manhattan(plot_data, cl, lab, base_size = 14, show_legend = TRUE)

  out_file <- file.path(out_dir, sprintf("manhattan_%s.png", cl))
  ggsave(out_file, p, width = 16, height = 6, dpi = 300)
  message(sprintf("  Saved %s", basename(out_file)))
}

# ── Combined multi-panel figure ──────────────────────────────────────────────

message("Generating combined multi-panel figure...")

# Build panel plots without individual legends (larger text + size-scaled markers)
panel_plots <- map2(clusters$cluster, clusters$cluster_label, function(cl, lab) {
  make_manhattan(plot_data, cl, lab, base_size = 18, show_legend = FALSE,
                 size_by_y = TRUE)
})

# Extract shared legend from a plot with legend enabled
legend_plot <- make_manhattan(plot_data, clusters$cluster[1],
                              clusters$cluster_label[1],
                              base_size = 18, show_legend = TRUE,
                              size_by_y = TRUE) +
  theme(legend.position = "bottom",
        legend.box = "horizontal") +
  guides(
    color = guide_legend(nrow = 2, override.aes = list(size = 3)),
    shape = guide_legend(nrow = 1, override.aes = list(size = 3))
  )
shared_legend <- get_legend(legend_plot)

# Layout: 2-column (left: K1-K5, right: K6-K10)
n_panels <- length(panel_plots)
n_left <- ceiling(n_panels / 2)
n_right <- n_panels - n_left

left_col  <- plot_grid(plotlist = panel_plots[1:n_left], ncol = 1,
                       labels = letters[1:n_left], label_size = 18)
right_plots <- panel_plots[(n_left + 1):n_panels]
if (n_right < n_left) right_plots <- c(right_plots, list(NULL))
right_col <- plot_grid(plotlist = right_plots, ncol = 1,
                       labels = c(letters[(n_left + 1):n_panels],
                                  if (n_right < n_left) "" else NULL),
                       label_size = 18)
panels    <- plot_grid(left_col, right_col, ncol = 2)

combined <- plot_grid(panels, shared_legend, ncol = 1, rel_heights = c(1, 0.03))

ggsave(file.path(out_dir, "manhattan_combined.png"), combined,
       width = 20, height = 28, dpi = 600)
message("  Saved manhattan_combined.png")

# ── 2x2 selected cluster figure ──────────────────────────────────────────────

message("Generating 2x2 selected cluster figure...")

# 2x2 layout (row-major): TL=Lpa(K1), TR=Glycemic(K10), BL=Blood Pressure-Stature(K5), BR=Adiponectin(K2)
selected_clusters <- c("K1", "K10", "K5", "K2")
panel_labels <- c("a", "b", "c", "d")
selected_plots <- map2(selected_clusters, panel_labels, function(cl, lbl) {
  make_manhattan(plot_data, cl, label_map[cl], base_size = 20, show_legend = FALSE,
                 size_by_y = TRUE) +
    labs(tag = lbl) +
    theme(plot.tag = element_text(size = 22, face = "bold"),
          plot.tag.position = c(0.02, 0.98))
})

legend_2x2 <- get_legend(
  make_manhattan(plot_data, "K4", label_map["K4"], base_size = 20, show_legend = TRUE,
                 size_by_y = TRUE) +
    theme(legend.position = "bottom", legend.box = "horizontal") +
    guides(
      color = guide_legend(nrow = 2, override.aes = list(size = 3)),
      shape = guide_legend(nrow = 1, override.aes = list(size = 3))
    )
)

grid_2x2 <- plot_grid(
  selected_plots[[1]], selected_plots[[2]],
  selected_plots[[3]], selected_plots[[4]],
  ncol = 2, nrow = 2
)

fig_2x2 <- plot_grid(grid_2x2, legend_2x2, ncol = 1, rel_heights = c(1, 0.06))

ggsave(file.path(out_dir, "manhattan_2x2_selected.png"), fig_2x2,
       width = 20, height = 14, dpi = 600)
message("  Saved manhattan_2x2_selected.png")

message("Done.")
