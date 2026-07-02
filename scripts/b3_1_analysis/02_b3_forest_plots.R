#!/usr/bin/env Rscript
# 02_b3_forest_plots.R
# Per-cluster forest plots and summary heatmap for B3 biomarker associations.
# Forest plots show All/Male/Female betas with trait category grouping.
#
# Usage:
#   Rscript scripts/b3_1_analysis/02_b3_forest_plots.R
#   Rscript scripts/b3_1_analysis/02_b3_forest_plots.R --config config/b3_1_config.yaml

library(data.table)
library(ggplot2)
library(ggh4x)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b3_1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

hard_assignment <- "--hard-assignment" %in% args
selected_clusters <- NULL
if ("--selected-clusters" %in% args) {
  selected_clusters <- trimws(strsplit(args[which(args == "--selected-clusters") + 1], ",")[[1]])
}

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir

if (hard_assignment) {
  results_dir <- "results/b3_1_analysis/prs_hard_assignment"
  cat("  ** HARD ASSIGNMENT MODE **\n")
  cat(sprintf("  Results dir: %s\n\n", results_dir))
}

cat("=== B3 Step 2: Forest Plots ===\n")

# ============================================================
# 1. Load results
# ============================================================
cat("\n--- Loading results ---\n")
all_results <- fread(file.path(results_dir, "biomarker_association_all.csv"))
indiv <- all_results[model_type == "individual"]
cat(sprintf("  Individual model results: %d rows\n", nrow(indiv)))

# ============================================================
# 2. Apply cluster labels
# ============================================================
cluster_labels <- unlist(cfg$cluster_labels)

indiv[, cluster := gsub("PRS_", "", predictor)]
for (k in names(cluster_labels)) {
  indiv[cluster == k, cluster_label := cluster_labels[[k]]]
}

# Facet column order: desired META display order (Glycemic -> Lpa); any clusters
# not in the list fall back to natural K order at the end.
all_k <- paste0("K", 1:length(cluster_labels))
desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
k_levels <- c(desired_k_order[desired_k_order %in% all_k],
              setdiff(all_k, desired_k_order))
label_levels <- unname(cluster_labels[k_levels])
indiv[, cluster_label := factor(cluster_label, levels = label_levels)]

# Trait ordering within categories (bottom to top per panel)
trait_order <- c("HbA1c", "Fasting Glucose",
                 "Lpa", "non-HDL Cholesterol", "Triglycerides", "LDL",
                 "WHtR", "Waist Circumference", "BMI",
                 "CRP", "eGFR", "SBP")
indiv[, trait_display := factor(trait_display, levels = rev(trait_order))]

# Trait category assignment
category_map <- c(
  "HbA1c" = "Glycemic", "Fasting Glucose" = "Glycemic",
  "Lpa" = "Lipid", "non-HDL Cholesterol" = "Lipid",
  "Triglycerides" = "Lipid", "LDL" = "Lipid",
  "WHtR" = "Anthropometric", "Waist Circumference" = "Anthropometric",
  "BMI" = "Anthropometric",
  "CRP" = "Other", "eGFR" = "Other", "SBP" = "Other"
)
category_order <- c("Glycemic", "Lipid", "Anthropometric", "Other")
category_colors <- c("Glycemic" = "#E67E22", "Lipid" = "#2980B9",
                      "Anthropometric" = "#27AE60", "Other" = "#8E44AD")

category_strips <- strip_themed(
  background_y = lapply(category_colors[category_order], function(col) {
    element_rect(fill = col, color = NA)
  })
)
indiv[, trait_category := factor(category_map[as.character(trait_display)],
                                  levels = category_order)]

# Exclude PREVENT-ASCVD from plots
indiv <- indiv[!is.na(trait_category)]

# Sex group factor
indiv[, sex_group := factor(sex_group, levels = c("All", "Male", "Female"))]

# Significance stars (for All group annotation only)
indiv[, sig_stars := fifelse(fdr_q < 0.001, "***",
                    fifelse(fdr_q < 0.01, "**",
                    fifelse(fdr_q < 0.05, "*",
                    fifelse(p_value < 0.05, "†", ""))))]

sex_colors <- c("All" = "#2C3E50", "Male" = "#1B9E77", "Female" = "#D95F02")

# Per-cluster colors, keyed by cluster label. Mirrors all_colors in
# scripts/b2_analysis/b2_1_analysis.R (the cumhaz_top33 figures) so the selected
# clusters keep the same color identity across the B-series figures.
cluster_palette <- c(
  "Lpa" = "#E31A1C", "Adiponectin" = "#FDBF6F", "Platelet" = "#33A02C",
  "SHBG" = "#6A3D9A", "Blood Pressure-Stature" = "#FF7F00", "Metabolic" = "#FB9A99",
  "Triglycerides-HDL" = "#A6CEE3", "ALP-LDL" = "#B2DF8A", "Obesity" = "#B15928",
  "Glycemic" = "#1F78B4"
)

# ============================================================
# 3. Per-cluster forest plots
# ============================================================
cat("\n--- Creating per-cluster forest plots ---\n")

for (k in k_levels) {
  k_label <- cluster_labels[[k]]
  k_data <- indiv[cluster == k]

  if (nrow(k_data) == 0) next

  # Stars only on All group
  k_data[, star_label := fifelse(sex_group == "All", sig_stars, "")]

  p <- ggplot(k_data, aes(x = beta, y = trait_display, color = sex_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   position = position_dodge(width = 0.6),
                   height = 0.2, linewidth = 0.6) +
    geom_text(aes(label = star_label),
              position = position_dodge(width = 0.6),
              vjust = -0.8, size = 3.5, show.legend = FALSE) +
    scale_color_manual(values = sex_colors, name = "Sex") +
    facet_grid2(trait_category ~ ., scales = "free_y", space = "free_y",
                switch = "y", strip = category_strips) +
    labs(
      x = "Beta per SD PRS (95% CI)",
      y = NULL,
      title = sprintf("Cluster %s: %s", k, k_label),
      subtitle = sprintf("EUR validation (N range: %s–%s)",
                         format(min(k_data$n_obs), big.mark = ","),
                         format(max(k_data$n_obs), big.mark = ","))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, color = "grey40"),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 10, color = "white"),
      panel.spacing = unit(0.3, "lines")
    )

  safe_label <- gsub("[^A-Za-z0-9]", "", k_label)
  out_file <- file.path(results_dir, sprintf("forest_%s_%s.png", k, safe_label))
  ggsave(out_file, p, width = 10, height = 8, dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(out_file)))
}

# ============================================================
# 4. Combined multi-panel forest plot
# ============================================================
cat("\n--- Creating combined forest plot ---\n")

indiv[, star_label := fifelse(sex_group == "All", sig_stars, "")]

p_combined <- ggplot(indiv, aes(x = beta, y = trait_display, color = sex_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.6), size = 1.5) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 position = position_dodge(width = 0.6),
                 height = 0.15, linewidth = 0.4) +
  scale_color_manual(values = sex_colors, name = "Sex") +
  facet_grid2(trait_category ~ cluster_label,
              scales = "free", space = "free_y",
              switch = "y", strip = category_strips) +
  labs(
    x = "Beta per SD PRS (95% CI)",
    y = NULL,
    title = "Cluster PRS Associations with Continuous Biomarker Traits"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    strip.text.x = element_text(face = "bold", size = 9),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 9, color = "white"),
    strip.placement = "outside",
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.3, "lines"),
    panel.spacing.y = unit(0.3, "lines")
  )

ggsave(file.path(results_dir, "forest_combined.png"),
       p_combined, width = 20, height = 10, dpi = 300)
cat("  Saved: forest_combined.png\n")

# ============================================================
# 4b. Selected-clusters combined forest plot (optional)
# ============================================================
if (!is.null(selected_clusters)) {
  cat(sprintf("\n--- Creating selected-clusters forest plot: %s ---\n",
              paste(selected_clusters, collapse = ", ")))

  # Shared theme for both selected-cluster figures (single column + sex-split)
  selected_theme <- theme_minimal(base_size = 16) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 14, color = "grey40"),
      strip.text.x = element_text(face = "bold", size = 14),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 14, color = "white"),
      strip.placement = "outside",
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      panel.grid.major.y = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.3, "lines"),
      panel.spacing.y = unit(0.3, "lines")
    )

  # --- Main figure: single column, "All" estimate only, colored by cluster ---
  indiv_sel <- indiv[cluster_label %in% selected_clusters & sex_group == "All"]
  indiv_sel[, cluster_label := factor(cluster_label, levels = selected_clusters)]

  p_selected <- ggplot(indiv_sel, aes(x = beta, y = trait_display, color = cluster_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.6), size = 4) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   position = position_dodge(width = 0.6),
                   height = 0.3, linewidth = 1.1) +
    scale_color_manual(values = cluster_palette, name = "cPGS") +
    facet_grid2(trait_category ~ ., scales = "free_y", space = "free_y",
                switch = "y", strip = category_strips) +
    labs(
      x = "cPGS Beta/SD",
      y = NULL,
      # Wrapped to two lines so the title fits the narrow single-column width
      title = "Cluster PRS Associations with\nContinuous Biomarker Traits"
    ) +
    selected_theme

  ggsave(file.path(results_dir, "forest_combined_selected.png"),
         p_selected, width = 9, height = 12, dpi = 300)
  cat("  Saved: forest_combined_selected.png\n")

  # --- Companion figure: sex-stratified (Male | Female), faceted by sex ---
  indiv_sex <- indiv[cluster_label %in% selected_clusters &
                       sex_group %in% c("Male", "Female")]
  indiv_sex[, cluster_label := factor(cluster_label, levels = selected_clusters)]
  indiv_sex[, sex_group := factor(sex_group, levels = c("Male", "Female"))]

  p_sex <- ggplot(indiv_sex, aes(x = beta, y = trait_display, color = cluster_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.6), size = 4) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   position = position_dodge(width = 0.6),
                   height = 0.3, linewidth = 1.1) +
    scale_color_manual(values = cluster_palette, name = "cPGS") +
    facet_grid2(trait_category ~ sex_group, scales = "free", space = "free_y",
                switch = "y", strip = category_strips) +
    labs(
      x = "cPGS Beta/SD",
      y = NULL,
      title = "Cluster PRS Associations with Continuous Biomarker Traits",
      subtitle = "Sex-stratified (Male vs Female)"
    ) +
    selected_theme

  ggsave(file.path(results_dir, "forest_combined_selected_sex.png"),
         p_sex, width = 13, height = 12, dpi = 300)
  cat("  Saved: forest_combined_selected_sex.png\n")
}

# ============================================================
# 5. Summary heatmap (clusters x traits, All only)
# ============================================================
cat("\n--- Creating summary heatmap ---\n")

hm_data <- indiv[sex_group == "All",
                 .(cluster_label, trait_display, beta, fdr_q, p_value, sig_stars)]

beta_max <- max(abs(hm_data$beta), na.rm = TRUE)

p_heatmap <- ggplot(hm_data, aes(x = cluster_label, y = trait_display)) +
  geom_tile(aes(fill = beta), color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_stars), size = 4, color = "black") +
  scale_fill_gradient2(
    low = "#3498DB", mid = "white", high = "#E74C3C",
    midpoint = 0, limits = c(-beta_max, beta_max),
    name = "Beta per SD"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Cluster PRS vs Biomarker Associations",
    subtitle = "* FDR<0.05  ** FDR<0.01  *** FDR<0.001  † nominal P<0.05"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 0.7)

ggsave(file.path(results_dir, "heatmap_cluster_biomarker.png"),
       p_heatmap, width = 12, height = 8, dpi = 300)
cat("  Saved: heatmap_cluster_biomarker.png\n")

# ============================================================
# 6. Individual vs Joint model comparison (All only)
# ============================================================
cat("\n--- Creating individual vs joint comparison ---\n")

joint <- all_results[model_type == "joint" & sex_group == "All"]
joint[, cluster := gsub("PRS_", "", predictor)]
for (k in names(cluster_labels)) {
  joint[cluster == k, cluster_label := cluster_labels[[k]]]
}
joint[, cluster_label := factor(cluster_label, levels = label_levels)]
joint[, trait_display := factor(trait_display, levels = rev(trait_order))]

indiv_all <- indiv[sex_group == "All"]

comparison <- merge(
  indiv_all[, .(trait, predictor, beta_indiv = beta,
            ci_lower_indiv = ci_lower, ci_upper_indiv = ci_upper,
            fdr_q_indiv = fdr_q)],
  joint[, .(trait, predictor, beta_joint = beta,
            ci_lower_joint = ci_lower, ci_upper_joint = ci_upper,
            fdr_q_joint = fdr_q)],
  by = c("trait", "predictor")
)

comparison_sig <- comparison[fdr_q_indiv < 0.05]

if (nrow(comparison_sig) > 0) {
  comparison_sig[, cluster := gsub("PRS_", "", predictor)]
  for (k in names(cluster_labels)) {
    comparison_sig[cluster == k, cluster_label := cluster_labels[[k]]]
  }

  comp_long <- rbindlist(list(
    comparison_sig[, .(trait, cluster_label, model = "Individual",
                       beta = beta_indiv, ci_lower = ci_lower_indiv,
                       ci_upper = ci_upper_indiv)],
    comparison_sig[, .(trait, cluster_label, model = "Joint",
                       beta = beta_joint, ci_lower = ci_lower_joint,
                       ci_upper = ci_upper_joint)]
  ))
  comp_long[, label := paste(cluster_label, trait, sep = " : ")]
  comp_long[, label := factor(label, levels = rev(unique(label)))]

  p_comp <- ggplot(comp_long, aes(x = beta, y = label, color = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.5), size = 2.5) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   position = position_dodge(width = 0.5),
                   height = 0.3, linewidth = 0.6) +
    scale_color_manual(values = c("Individual" = "#E74C3C", "Joint" = "#3498DB"),
                       name = "Model") +
    labs(
      x = "Beta per SD PRS (95% CI)",
      y = NULL,
      title = "Individual vs Joint Model: FDR-Significant Associations"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 15),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey92"),
      panel.grid.minor = element_blank()
    )

  n_sig <- length(unique(comp_long$label))
  plot_height <- max(6, 0.5 * n_sig + 2)

  ggsave(file.path(results_dir, "forest_individual_vs_joint.png"),
         p_comp, width = 12, height = plot_height, dpi = 300)
  cat("  Saved: forest_individual_vs_joint.png\n")
} else {
  cat("  No FDR-significant individual results — skipping comparison plot.\n")
}

cat("\n=== Done ===\n")
