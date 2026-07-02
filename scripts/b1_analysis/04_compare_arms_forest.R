#!/usr/bin/env Rscript
# 04_compare_arms_forest.R
# Side-by-side forest plots comparing two B1 arms.
# Loads association_results_all.csv + genome_wide_prs_results.csv from each arm.
# Cluster labels are read from each arm's config YAML.
#
# Usage:
#   Rscript scripts/b1_analysis/04_compare_arms_forest.R \
#     --arm1-dir results/b1_comparison/core_published \
#     --arm2-dir results/b1_comparison/published \
#     --arm1-config config/b1_comparison_core_published.yaml \
#     --arm2-config config/b1_comparison_published.yaml \
#     --arm1-name "Core Published" --arm2-name "Published" \
#     --output-dir results/b1_comparison/figures_core_vs_full_published

library(data.table)
library(ggplot2)
library(yaml)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0) args[idx + 1] else default
}

arm1_dir    <- get_arg("--arm1-dir", "results/b1_comparison/splitukb")
arm2_dir    <- get_arg("--arm2-dir", "results/b1_comparison/published")
arm1_config <- get_arg("--arm1-config")
arm2_config <- get_arg("--arm2-config")
arm1_name   <- get_arg("--arm1-name", "Arm 1")
arm2_name   <- get_arg("--arm2-name", "Arm 2")
output_dir  <- get_arg("--output-dir", "results/b1_comparison/figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load cluster labels from configs ---
load_labels <- function(config_path) {
  if (is.null(config_path)) return(NULL)
  cfg <- read_yaml(config_path)
  if (!is.null(cfg$cluster_labels)) unlist(cfg$cluster_labels) else NULL
}

arm1_labels <- load_labels(arm1_config)
arm2_labels <- load_labels(arm2_config)

cat(sprintf("=== B1 Comparison: %s vs %s ===\n", arm1_name, arm2_name))
cat(sprintf("  Arm 1 dir:    %s\n", arm1_dir))
cat(sprintf("  Arm 2 dir:    %s\n", arm2_dir))
cat(sprintf("  Arm 1 labels: %s\n",
            if (!is.null(arm1_labels)) paste(arm1_labels, collapse = ", ") else "auto"))
cat(sprintf("  Arm 2 labels: %s\n",
            if (!is.null(arm2_labels)) paste(arm2_labels, collapse = ", ") else "auto"))
cat(sprintf("  Output dir:   %s\n", output_dir))

# --- Load results from both arms ---
load_arm <- function(arm_dir, arm_label, cluster_labels) {
  cluster_res <- fread(file.path(arm_dir, "association_results_all.csv"))
  cluster_res <- cluster_res[model_type == "individual"]
  cluster_res[, predictor := gsub("PRS_", "", predictor)]

  gw_res <- fread(file.path(arm_dir, "genome_wide_prs_results.csv"))

  combined <- rbindlist(list(
    cluster_res[, .(group, outcome, predictor, OR, CI_lower, CI_upper, p_value)],
    gw_res[, .(group, outcome, predictor, OR, CI_lower, CI_upper, p_value)]
  ))
  combined[, arm := arm_label]

  # Detect K clusters and apply semantic labels
  k_preds <- sort(unique(grep("^K\\d+$", combined$predictor, value = TRUE)))
  cat(sprintf("  %s: %d cluster results, %d GW results, K=%d\n",
              arm_label, nrow(cluster_res), nrow(gw_res), length(k_preds)))

  # Use provided labels or fall back to generic Kn names
  if (!is.null(cluster_labels)) {
    for (k in names(cluster_labels)) {
      combined[predictor == k, predictor := cluster_labels[[k]]]
    }
  }

  combined
}

arm1_data <- load_arm(arm1_dir, arm1_name, arm1_labels)
arm2_data <- load_arm(arm2_dir, arm2_name, arm2_labels)
plot_data <- rbindlist(list(arm1_data, arm2_data))

# --- Predictor ordering ---
# Arm2 (full) labels first, then any arm1-only labels
arm2_label_vals <- if (!is.null(arm2_labels)) unname(arm2_labels) else character(0)
arm1_label_vals <- if (!is.null(arm1_labels)) unname(arm1_labels) else character(0)
all_labels_ordered <- c(arm2_label_vals,
                        setdiff(arm1_label_vals, arm2_label_vals))
present_labels <- intersect(all_labels_ordered,
                            unique(as.character(plot_data$predictor)))

# Remove genome-wide PRS from comparison plots (cluster-only)
plot_data <- plot_data[!predictor %in% c("GW T2D", "GW CAD")]
plot_data[, predictor := factor(predictor, levels = present_labels)]
plot_data[, arm := factor(arm, levels = c(arm1_name, arm2_name))]

# Outcome labels
outcome_labels <- c("T2D" = "Type 2 Diabetes", "CAD" = "Coronary Artery Disease",
                     "Synthetic" = "T2D or CAD")
plot_data[, outcome_label := outcome_labels[outcome]]
plot_data[, outcome_label := factor(outcome_label, levels = outcome_labels)]

# Ancestry labels
ancestry_labels <- c("afr_validation" = "AFR", "eas_validation" = "EAS",
                      "sas_validation" = "SAS", "eur_validation" = "EUR")
plot_data[, ancestry := ancestry_labels[group]]

# Significance
plot_data[, sig := fifelse(p_value < 0.001, "***",
                  fifelse(p_value < 0.01, "**",
                  fifelse(p_value < 0.05, "*", "")))]

# ============================================================
# EUR validation: side-by-side comparison
# ============================================================
cat("\n--- Creating EUR comparison forest plot ---\n")

outcome_colors <- c("Type 2 Diabetes" = "#E74C3C",
                     "Coronary Artery Disease" = "#3498DB",
                     "T2D or CAD" = "#9B59B6")

arm_shapes <- setNames(c(16, 17), c(arm1_name, arm2_name))

eur_data <- plot_data[group == "eur_validation"]
eur_data <- eur_data[!is.na(predictor)]

p_eur <- ggplot(eur_data, aes(x = predictor, y = OR, color = outcome_label,
                               shape = arm)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.7), width = 0.25, linewidth = 0.6) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  scale_shape_manual(values = arm_shapes, name = "Arm") +
  labs(
    x = "PRS Predictor",
    y = "Odds Ratio per SD (95% CI)",
    title = sprintf("EUR Validation: %s vs %s Cluster PRS", arm1_name, arm2_name)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

ggsave(file.path(output_dir, "comparison_forest_eur.png"),
       p_eur, width = 16, height = 8, dpi = 300, bg = "white")
cat("  Saved: comparison_forest_eur.png\n")

# ============================================================
# Non-EUR: side-by-side comparison faceted by ancestry
# ============================================================
cat("\n--- Creating non-EUR comparison forest plot ---\n")

noneur_data <- plot_data[group %in% c("afr_validation", "eas_validation", "sas_validation")]
noneur_data <- noneur_data[!is.na(predictor)]
noneur_data[, ancestry := factor(ancestry, levels = c("AFR", "EAS", "SAS"))]

p_noneur <- ggplot(noneur_data, aes(x = predictor, y = OR, color = outcome_label,
                                     shape = arm)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.7), size = 2.5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.7), width = 0.25, linewidth = 0.5) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  scale_shape_manual(values = arm_shapes, name = "Arm") +
  facet_wrap(~ancestry, ncol = 1, scales = "free_y") +
  labs(
    x = "PRS Predictor",
    y = "Odds Ratio per SD (95% CI)",
    title = sprintf("Cross-Ancestry: %s vs %s Cluster PRS", arm1_name, arm2_name)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

ggsave(file.path(output_dir, "comparison_forest_noneur.png"),
       p_noneur, width = 16, height = 14, dpi = 300, bg = "white")
cat("  Saved: comparison_forest_noneur.png\n")

# ============================================================
# Summary table: best cluster per outcome per arm
# ============================================================
cat("\n--- Summary: Best cluster PRS per outcome per arm ---\n")

cluster_only <- plot_data[!as.character(predictor) %in% c("GW T2D", "GW CAD") &
                          group == "eur_validation"]
best_per <- cluster_only[, .SD[which.min(p_value)],
                         by = .(arm, outcome)]

cat("\n  Best cluster PRS (EUR validation, lowest p-value):\n")
print(best_per[, .(arm, outcome, predictor, OR, CI_lower, CI_upper, p_value)])

# Save combined data for downstream use
fwrite(plot_data, file.path(output_dir, "comparison_plot_data.csv"))
cat(sprintf("\n  Plot data saved: %s\n", file.path(output_dir, "comparison_plot_data.csv")))

cat("\n=== Done ===\n")
