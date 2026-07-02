#!/usr/bin/env Rscript
# 06_compare_pathway_vs_cluster.R
# Select top 9 pathways by Synthetic OR/SD (validation) and create comparison
# forest plot against bNMF cluster PRS from b1.
#
# Usage:
#   Rscript scripts/b1_2_analysis/06_compare_pathway_vs_cluster.R
#   Rscript scripts/b1_2_analysis/06_compare_pathway_vs_cluster.R --config config/b1_2_config.yaml

library(data.table)
library(yaml)
library(ggplot2)

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
source(file.path(base_dir, "scripts/a1_analysis/figure_utils.R"))

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_2_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}

cfg <- read_yaml(config_path)
results_dir    <- cfg$results_dir
n_top          <- cfg$analysis$n_top_pathways
cluster_labels <- unlist(cfg$cluster_labels)

cat("=== B1.2 Step 6: Compare Pathway vs Cluster PRS ===\n")

# --- Load pathway association results ---
cat("\n--- Loading pathway results ---\n")
pw_results <- fread(file.path(results_dir, "pathway_association_results.csv"))
cat(sprintf("  Pathway results: %d rows (%d pathways x %d outcomes)\n",
            nrow(pw_results),
            length(unique(pw_results$pathway)),
            length(unique(pw_results$outcome))))

# --- Select top pathways by Synthetic validation OR ---
cat(sprintf("\n--- Selecting top %d pathways by Synthetic OR (validation) ---\n", n_top))
synth_pw <- pw_results[outcome == "Synthetic"][order(-OR_val)]
top_pathways <- synth_pw$pathway[1:min(n_top, nrow(synth_pw))]

cat(sprintf("  Top %d pathways:\n", n_top))
for (i in seq_along(top_pathways)) {
  row <- synth_pw[pathway == top_pathways[i]]
  cat(sprintf("    %d. %s: OR=%.3f [%.3f-%.3f], P=%.2e\n",
              i, top_pathways[i], row$OR_val,
              row$CI_lower_val, row$CI_upper_val, row$p_value_val))
}

# --- Get pathway results for all outcomes ---
top_all <- pw_results[pathway %in% top_pathways]
top_all[, short_name := gsub("_", " ", pathway)]
top_all[, short_name := fifelse(nchar(short_name) > 40,
                                paste0(substr(short_name, 1, 37), "..."),
                                short_name)]

# --- Save detailed top pathway results ---
top_detail <- file.path(results_dir, sprintf("top%d_pathway_details.csv", n_top))
fwrite(top_all, top_detail)
cat(sprintf("\n  Saved: %s\n", top_detail))

# --- Load b1 cluster PRS results ---
cat("\n--- Loading b1 cluster PRS results ---\n")
b1_results <- fread(cfg$b1_results)
cluster_val <- b1_results[group == "eur_validation" & model_type == "individual"]
cat(sprintf("  b1 EUR validation individual results: %d rows\n", nrow(cluster_val)))

# --- Prepare plot data ---
cat("\n--- Preparing comparison plot ---\n")

# Cluster PRS data — use readable labels from config
cluster_plot <- cluster_val[, .(
  predictor = cluster_labels[gsub("PRS_", "", predictor)],
  outcome = outcome,
  OR = OR,
  CI_lower = CI_lower,
  CI_upper = CI_upper,
  p_value = p_value,
  type = "bNMF Cluster"
)]

# Pathway PRS data — label as P1..P9
pathway_label_map <- data.table(
  pathway = top_pathways,
  label = paste0("P", seq_along(top_pathways))
)

pw_plot <- top_all[, .(
  predictor = pathway,
  outcome = outcome,
  OR = OR_val,
  CI_lower = CI_lower_val,
  CI_upper = CI_upper_val,
  p_value = p_value_val,
  type = "Pathway (PRSet)"
)]
pw_plot <- merge(pw_plot, pathway_label_map, by.x = "predictor", by.y = "pathway")
pw_plot[, predictor := label]
pw_plot[, label := NULL]

# Combine
plot_data <- rbindlist(list(cluster_plot, pw_plot), fill = TRUE)

# Factor ordering: cluster labels then P1-P9
cluster_level_order <- cluster_labels  # named vector K1-K9
pathway_levels <- paste0("P", seq_len(n_top))
pred_levels <- c(cluster_level_order, pathway_levels)
plot_data[, predictor := factor(predictor, levels = pred_levels)]

# Outcome labels and colors
outcome_labels <- c("T2D" = "Type 2 Diabetes", "CAD" = "Coronary Artery Disease",
                     "Synthetic" = "T2D or CAD")
plot_data[, outcome_label := outcome_labels[outcome]]
plot_data[, outcome_label := factor(outcome_label, levels = outcome_labels)]

outcome_colors <- c("Type 2 Diabetes" = "#E74C3C",
                     "Coronary Artery Disease" = "#3498DB",
                     "T2D or CAD" = "#9B59B6")

n_clusters <- length(cluster_labels)
n_total <- n_clusters + n_top

# --- Forest plot ---
p <- ggplot(plot_data, aes(x = predictor, y = OR, color = outcome_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.6), width = 0.25, linewidth = 0.7) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  # Shade pathway region
  annotate("rect", xmin = n_clusters + 0.5, xmax = n_total + 0.5,
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey50") +
  # Section labels
  annotate("text", x = (n_clusters + 1) / 2, y = Inf,
           label = "bNMF Clusters", vjust = 1.5, fontface = "bold", size = 4) +
  annotate("text", x = n_clusters + (n_top + 1) / 2, y = Inf,
           label = sprintf("Top %d Pathways", n_top), vjust = 1.5,
           fontface = "bold", size = 4) +
  labs(
    x = "PRS Predictor",
    y = "Odds Ratio per SD (95% CI)",
    title = "EUR Validation: bNMF Cluster PRS vs Top Pathway PRS"
  ) +
  theme_big_text(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

out_plot <- file.path(results_dir, "forest_cluster_vs_pathway.png")
ggsave(out_plot, p, width = 18, height = 7, dpi = 600)
cat(sprintf("  Saved: %s\n", out_plot))

# --- Summary table ---
cat("\n--- Creating summary table ---\n")

# Build pathway legend (P1 = full name)
pw_legend <- data.table(
  label = paste0("P", seq_along(top_pathways)),
  full_pathway_name = top_pathways
)

summary_dt <- dcast(plot_data, predictor + type ~ outcome,
                    value.var = c("OR", "CI_lower", "CI_upper", "p_value"))
summary_dt <- merge(summary_dt, pw_legend, by.x = "predictor", by.y = "label", all.x = TRUE)

out_summary <- file.path(results_dir, "comparison_summary.csv")
fwrite(summary_dt, out_summary)
cat(sprintf("  Saved: %s\n", out_summary))

# --- Print comparison ---
cat("\n=== Comparison Summary ===\n")
cat("\nbNMF Cluster PRS (EUR validation):\n")
print(cluster_plot[, .(predictor, outcome, OR, CI_lower, CI_upper, p_value)])

cat(sprintf("\nTop %d Pathway PRS (EUR validation):\n", n_top))
for (i in seq_along(top_pathways)) {
  cat(sprintf("  P%d = %s\n", i, top_pathways[i]))
}
print(pw_plot[, .(predictor, outcome, OR, CI_lower, CI_upper, p_value)])

cat("\n=== Step 6 complete ===\n")
