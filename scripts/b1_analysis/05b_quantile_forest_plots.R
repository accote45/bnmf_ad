#!/usr/bin/env Rscript
# 05b_quantile_forest_plots.R
# Top-10% forest plots (EUR + cross-ancestry) for the quantile cluster-PGS
# analysis. Reads the association table written by 05_quantile_prs_analysis.R —
# no regressions are run here, so figure tweaks (labels, colors, sizes) don't
# require rerunning the stats.
#
# Both the cluster PGS and the genome-wide CAD/T2D PGS are plotted as the SAME
# top-decile model (top 10% vs remaining 90%) so the comparison is like-for-like.
#
# Inputs (under <results_dir>):
#   quantile_prs_associations.csv   (from 05_quantile_prs_analysis.R; includes
#                                    GW T2D / GW CAD top-decile rows)
#
# Usage: Rscript scripts/b1_analysis/05b_quantile_forest_plots.R [--config ...] [--hard-assignment]

library(data.table)
library(ggplot2)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}
hard_assignment <- "--hard-assignment" %in% args

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
if (hard_assignment) {
  results_dir <- file.path(results_dir, "prs_hard_assignment")
  cat("  ** HARD ASSIGNMENT MODE **\n")
  cat(sprintf("  Output dir: %s\n", results_dir))
}

cat("=== B1 Step 5b: Quantile PRS Forest Plots ===\n")

# ============================================================
# Load association results (no regressions here)
# ============================================================
assoc_path <- file.path(results_dir, "quantile_prs_associations.csv")
if (!file.exists(assoc_path)) {
  stop(sprintf("%s not found — run 05_quantile_prs_analysis.R first.", assoc_path))
}
results <- fread(assoc_path)
cat(sprintf("  Loaded %d association rows from %s\n", nrow(results), basename(assoc_path)))

# Number of cluster predictors (K1..Kn), inferred from the association table
n_clusters <- length(grep("^K[0-9]+$", unique(results$cluster), value = TRUE))
cat(sprintf("  Clusters: %d\n", n_clusters))

# Cluster labels (config-driven; fall back to raw K labels)
cluster_labels <- if (!is.null(cfg$cluster_labels)) unlist(cfg$cluster_labels) else setNames(paste0("K", 1:n_clusters), paste0("K", 1:n_clusters))

# ============================================================
# Assemble plot data: cluster + genome-wide PGS, all top-decile
# ============================================================
outcome_colors <- c("Type 2 Diabetes" = "#E74C3C",
                     "Coronary Artery Disease" = "#3498DB",
                     "T2D and CAD" = "#9B59B6")
outcome_labels <- c("T2D" = "Type 2 Diabetes", "CAD" = "Coronary Artery Disease",
                     "Synthetic" = "T2D and CAD")
ancestry_labels <- c("afr_validation" = "AFR", "eas_validation" = "EAS",
                      "sas_validation" = "SAS", "eur_validation" = "EUR")

# Top-decile (top 10% vs remaining 90%) ORs for every predictor. The GW CAD/T2D
# PGS rows (cluster_label "GW T2D"/"GW CAD") are already top-decile models in the
# quantile association table, so we plot those rather than the per-SD GW ORs —
# matching the cluster PGS so the contrast is like-for-like.
top10 <- results[quantile == "top_10pct"]
top10[, predictor := cluster_label]

plot_data <- top10[, .(group, outcome, predictor, OR, CI_lower, CI_upper, p_value)]

plot_data[, outcome_label := factor(outcome_labels[outcome], levels = outcome_labels)]
plot_data[, ancestry := ancestry_labels[group]]

# Predictor ordering: GW T2D (leftmost), clusters in desired order, GW CAD (rightmost)
k_preds <- paste0("K", 1:n_clusters)
# Desired META cluster display order: Glycemic, Obesity, SHBG, Adiponectin,
# Triglycerides-HDL, ALP-LDL, Metabolic, Platelet, Blood Pressure-Stature, Lpa
desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
k_ordered <- c(desired_k_order[desired_k_order %in% k_preds],
               sort(setdiff(k_preds, desired_k_order)))
pred_levels <- c("GW T2D", unname(cluster_labels[k_ordered]), "GW CAD")
plot_data[, predictor := factor(predictor, levels = pred_levels)]
# Drop predictors outside these levels (combined-variant GW PGS, which are only
# shown in the supplementary figures) so they don't render as an "NA" column.
plot_data <- plot_data[!is.na(predictor)]

# Grey shading behind the two genome-wide predictors (now leftmost + rightmost)
gw_band <- data.frame(xmin = c(0.5, n_clusters + 1.5),
                      xmax = c(1.5, n_clusters + 2.5),
                      ymin = -Inf, ymax = Inf)

# ============================================================
# Forest plots: top 10%
# ============================================================
cat("\n--- Creating top 10% forest plots ---\n")

# --- EUR validation forest plot ---
eur_data <- plot_data[group == "eur_validation"]

p_eur <- ggplot(eur_data, aes(x = predictor, y = OR, color = outcome_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.6), width = 0.25, linewidth = 0.7) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  geom_rect(data = gw_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, alpha = 0.08, fill = "grey50") +
  labs(
    x = "PGS Predictor",
    y = "OR in Top Decile",
    title = "Top 10% Cluster PGS vs. Genome-Wide PGS (top 10%)"
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
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "forest_eur_top10pct.png"),
       p_eur, width = 12, height = 7, dpi = 300)
cat("  Saved: forest_eur_top10pct.png\n")

# --- Non-EUR combined forest plot ---
noneur_data <- plot_data[group %in% c("afr_validation", "eas_validation", "sas_validation")]
noneur_data[, ancestry := factor(ancestry, levels = c("AFR", "EAS", "SAS"))]

p_noneur <- ggplot(noneur_data, aes(x = predictor, y = OR, color = outcome_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.6), width = 0.25, linewidth = 0.6) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  geom_rect(data = gw_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, alpha = 0.08, fill = "grey50") +
  facet_wrap(~ancestry, ncol = 1, scales = "free_y") +
  labs(
    x = "PRS Predictor",
    y = "OR in Top Decile",
    title = "Cross-Ancestry: Top 10% Cluster PRS vs Genome-Wide PRS (top 10%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "forest_noneur_top10pct.png"),
       p_noneur, width = 12, height = 12, dpi = 300)
cat("  Saved: forest_noneur_top10pct.png\n")

# ============================================================
# Supplementary: top 10% forests including combined-variant GW PGS
# ============================================================
cat("\n--- Creating supplementary top 10% forests (incl. combined-variant GW PGS) ---\n")

plot_supp <- top10[, .(group, outcome, predictor, OR, CI_lower, CI_upper, p_value)]
plot_supp[, outcome_label := factor(outcome_labels[outcome], levels = outcome_labels)]
plot_supp[, ancestry := ancestry_labels[group]]

# Order: GW T2D pair (leftmost), clusters, GW CAD pair (rightmost) — each
# trait-specific GW PGS immediately followed by its combined-variant version.
pred_levels_supp <- c("GW T2D", "GW T2D (combined variants)",
                      unname(cluster_labels[k_ordered]),
                      "GW CAD", "GW CAD (combined variants)")
plot_supp[, predictor := factor(predictor, levels = pred_levels_supp)]
plot_supp <- plot_supp[!is.na(predictor)]

# Grey bands behind the two GW pairs (leftmost two positions + rightmost two).
gw_band_supp <- data.frame(xmin = c(0.5, n_clusters + 2.5),
                           xmax = c(2.5, n_clusters + 4.5),
                           ymin = -Inf, ymax = Inf)

supp_top10_forest <- function(df, facet, title, fname, w, h) {
  p <- ggplot(df, aes(x = predictor, y = OR, color = outcome_label)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_rect(data = gw_band_supp,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, alpha = 0.08, fill = "grey50") +
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.6), width = 0.25, linewidth = 0.6) +
    scale_color_manual(values = outcome_colors, name = "Outcome") +
    labs(x = "PGS Predictor", y = "OR in Top Decile", title = title) +
    theme_minimal(base_size = 14) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 11),
          legend.position = "bottom", legend.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
  if (facet) {
    p <- p + facet_wrap(~ancestry, ncol = 1, scales = "free_y") +
      theme(strip.text = element_text(face = "bold", size = 14))
  }
  ggsave(file.path(results_dir, fname), p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

supp_top10_forest(plot_supp[group == "eur_validation"], FALSE,
  "Top 10% (supp.): Cluster vs Trait-specific & Combined-variant GW PGS",
  "forest_eur_top10pct_supp.png", 14, 7)

sn <- plot_supp[group %in% c("afr_validation", "eas_validation", "sas_validation")]
sn[, ancestry := factor(ancestry, levels = c("AFR", "EAS", "SAS"))]
supp_top10_forest(sn, TRUE,
  "Cross-Ancestry Top 10% (supp.): Cluster vs Trait-specific & Combined-variant GW PGS",
  "forest_noneur_top10pct_supp.png", 14, 12)

cat("\n=== Done ===\n")
