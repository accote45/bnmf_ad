#!/usr/bin/env Rscript
# 07_decile_event_rate.R
# Decile-stratified event rate plots for each bNMF cluster PRS.
# For each cluster, plots T2D / CAD / T2D+CAD event rates across PGS deciles.
#
# Usage:
#   Rscript scripts/b1_analysis/07_decile_event_rate.R --config config/b1_config.yaml --hard-assignment

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
}

cat("=== B1 Step 7: Decile Event Rate Plots ===\n")
cat(sprintf("  Output dir: %s\n", results_dir))

# ============================================================
# Load and merge data
# ============================================================
cat("\n--- Loading data ---\n")

prs_all <- fread(file.path(results_dir, "prs/cluster_prs_all.tsv"))
pheno   <- fread(cfg$phenotypes$phenotype_file)

prs_cols <- grep("^PRS_K", colnames(prs_all), value = TRUE)
n_clusters <- length(prs_cols)
cat(sprintf("  PRS: %d individuals, %d clusters\n", nrow(prs_all), n_clusters))

dat <- merge(prs_all, pheno, by = c("FID", "IID"), all.x = TRUE)
dat[, Synthetic := fifelse(T2D == 1 & CAD == 1, 1L,
                   fifelse(T2D == 0 | CAD == 0, 0L, NA_integer_))]

# Filter to EUR validation
eur_val <- dat[group == "eur_validation"]
cat(sprintf("  EUR validation: %d individuals\n", nrow(eur_val)))
cat(sprintf("  T2D cases: %d, CAD cases: %d, T2D+CAD cases: %d\n",
            sum(eur_val$T2D == 1, na.rm = TRUE),
            sum(eur_val$CAD == 1, na.rm = TRUE),
            sum(eur_val$Synthetic == 1, na.rm = TRUE)))

# ============================================================
# Compute decile event rates
# ============================================================
cat("\n--- Computing decile event rates ---\n")

# data.table ntile equivalent
dt_ntile <- function(x, n = 10L) {
  ranks <- frank(x, ties.method = "random", na.last = "keep")
  as.integer(ceiling(ranks / (sum(!is.na(x)) / n)))
}

outcomes <- c("T2D", "CAD", "Synthetic")
outcome_labels <- c("T2D" = "Type 2 Diabetes",
                     "CAD" = "Coronary Artery Disease",
                     "Synthetic" = "T2D and CAD")

# Cluster labels from config
if (!is.null(cfg$cluster_labels)) {
  cluster_labels <- unlist(cfg$cluster_labels)
} else {
  cluster_labels <- setNames(gsub("PRS_", "", prs_cols), gsub("PRS_", "", prs_cols))
}

results_list <- list()

for (prs_col in prs_cols) {
  k_id <- gsub("PRS_", "", prs_col)
  cluster_name <- cluster_labels[[k_id]]

  eur_val[, decile := dt_ntile(get(prs_col), 10L)]

  for (outcome in outcomes) {
    agg <- eur_val[!is.na(get(outcome)) & !is.na(decile),
                   .(event_rate = mean(get(outcome))),
                   by = decile]
    agg[, outcome_label := outcome_labels[[outcome]]]
    agg[, cluster := cluster_name]
    agg[, cluster_k := k_id]
    results_list[[length(results_list) + 1]] <- agg
  }

  eur_val[, decile := NULL]
}

plot_dt <- rbindlist(results_list)

# Factor ordering
k_ids <- gsub("PRS_", "", prs_cols)
cluster_order <- unname(cluster_labels[k_ids])
plot_dt[, cluster := factor(cluster, levels = cluster_order)]
plot_dt[, outcome_label := factor(outcome_label, levels = outcome_labels)]

cat(sprintf("  Computed event rates: %d rows\n", nrow(plot_dt)))

# ============================================================
# Plot: 5x2 grid of decile event rate scatterplots
# ============================================================
cat("\n--- Creating decile event rate figure ---\n")

outcome_colors <- c("Type 2 Diabetes" = "#E74C3C",
                     "Coronary Artery Disease" = "#3498DB",
                     "T2D and CAD" = "#9B59B6")

p <- ggplot(plot_dt, aes(x = decile, y = event_rate,
                          color = outcome_label, group = outcome_label)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~cluster, ncol = 2, nrow = 5, scales = "free_y") +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "PGS Decile", y = "Event Rate") +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

out_path <- file.path(results_dir, "decile_event_rate.png")
ggsave(out_path, p, width = 10, height = 16, dpi = 300)
cat(sprintf("  Saved: %s\n", out_path))

cat("\n=== Done ===\n")
