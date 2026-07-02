#!/usr/bin/env Rscript
# 03_cumulative_hazard_plots.R
# Generate cumulative hazard curves for each cluster PRS.
# One plot per direction (CAD-first->T2D, T2D-first->CAD, Combined).
# EUR validation set only.
#
# Usage:
#   Rscript scripts/b2_analysis/03_cumulative_hazard_plots.R --config config/b2_config.yaml
#   Rscript scripts/b2_analysis/03_cumulative_hazard_plots.R --config config/b2_1_config.yaml

library(data.table)
library(yaml)
library(survival)
library(ggplot2)

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") {
    config_path <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}
if (is.null(config_path)) stop("--config is required. Usage: Rscript 03_cumulative_hazard_plots.R --config <config.yaml>")

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
clusters <- cfg$analysis$clusters

cat("=== B2 Step 3: Cumulative Hazard Plots ===\n")
cat(sprintf("  Config: %s\n\n", config_path))

# --- Load data ---
cat("--- Loading data ---\n")

# Cluster assignments (validation)
assign_dt <- fread(file.path(results_dir, "cluster_assignments_validation.tsv"))
cat(sprintf("  Validation cluster assignments: %d individuals\n", nrow(assign_dt)))

# Survival phenotypes (for any additional fields)
surv <- fread(file.path(results_dir, "survival_phenotypes.tsv"))

# Merge to get full survival info
dt <- merge(assign_dt, surv[, .(IID, first_dx_date, second_dx_date)], by = "IID", all.x = TRUE)
cat(sprintf("  Merged: %d individuals\n", nrow(dt)))

# Filter to assigned individuals only
dt_assigned <- dt[!is.na(assigned_cluster)]
cat(sprintf("  Assigned individuals: %d\n", nrow(dt_assigned)))

# --- Cluster colors ---
cluster_colors <- c(
  K1 = "#E74C3C", K2 = "#3498DB", K3 = "#2ECC71",
  K4 = "#F39C12", K5 = "#9B59B6", K6 = "#1ABC9C"
)

# --- Function to extract cumulative hazard + CI from survfit ---
extract_cumhaz <- function(fit) {
  # Extract data from survfit object including confidence intervals
  # For Nelson-Aalen, cumhaz = -log(surv), CI from surv CI
  strata_names <- names(fit$strata)
  n_strata <- length(fit$strata)

  plot_list <- list()
  idx <- 0
  for (s in seq_len(n_strata)) {
    n_times <- fit$strata[s]
    rows <- (idx + 1):(idx + n_times)
    cluster_label <- gsub("assigned_cluster=", "", strata_names[s])

    # Cumulative hazard CI from survival CI: H = -log(S)
    # lower CI of H = -log(upper CI of S), upper CI of H = -log(lower CI of S)
    cumhaz_lower <- -log(fit$upper[rows])
    cumhaz_upper <- -log(fit$lower[rows])

    plot_list[[s]] <- data.table(
      time = fit$time[rows],
      cumhaz = fit$cumhaz[rows],
      cumhaz_lower = cumhaz_lower,
      cumhaz_upper = cumhaz_upper,
      cluster = cluster_label
    )
    idx <- idx + n_times
  }
  rbindlist(plot_list)
}

# --- Generate plots for each direction ---
directions <- list(
  list(dx = "CAD", label = "CAD-first: Cumulative hazard of developing T2D",
       file = "cumhaz_cad_to_t2d.png"),
  list(dx = "T2D", label = "T2D-first: Cumulative hazard of developing CAD",
       file = "cumhaz_t2d_to_cad.png"),
  list(dx = "Combined", label = "Combined: Cumulative hazard of developing comorbidity",
       file = "cumhaz_combined.png")
)

for (dir_info in directions) {
  cat(sprintf("\n--- %s ---\n", dir_info$label))

  if (dir_info$dx == "Combined") {
    sub <- copy(dt_assigned)
  } else {
    sub <- dt_assigned[first_dx == dir_info$dx]
  }
  cat(sprintf("  Individuals: %d, Events: %d\n", nrow(sub), sum(sub$event)))

  # Check per-cluster event counts
  for (k in clusters) {
    nk <- sum(sub$assigned_cluster == k)
    ek <- sum(sub$assigned_cluster == k & sub$event == 1)
    cat(sprintf("    %s: n=%d, events=%d\n", k, nk, ek))
  }

  if (nrow(sub) == 0 || sum(sub$event) == 0) {
    cat("  SKIPPING: No events in this direction\n")
    next
  }

  # Ensure cluster is a factor with correct levels
  sub[, assigned_cluster := factor(assigned_cluster, levels = clusters)]

  # Remove clusters with 0 individuals
  present_clusters <- intersect(clusters, unique(as.character(sub$assigned_cluster)))

  # Fit Nelson-Aalen estimator
  fit <- survfit(Surv(time_years, event) ~ assigned_cluster,
                 data = sub, type = "fleming-harrington")

  # Extract cumulative hazard data
  plot_dt <- extract_cumhaz(fit)
  plot_dt[, cluster := factor(cluster, levels = clusters)]

  # Add time=0, cumhaz=0 for each cluster (origin)
  origins <- data.table(
    time = 0,
    cumhaz = 0,
    cumhaz_lower = 0,
    cumhaz_upper = 0,
    cluster = factor(present_clusters, levels = clusters)
  )
  plot_dt <- rbind(origins, plot_dt)

  # Build fill colors (lighter versions for CI ribbon)
  cluster_fills <- adjustcolor(cluster_colors, alpha.f = 0.15)
  names(cluster_fills) <- names(cluster_colors)

  # Plot with confidence intervals
  p <- ggplot(plot_dt, aes(x = time, color = cluster, fill = cluster)) +
    geom_ribbon(aes(ymin = cumhaz_lower, ymax = cumhaz_upper),
                alpha = 0.15, linewidth = 0, show.legend = FALSE) +
    geom_step(aes(y = cumhaz), linewidth = 1) +
    scale_color_manual(values = cluster_colors, name = "Cluster",
                       drop = FALSE) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    labs(
      x = "Years from first diagnosis",
      y = "Cumulative hazard",
      title = dir_info$label
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(0, 20), ylim = c(0, 0.5))

  out_file <- file.path(results_dir, dir_info$file)
  ggsave(out_file, p, width = 10, height = 6, dpi = 300)
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("\n=== Step 3 complete ===\n")
