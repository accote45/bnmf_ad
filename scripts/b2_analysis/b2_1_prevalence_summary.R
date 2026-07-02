#!/usr/bin/env Rscript
# b2_1_prevalence_summary.R
# Comorbidity prevalence and person-years of risk per cluster assignment.
# EUR validation sample only.
#
# Usage:
#   Rscript scripts/b2_analysis/b2_1_prevalence_summary.R --config config/b2_1_config.yaml

library(data.table)
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--config") { config_path <- args[i + 1]; i <- i + 2 }
  else { i <- i + 1 }
}
if (is.null(config_path)) stop("--config is required.")

cfg <- read_yaml(config_path)
results_dir  <- cfg$results_dir
cluster_label_map <- unlist(cfg$cluster_labels)
clusters     <- cfg$analysis$clusters
thresholds   <- cfg$analysis$thresholds
thresh_labels <- cfg$analysis$threshold_labels
sample_group <- cfg$analysis$sample_group
prs_cols     <- paste0("PRS_", clusters)
all_levels   <- c(clusters, "None")

cat("=== B2.1 Prevalence & Person-Years Summary ===\n")
cat(sprintf("  Sample group: %s\n\n", sample_group))

# --- Load and merge ---
prs  <- fread(cfg$b1_results$prs_file)
surv <- fread(cfg$survival_phenotypes)
prs  <- prs[group == sample_group]

dt <- merge(prs[, c("FID", "IID", prs_cols), with = FALSE],
            surv[, .(FID, IID, first_dx, age_at_first_dx, time_years, event)],
            by = c("FID", "IID"))

cat(sprintf("  Individuals with T2D or CAD: %d\n", nrow(dt)))
cat(sprintf("  Comorbidity events: %d (%.1f%%)\n\n", sum(dt$event), 100 * mean(dt$event)))

# --- Percentile ranks ---
pctile_cols <- paste0("pctile_", clusters)
for (k in seq_along(prs_cols)) {
  dt[, (pctile_cols[k]) := frank(get(prs_cols[k])) / .N]
}

# --- Compute per threshold ---
results_list <- list()

for (t_idx in seq_along(thresholds)) {
  thresh <- thresholds[t_idx]
  t_label <- thresh_labels[t_idx]

  dt[, assigned_cluster := {
    pctiles <- numeric(length(pctile_cols))
    for (j in seq_along(pctile_cols)) pctiles[j] <- get(pctile_cols[j])
    eligible <- which(pctiles > thresh)
    if (length(eligible) == 0) "None"
    else clusters[eligible[which.max(pctiles[eligible])]]
  }, by = seq_len(nrow(dt))]

  dt[, assigned_cluster := factor(assigned_cluster, levels = all_levels)]

  for (lev in all_levels) {
    sub <- dt[assigned_cluster == lev]
    results_list[[length(results_list) + 1]] <- data.table(
      threshold = t_label,
      cluster = lev,
      cluster_label = if (lev %in% names(cluster_label_map)) cluster_label_map[[lev]] else lev,
      n_total = nrow(sub),
      n_events = sum(sub$event),
      prevalence_pct = round(100 * mean(sub$event), 2),
      person_years = round(sum(sub$time_years), 1),
      mean_followup_years = round(mean(sub$time_years), 2),
      median_followup_years = round(median(sub$time_years), 2),
      incidence_rate_per_100py = round(100 * sum(sub$event) / sum(sub$time_years), 2)
    )
  }
}

out <- rbindlist(results_list)
out_file <- file.path(results_dir, "prevalence_person_years.csv")
fwrite(out, out_file)
cat(sprintf("Saved: %s (%d rows)\n\n", out_file, nrow(out)))

# Print table
cat(sprintf("%-8s %-6s %6s %6s %8s %12s %10s %12s\n",
            "Thresh", "Clust", "n", "events", "prev%", "person-yrs", "mean_FU", "IR/100py"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (r in seq_len(nrow(out))) {
  row <- out[r]
  cat(sprintf("%-8s %-6s %6d %6d %8.1f %12.1f %10.2f %12.2f\n",
              row$threshold, row$cluster, row$n_total, row$n_events,
              row$prevalence_pct, row$person_years, row$mean_followup_years,
              row$incidence_rate_per_100py))
}
