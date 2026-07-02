#!/usr/bin/env Rscript
# 01_assign_clusters.R
# Assign individuals to top-N% cluster PRS; tabulate comorbidity prevalence.
#
# An individual is assigned to a cluster if they are in the top N% of that
# cluster's PRS distribution (threshold from config). If they qualify for
# multiple clusters, they are assigned to the one where their percentile
# rank is highest.
#
# Usage:
#   Rscript scripts/b2_analysis/01_assign_clusters.R --config config/b2_config.yaml
#   Rscript scripts/b2_analysis/01_assign_clusters.R --config config/b2_1_config.yaml

library(data.table)
library(yaml)

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
if (is.null(config_path)) stop("--config is required. Usage: Rscript 01_assign_clusters.R --config <config.yaml>")

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
clusters <- cfg$analysis$clusters
pctile_thresh <- cfg$analysis$top_prs_percentile
prs_cols <- paste0("PRS_", clusters)

cat("=== B2 Step 1: Cluster Assignment & Prevalence ===\n")
cat(sprintf("  PRS type: %s\n", ifelse(is.null(cfg$analysis$prs_type), "weighted", cfg$analysis$prs_type)))
cat(sprintf("  Config: %s\n", config_path))
cat(sprintf("  Top percentile threshold: %.0f%%\n", 100 * (1 - pctile_thresh)))
cat(sprintf("  Clusters: %s\n\n", paste(clusters, collapse = ", ")))

# --- Load data ---
cat("--- Loading data ---\n")

# PRS scores
prs <- fread(cfg$b1_results$prs_file)
cat(sprintf("  PRS file: %d individuals\n", nrow(prs)))

# Survival phenotypes
surv <- fread(file.path(results_dir, "survival_phenotypes.tsv"))
cat(sprintf("  Survival phenotypes: %d individuals\n", nrow(surv)))

# Sample splits
eur_train_ids <- fread(file.path(cfg$b1_results$samples_dir, "eur_train.keep"))$V1
eur_val_ids <- fread(file.path(cfg$b1_results$samples_dir, "eur_validation.keep"))$V1
cat(sprintf("  EUR train IDs: %d\n", length(eur_train_ids)))
cat(sprintf("  EUR validation IDs: %d\n", length(eur_val_ids)))

# --- Merge PRS with survival data ---
dt <- merge(surv, prs[, c("IID", prs_cols), with = FALSE], by = "IID")
cat(sprintf("\n  Merged (surv + PRS): %d individuals\n", nrow(dt)))

# --- Process each split ---
process_split <- function(dt, split_ids, split_name) {
  cat(sprintf("\n--- Processing %s ---\n", split_name))

  sub <- dt[IID %in% split_ids]
  cat(sprintf("  Individuals with T2D/CAD + PRS: %d\n", nrow(sub)))

  if (nrow(sub) == 0) {
    cat("  WARNING: No individuals in this split with survival data!\n")
    return(NULL)
  }

  # Compute percentile ranks for each cluster PRS
  pctile_cols <- paste0("pctile_", clusters)
  for (k in seq_along(prs_cols)) {
    sub[, (pctile_cols[k]) := frank(get(prs_cols[k])) / .N]
  }

  # Assign to cluster: top N% (pctile > threshold), pick highest rank
  sub[, assigned_cluster := {
    pctiles <- c()
    for (pc in pctile_cols) pctiles <- c(pctiles, get(pc))
    eligible <- which(pctiles > pctile_thresh)
    if (length(eligible) == 0) {
      NA_character_
    } else {
      clusters[eligible[which.max(pctiles[eligible])]]
    }
  }, by = seq_len(nrow(sub))]

  # Report assignment counts
  n_assigned <- sum(!is.na(sub$assigned_cluster))
  n_unassigned <- sum(is.na(sub$assigned_cluster))
  cat(sprintf("  Assigned: %d (%.1f%%)\n", n_assigned,
              100 * n_assigned / nrow(sub)))
  cat(sprintf("  Unassigned (not in top %.0f%% of any cluster): %d (%.1f%%)\n",
              100 * (1 - pctile_thresh),
              n_unassigned, 100 * n_unassigned / nrow(sub)))

  cat("\n  Assignment by cluster:\n")
  for (k in clusters) {
    nk <- sum(sub$assigned_cluster == k, na.rm = TRUE)
    cat(sprintf("    %s: %d\n", k, nk))
  }

  # --- Comorbidity prevalence per cluster ---
  cat("\n  Comorbidity prevalence by cluster:\n")

  prev_rows <- list()
  for (k in clusters) {
    ksub <- sub[assigned_cluster == k]
    n_total <- nrow(ksub)
    n_cad_first <- sum(ksub$first_dx == "CAD")
    n_t2d_first <- sum(ksub$first_dx == "T2D")
    n_events <- sum(ksub$event == 1)
    n_events_cad_to_t2d <- sum(ksub$first_dx == "CAD" & ksub$event == 1)
    n_events_t2d_to_cad <- sum(ksub$first_dx == "T2D" & ksub$event == 1)
    prev <- if (n_total > 0) round(100 * n_events / n_total, 1) else NA
    prev_cad_t2d <- if (n_cad_first > 0) round(100 * n_events_cad_to_t2d / n_cad_first, 1) else NA
    prev_t2d_cad <- if (n_t2d_first > 0) round(100 * n_events_t2d_to_cad / n_t2d_first, 1) else NA

    cat(sprintf("    %s: n=%d, events=%d (%.1f%%), CAD→T2D=%d/%-5d (%.1f%%), T2D→CAD=%d/%-5d (%.1f%%)\n",
                k, n_total, n_events, prev,
                n_events_cad_to_t2d, n_cad_first, ifelse(is.na(prev_cad_t2d), 0, prev_cad_t2d),
                n_events_t2d_to_cad, n_t2d_first, ifelse(is.na(prev_t2d_cad), 0, prev_t2d_cad)))

    prev_rows[[k]] <- data.table(
      cluster = k,
      n_total = n_total,
      n_cad_first = n_cad_first,
      n_t2d_first = n_t2d_first,
      n_events_total = n_events,
      n_events_cad_to_t2d = n_events_cad_to_t2d,
      n_events_t2d_to_cad = n_events_t2d_to_cad,
      prevalence_overall_pct = prev,
      prevalence_cad_to_t2d_pct = prev_cad_t2d,
      prevalence_t2d_to_cad_pct = prev_t2d_cad
    )
  }
  prev_dt <- rbindlist(prev_rows)

  # Also add unassigned row
  unsub <- sub[is.na(assigned_cluster)]
  if (nrow(unsub) > 0) {
    prev_dt <- rbind(prev_dt, data.table(
      cluster = "Unassigned",
      n_total = nrow(unsub),
      n_cad_first = sum(unsub$first_dx == "CAD"),
      n_t2d_first = sum(unsub$first_dx == "T2D"),
      n_events_total = sum(unsub$event == 1),
      n_events_cad_to_t2d = sum(unsub$first_dx == "CAD" & unsub$event == 1),
      n_events_t2d_to_cad = sum(unsub$first_dx == "T2D" & unsub$event == 1),
      prevalence_overall_pct = round(100 * mean(unsub$event), 1),
      prevalence_cad_to_t2d_pct = round(100 * mean(unsub$event[unsub$first_dx == "CAD"]), 1),
      prevalence_t2d_to_cad_pct = round(100 * mean(unsub$event[unsub$first_dx == "T2D"]), 1)
    ))
  }

  # Write outputs
  assign_file <- file.path(results_dir, sprintf("cluster_assignments_%s.tsv", split_name))
  out_cols <- c("FID", "IID", "first_dx", "event", "time_years", "assigned_cluster",
                prs_cols, pctile_cols)
  fwrite(sub[, ..out_cols], assign_file, sep = "\t")
  cat(sprintf("\n  Saved assignments: %s\n", assign_file))

  prev_file <- file.path(results_dir, sprintf("cluster_prevalence_%s.tsv", split_name))
  fwrite(prev_dt, prev_file, sep = "\t")
  cat(sprintf("  Saved prevalence: %s\n", prev_file))

  return(sub)
}

# --- Run for EUR train and validation ---
train_dt <- process_split(dt, eur_train_ids, "train")
val_dt   <- process_split(dt, eur_val_ids, "validation")

cat("\n=== Step 1 complete ===\n")
