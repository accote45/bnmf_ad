#!/usr/bin/env Rscript
# b2_1_analysis.R
# Simplified B2.1 survival analysis: Cox models and cumulative hazard curves
# for META cluster PRS at 3 thresholds (top 33%, 20%, 10%).
#
# Outcome: T2D/CAD composite — age of onset = first diagnosis of either,
# event = second comorbidity, time = years between diagnoses.
#
# Cumulative hazard is covariate-adjusted via Cox model (age, sex, PC1-10).
#
# Usage:
#   Rscript scripts/b2_analysis/b2_1_analysis.R --config config/b2_1_config.yaml
#   Rscript scripts/b2_analysis/b2_1_analysis.R --config config/b2_1_config.yaml --hard-assignment

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
if (is.null(config_path)) stop("--config is required.")

cfg <- read_yaml(config_path)

hard_assignment <- "--hard-assignment" %in% args

if (hard_assignment) {
  cfg$b1_results$prs_file <- "results/b1_analysis/prs_hard_assignment/prs/cluster_prs_all.tsv"
  cfg$results_dir <- "results/b2_1_analysis_hard_assignment"
  cat("  ** HARD ASSIGNMENT MODE **\n")
  cat(sprintf("  PRS file: %s\n", cfg$b1_results$prs_file))
  cat(sprintf("  Output dir: %s\n\n", cfg$results_dir))
}

results_dir  <- cfg$results_dir
clusters     <- cfg$analysis$clusters
thresholds   <- cfg$analysis$thresholds
thresh_labels <- cfg$analysis$threshold_labels
sample_group <- cfg$analysis$sample_group
prs_cols     <- paste0("PRS_", clusters)

# Cluster display labels (e.g., K1 -> "SHBG")
cluster_labels <- unlist(cfg$cluster_labels)
# Build display names: "SHBG", "Insulin Sensitivity", ..., "None"
cluster_display <- setNames(
  c(cluster_labels[clusters], "None"),
  c(clusters, "None")
)

# Manual color palette mapped to display labels
all_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#6A3D9A", "#FF7F00",
                "#FB9A99", "#A6CEE3", "#B2DF8A", "#B15928", "#FDBF6F",
                "#CAB2D6")
cluster_colors <- setNames(
  all_colors[seq_len(length(clusters) + 1)],
  c(cluster_labels[clusters], "None")
)
# Genome-wide PGS comparison curves (overlaid on the top33 direction plots).
# Black distinguishes them from the Paired-style cluster palette and "None".
cluster_colors <- c(cluster_colors, "GW CAD" = "#000000", "GW T2D" = "#000000")

cat("=== B2.1 Survival Analysis ===\n")
cat(sprintf("  Clusters: %s\n", paste(clusters, collapse = ", ")))
cat(sprintf("  Thresholds: %s\n", paste(thresh_labels, collapse = ", ")))
cat(sprintf("  Sample group: %s\n\n", sample_group))

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# 1. Load and merge data
# =====================================================================

cat("--- Loading data ---\n")

prs  <- fread(cfg$b1_results$prs_file)
surv <- fread(cfg$survival_phenotypes)
covar <- fread(cfg$phenotypes$covariate_file)

# Genome-wide PGS (standalone CAD/T2D scores) for the comparison overlay.
gw_cad_col <- cfg$gw_prs$cad_col
gw_t2d_col <- cfg$gw_prs$t2d_col
gw_prs <- fread(cfg$gw_prs$file)
gw_prs <- gw_prs[group == sample_group]

# Filter PRS to target sample group
prs <- prs[group == sample_group]
cat(sprintf("  PRS (%s): %d individuals\n", sample_group, nrow(prs)))
cat(sprintf("  GW PRS (%s): %d individuals\n", sample_group, nrow(gw_prs)))
cat(sprintf("  Survival phenotypes: %d individuals\n", nrow(surv)))
cat(sprintf("  Covariates: %d individuals\n", nrow(covar)))

# Inner join: PRS + GW PRS + survival + covariates
dt <- merge(prs[, c("FID", "IID", prs_cols), with = FALSE],
            gw_prs[, c("FID", "IID", gw_cad_col, gw_t2d_col), with = FALSE],
            by = c("FID", "IID"))
dt <- merge(dt,
            surv[, .(FID, IID, first_dx, age_at_first_dx, time_years, event)],
            by = c("FID", "IID"))
dt <- merge(dt,
            covar[, .(IID, sex, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)],
            by = "IID")

cat(sprintf("  Merged (individuals with T2D or CAD diagnosis + PRS + covariates): %d\n", nrow(dt)))
cat(sprintf("  Events (developed comorbidity): %d (%.1f%%)\n",
            sum(dt$event), 100 * mean(dt$event)))

# Precompute age squared
dt[, age_at_first_dx2 := age_at_first_dx^2]

# =====================================================================
# 2. Compute percentile ranks (once, reused across thresholds)
# =====================================================================

cat("\n--- Computing percentile ranks ---\n")
pctile_cols <- paste0("pctile_", clusters)
for (k in seq_along(prs_cols)) {
  dt[, (pctile_cols[k]) := frank(get(prs_cols[k])) / .N]
}
# Genome-wide PGS percentiles (computed over the same merged sample as clusters);
# reused to define the GW top-33% comparison group within each direction subset.
dt[, pctile_gw_cad := frank(get(gw_cad_col)) / .N]
dt[, pctile_gw_t2d := frank(get(gw_t2d_col)) / .N]
# =====================================================================
# 3. Cluster assignment, Cox models, and plots per threshold
# =====================================================================

covar_terms <- c("age_at_first_dx", "age_at_first_dx2", "sex",
                 paste0("PC", 1:10))

all_levels <- c(clusters, "None")

# ---------------------------------------------------------------------
# Helper: covariate-adjusted cumulative hazard plot for a data subset.
# Fits a Cox model on dt_sub, predicts the adjusted cumulative hazard at
# median covariates for each cluster level present in the subset, and saves
# the curves to out_png. Used for the pooled (combined) plot at every
# threshold and for the per-direction (T2D-first / CAD-first) top33 plots.
# Reuses the enclosing-scope objects: covar_terms, all_levels,
# cluster_display, cluster_colors.
# ---------------------------------------------------------------------
plot_adjusted_cumhaz <- function(dt_sub, title, ylab, out_png, show_ci = TRUE,
                                 gw_pctile_col = NULL, gw_label = NULL,
                                 gw_thresh = NULL,
                                 xlab = "Years from first diagnosis") {
  # Keep only cluster levels with at least one individual in this subset, so
  # coxph/survfit do not choke on empty factor levels in direction subsets.
  present <- all_levels[vapply(all_levels,
                               function(l) sum(dt_sub$assigned_cluster == l) > 0,
                               logical(1))]
  if (!("None" %in% present)) {
    cat(sprintf("    WARNING: no 'None' reference individuals for %s; skipping\n", out_png))
    return(invisible(NULL))
  }

  ds <- copy(dt_sub)
  ds[, assigned_cluster := relevel(factor(as.character(assigned_cluster),
                                          levels = present), ref = "None")]

  form <- as.formula(paste(
    "Surv(time_years, event) ~",
    paste(c("assigned_cluster", covar_terms), collapse = " + ")
  ))
  fit <- tryCatch(coxph(form, data = ds), error = function(e) {
    cat(sprintf("    ERROR fitting Cox for %s: %s\n", out_png, e$message)); NULL
  })
  if (is.null(fit)) return(invisible(NULL))

  # Representative newdata: median covariates, vary only assigned_cluster
  median_covars <- ds[, lapply(.SD, median, na.rm = TRUE),
                      .SDcols = c("age_at_first_dx", "age_at_first_dx2", "sex",
                                  paste0("PC", 1:10))]
  newdata_list <- list()
  for (lev in present) {
    nd <- copy(median_covars)
    nd[, assigned_cluster := factor(lev, levels = present)]
    newdata_list[[lev]] <- nd
  }
  newdata <- rbindlist(newdata_list)
  newdata[, assigned_cluster := relevel(factor(assigned_cluster, levels = present),
                                        ref = "None")]

  fit_surv <- survfit(fit, newdata = newdata)

  # Extract per-cluster cumulative hazard curves (column i == present[i])
  plot_list <- list()
  for (i in seq_along(present)) {
    plot_list[[i]] <- data.table(
      time = fit_surv$time,
      cumhaz = fit_surv$cumhaz[, i],
      cumhaz_lower = -log(pmax(fit_surv$upper[, i], 1e-10)),
      cumhaz_upper = -log(pmax(fit_surv$lower[, i], 1e-10)),
      group = present[i]
    )
  }
  plot_dt <- rbindlist(plot_list)
  plot_dt[, group := factor(group, levels = present)]

  # Add origins (time 0, cumhaz 0)
  origins <- data.table(
    time = 0, cumhaz = 0, cumhaz_lower = 0, cumhaz_upper = 0,
    group = factor(present, levels = present)
  )
  plot_dt <- rbind(origins, plot_dt)

  # Relabel with biological names for plotting
  display_levels <- cluster_display[present]
  plot_dt[, group := factor(cluster_display[as.character(group)], levels = display_levels)]

  # Optional genome-wide PGS comparison curve. The GW top-33% group is NOT
  # mutually exclusive with the clusters, so it is fit from a separate Cox model
  # (binary GW-top vs rest) on the same subset and overlaid; the cluster curves
  # are left unchanged. Convert group to character so the new label rbinds and
  # the legend-ordering step below re-factors all groups together.
  plot_dt[, group := as.character(group)]
  if (!is.null(gw_pctile_col)) {
    ds_gw <- copy(dt_sub)
    ds_gw[, gw_grp := factor(ifelse(get(gw_pctile_col) > gw_thresh, "top", "rest"),
                             levels = c("rest", "top"))]
    form_gw <- as.formula(paste(
      "Surv(time_years, event) ~",
      paste(c("gw_grp", covar_terms), collapse = " + ")
    ))
    fit_gw <- tryCatch(coxph(form_gw, data = ds_gw), error = function(e) {
      cat(sprintf("    ERROR fitting GW Cox for %s: %s\n", out_png, e$message)); NULL
    })
    if (!is.null(fit_gw)) {
      # Adjusted cumulative hazard for the GW top-33% group at median covariates
      nd_gw <- copy(median_covars)
      nd_gw[, gw_grp := factor("top", levels = c("rest", "top"))]
      fs_gw <- survfit(fit_gw, newdata = nd_gw)
      gw_curve <- data.table(
        time = c(0, fs_gw$time),
        cumhaz = c(0, fs_gw$cumhaz),
        cumhaz_lower = 0, cumhaz_upper = 0,
        group = gw_label
      )
      plot_dt <- rbind(plot_dt, gw_curve)
    }
  }

  # Order the legend to match the visual order of the curves (descending final
  # cumulative hazard), so the legend reads top-to-bottom like the curves.
  curve_order <- plot_dt[, .(maxch = max(cumhaz, na.rm = TRUE)), by = group][
    order(-maxch), as.character(group)]
  plot_dt[, group := factor(as.character(group), levels = curve_order)]

  p <- ggplot(plot_dt, aes(x = time, color = group, fill = group))
  if (show_ci) {
    p <- p + geom_ribbon(aes(ymin = cumhaz_lower, ymax = cumhaz_upper),
                         alpha = 0.12, linewidth = 0, show.legend = FALSE)
  }
  p <- p +
    geom_step(aes(y = cumhaz), linewidth = 0.8) +
    scale_color_manual(values = cluster_colors, name = "PGS", drop = FALSE) +
    scale_fill_manual(values = cluster_colors, guide = "none", drop = FALSE) +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(xlim = c(0, 25), ylim = c(0, 0.8))

  ggsave(out_png, p, width = 10, height = 6, dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s\n", out_png))
  invisible(NULL)
}

cox_results_all <- list()

for (t_idx in seq_along(thresholds)) {
  thresh <- thresholds[t_idx]
  t_label <- thresh_labels[t_idx]
  pct_label <- sprintf("top %d%%", round(100 * (1 - thresh)))

  cat(sprintf("\n========== Threshold: %s (percentile > %.2f) ==========\n",
              pct_label, thresh))

  # --- Cluster assignment ---
  dt[, assigned_cluster := {
    pctiles <- numeric(length(pctile_cols))
    for (j in seq_along(pctile_cols)) pctiles[j] <- get(pctile_cols[j])
    eligible <- which(pctiles > thresh)
    if (length(eligible) == 0) {
      "None"
    } else {
      clusters[eligible[which.max(pctiles[eligible])]]
    }
  }, by = seq_len(nrow(dt))]

  dt[, assigned_cluster := factor(assigned_cluster, levels = all_levels)]

  # Report counts
  cat("\n  Cluster assignment counts:\n")
  assign_tab <- table(dt$assigned_cluster)
  for (lev in all_levels) {
    cat(sprintf("    %s: %d (%.1f%%)\n", lev, assign_tab[lev],
                100 * assign_tab[lev] / nrow(dt)))
  }

  # --- Cox model: single categorical model (None = reference) ---
  cat(sprintf("\n  --- Cox model (%s) ---\n", pct_label))

  # Set None as reference level
  dt[, assigned_cluster := relevel(assigned_cluster, ref = "None")]

  form <- as.formula(paste(
    "Surv(time_years, event) ~",
    paste(c("assigned_cluster", covar_terms), collapse = " + ")
  ))

  fit <- tryCatch(coxph(form, data = dt), error = function(e) {
    cat(sprintf("    ERROR: %s\n", e$message)); NULL
  })

  if (!is.null(fit)) {
    s <- summary(fit)
    conc <- s$concordance["C"]

    # Extract HR for each cluster vs None
    for (k in clusters) {
      coef_name <- paste0("assigned_cluster", k)
      if (!(coef_name %in% rownames(s$coefficients))) next

      n_in <- sum(dt$assigned_cluster == k)
      n_events_in <- sum(dt$assigned_cluster == k & dt$event == 1)
      hr <- exp(coef(fit)[coef_name])
      ci <- exp(confint(fit)[coef_name, ])
      pval <- s$coefficients[coef_name, "Pr(>|z|)"]

      ph_p <- tryCatch({
        zph <- cox.zph(fit)
        zph$table[coef_name, "p"]
      }, error = function(e) NA_real_)

      cox_results_all[[length(cox_results_all) + 1]] <- data.table(
        threshold = t_label,
        cluster = k,
        n_cluster = n_in,
        n_events_cluster = n_events_in,
        n_total = nrow(dt),
        n_events_total = sum(dt$event),
        HR = round(hr, 4),
        CI_lower = round(ci[1], 4),
        CI_upper = round(ci[2], 4),
        p_value = pval,
        c_index = round(conc, 4),
        ph_test_p = round(ph_p, 4)
      )

      cat(sprintf("    %s vs None: n=%d, events=%d, HR=%.3f (%.3f-%.3f), p=%.2e\n",
                  k, n_in, n_events_in, hr, ci[1], ci[2], pval))
    }
    cat(sprintf("    C-index: %.3f\n", conc))
  }

  # --- Adjusted cumulative hazard plot (from single model) ---
  cat(sprintf("\n  --- Adjusted cumulative hazard plot (%s) ---\n", pct_label))

  # Restore factor levels for plotting
  dt[, assigned_cluster := factor(assigned_cluster, levels = all_levels)]

  # Pooled (both comorbidity directions combined) — unchanged behavior
  plot_adjusted_cumhaz(
    dt,
    title   = sprintf("Adjusted cumulative hazard by cluster PGS (%s)", pct_label),
    ylab    = "Adjusted Cumulative Hazard\nof Developing Second Diagnosis",
    out_png = file.path(results_dir, sprintf("cumhaz_%s.png", t_label))
  )

  # --- Direction-split cumulative hazard (top33 and top10) ---
  # Cluster assignment is by PRS percentile and direction-independent, so we
  # split the already-assigned sample by which disease was diagnosed first.
  if (t_label %in% c("top33", "top10")) {
    directions <- list(
      # gw_pctile_col / gw_label overlay the outcome-matched genome-wide PGS
      # top curve at this threshold: GW CAD on T2D-first->CAD, GW T2D on CAD-first->T2D.
      list(dx = "T2D",
           xlab = "Years since Diagnosis of T2D",
           ylab = "Hazards of CAD among T2D Patients",
           title = sprintf("Adjusted cumulative hazard by cluster PGS (%s): T2D-first -> CAD", pct_label),
           file = sprintf("cumhaz_%s_t2d_to_cad.png", t_label),
           gw_pctile_col = "pctile_gw_cad", gw_label = "GW CAD"),
      list(dx = "CAD",
           xlab = "Years since Diagnosis of CAD",
           ylab = "Hazards of T2D among CAD Patients",
           title = sprintf("Adjusted cumulative hazard by cluster PGS (%s): CAD-first -> T2D", pct_label),
           file = sprintf("cumhaz_%s_cad_to_t2d.png", t_label),
           gw_pctile_col = "pctile_gw_t2d", gw_label = "GW T2D")
    )

    for (di in directions) {
      sub <- dt[first_dx == di$dx]
      cat(sprintf("\n  --- Direction %s-first (%s): n=%d, events=%d ---\n",
                  di$dx, pct_label, nrow(sub), sum(sub$event)))
      for (lev in all_levels) {
        nk <- sum(sub$assigned_cluster == lev)
        ek <- sum(sub$assigned_cluster == lev & sub$event == 1)
        cat(sprintf("    %s: n=%d, events=%d\n", lev, nk, ek))
      }
      # GW top comparison group at this threshold (not mutually exclusive with clusters)
      gw_top <- sub[[di$gw_pctile_col]] > thresh
      cat(sprintf("    %s: n=%d, events=%d\n", di$gw_label,
                  sum(gw_top), sum(gw_top & sub$event == 1)))
      plot_adjusted_cumhaz(
        sub,
        title   = di$title,
        ylab    = di$ylab,
        xlab    = di$xlab,
        out_png = file.path(results_dir, di$file),
        show_ci = FALSE,
        gw_pctile_col = di$gw_pctile_col,
        gw_label = di$gw_label,
        gw_thresh = thresh
      )
    }
  }
}

# =====================================================================
# 4. Save Cox results
# =====================================================================

cox_dt <- rbindlist(cox_results_all)
cox_file <- file.path(results_dir, "cox_results.csv")
fwrite(cox_dt, cox_file)
cat(sprintf("\n  Saved Cox results: %s (%d rows)\n", cox_file, nrow(cox_dt)))

# --- Print summary table ---
cat("\n--- Cox Results Summary ---\n")
cat(sprintf("%-8s %-8s %6s %6s %8s %15s %12s %8s\n",
            "Thresh", "Cluster", "n", "events", "HR", "95% CI", "p-value", "C-index"))
cat(paste(rep("-", 85), collapse = ""), "\n")
for (r in seq_len(nrow(cox_dt))) {
  row <- cox_dt[r]
  cat(sprintf("%-8s %-8s %6d %6d %8.3f (%6.3f-%6.3f) %12.2e %8.3f\n",
              row$threshold, row$cluster, row$n_cluster, row$n_events_cluster,
              row$HR, row$CI_lower, row$CI_upper, row$p_value, row$c_index))
}

cat("\n=== B2.1 Analysis complete ===\n")
