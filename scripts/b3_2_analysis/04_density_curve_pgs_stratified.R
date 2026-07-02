#!/usr/bin/env Rscript
# 04_density_curve_pgs_stratified.R
# Line plots showing cluster PGS decile (x-axis) vs a biomarker/risk (y-axis),
# with the sample stratified by genome-wide PRS (top 10% / all / bottom 10%).
#
#   Soft (default, META clusters) -> 2x2 grid:
#     TL) Glycemic (K10) PGS decile vs HbA1c (%),     stratified by GW T2D PRS
#     TR) Blood Pressure-Stature (K5) vs PCE risk (%), stratified by GW CAD PRS
#     BL) Blood Pressure-Stature (K5) vs LDL (mmol/L), stratified by GW CAD PRS
#     BR) Blood Pressure-Stature (K5) vs BMI (kg/m2),  stratified by GW CAD PRS
#
#   Hard assignment (--hard-assignment) -> original 1x2 grid:
#     A) Glucose 2 (K9) PGS decile vs HbA1c (%),  stratified by GW T2D PRS
#     B) Lipid (K4) PGS decile vs PCE risk (%),   stratified by GW CAD PRS
#
# Usage:
#   Rscript scripts/b3_2_analysis/04_density_curve_pgs_stratified.R
#   Rscript scripts/b3_2_analysis/04_density_curve_pgs_stratified.R --hard-assignment

library(data.table)
library(ggplot2)
library(yaml)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b3_2_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}
hard_assignment <- "--hard-assignment" %in% args
combined_mode   <- "--combined" %in% args
cad_cluster_grid <- "--cad-cluster-grid" %in% args
t2d_cluster_grid <- "--t2d-cluster-grid" %in% args

# Fixed seed so the random tie-breaking in decile assignment is reproducible
# (matters for any downstream counts derived from these deciles).
set.seed(1)

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir

# Soft (default) vs hard-assignment cluster PRS. The flag swaps both the PRS
# source file and the output directory; soft cluster PGS is the default.
if (hard_assignment) {
  cfg$prs_file <- "results/b1_analysis/prs_hard_assignment/prs/cluster_prs_all.tsv"
  results_dir <- file.path(results_dir, "prs_hard_assignment")
} else {
  cfg$prs_file <- "results/b1_analysis/prs/cluster_prs_all.tsv"
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B3.2: Cluster x GW-PRS Interaction Line Plots ===\n")
cat(sprintf("  Config:   %s\n", config_path))
cat(sprintf("  PRS file: %s\n", cfg$prs_file))
cat(sprintf("  Output:   %s\n", results_dir))
if (hard_assignment) cat("  ** HARD ASSIGNMENT MODE **\n")
if (combined_mode)   cat("  ** COMBINED-SCORE MODE (GW PRS vs GW PRS + cluster PGS) **\n")
if (cad_cluster_grid) cat("  ** CAD CLUSTER-GRID MODE (PCE risk: GW CAD vs GW CAD + each cluster cPGS) **\n")
if (t2d_cluster_grid) cat("  ** T2D CLUSTER-GRID MODE (HbA1c: GW T2D vs GW T2D + each cluster cPGS) **\n")

# ============================================================
# Load and merge data
# ============================================================
cat("\n--- Loading data ---\n")

prs_all <- fread(cfg$prs_file)
biomarker <- fread(cfg$biomarker_file,
                   select = c("FID", "IID", "HbA1c_corrected", "PCE_ASCVD",
                              "LDL_corrected", "BMI"))

dat <- merge(prs_all, biomarker, by = c("FID", "IID"), all.x = TRUE)

# Genome-wide PRS (independent of soft/hard cluster assignment). Uses the most
# recent trait-specific GW scores (GW_T2D = 768 variants, GW_CAD = 141 variants)
# from the B1 pipeline; `group` is dropped to avoid clashing with prs_all's.
gw_prs <- fread("results/b1_analysis/prs/genomewide_prs_all.tsv")
setnames(gw_prs, c("GW_T2D", "GW_CAD"), c("GW_T2D_PRS", "GW_CAD_PRS"))

dat <- merge(dat, gw_prs[, .(FID, IID, GW_T2D_PRS, GW_CAD_PRS)],
             by = c("FID", "IID"), all.x = TRUE)

eur_val <- dat[group == "eur_validation"]
cat(sprintf("  EUR validation: %d individuals\n", nrow(eur_val)))

# Unit conversion: HbA1c (mmol/mol) -> HbA1c (%)
eur_val[, HbA1c_pct := 0.09148 * HbA1c_corrected + 2.152]

dt_ntile <- function(x, n = 10L) {
  ranks <- frank(x, ties.method = "random", na.last = "keep")
  as.integer(ceiling(ranks / (sum(!is.na(x)) / n)))
}

# ============================================================
# GW PRS top/bottom 10% cutoffs
# ============================================================
t2d_q90 <- quantile(eur_val$GW_T2D_PRS, 0.90, na.rm = TRUE)
t2d_q10 <- quantile(eur_val$GW_T2D_PRS, 0.10, na.rm = TRUE)
cad_q90 <- quantile(eur_val$GW_CAD_PRS, 0.90, na.rm = TRUE)
cad_q10 <- quantile(eur_val$GW_CAD_PRS, 0.10, na.rm = TRUE)

cat(sprintf("  GW T2D PRS: bottom 10%% <= %.4f, top 10%% >= %.4f\n", t2d_q10, t2d_q90))
cat(sprintf("  GW CAD PRS: bottom 10%% <= %.4f, top 10%% >= %.4f\n", cad_q10, cad_q90))

# Cluster PGS top/bottom 10% cutoffs (used as the stratifier in soft mode)
k10_q90 <- quantile(eur_val$PRS_K10, 0.90, na.rm = TRUE)  # Glycemic
k10_q10 <- quantile(eur_val$PRS_K10, 0.10, na.rm = TRUE)
k5_q90  <- quantile(eur_val$PRS_K5,  0.90, na.rm = TRUE)   # Blood Pressure-Stature
k5_q10  <- quantile(eur_val$PRS_K5,  0.10, na.rm = TRUE)

cat(sprintf("  Glycemic (K10) PGS: bottom 10%% <= %.4f, top 10%% >= %.4f\n", k10_q10, k10_q90))
cat(sprintf("  BP-Stature (K5) PGS: bottom 10%% <= %.4f, top 10%% >= %.4f\n", k5_q10, k5_q90))

# ============================================================
# Aggregation: mean biomarker per decile of `decile_col`, within strata
# defined by the top/bottom 10% of `strat_col`
# ============================================================
# decile_col defines the x-axis decile (1-10); strat_col's extremes define three
# samples: all individuals, top 10% of strat_col, and bottom 10%. base_label
# names the x-axis variable (and the "all" line); strat_label names the stratifier.
build_strata_agg <- function(d, decile_col, value_col, strat_col,
                             strat_q90, strat_q10, base_label, strat_label) {
  dd <- d[!is.na(get(decile_col)) & !is.na(get(value_col)) & !is.na(get(strat_col))]
  dd[, decile := dt_ntile(get(decile_col), 10L)]

  lbl_top <- sprintf("%s + Top 10%% %s", base_label, strat_label)
  lbl_bot <- sprintf("%s + Bottom 10%% %s", base_label, strat_label)

  agg_all <- dd[, .(mean_val = mean(get(value_col)), n = .N), by = decile]
  agg_all[, stratum := base_label]
  agg_top <- dd[get(strat_col) >= strat_q90, .(mean_val = mean(get(value_col)), n = .N), by = decile]
  agg_top[, stratum := lbl_top]
  agg_bot <- dd[get(strat_col) <= strat_q10, .(mean_val = mean(get(value_col)), n = .N), by = decile]
  agg_bot[, stratum := lbl_bot]

  agg <- rbindlist(list(agg_all, agg_top, agg_bot))
  # Order: darkest (top) -> intermediate (all) -> lightest (bottom)
  agg[, stratum := factor(stratum, levels = c(lbl_top, base_label, lbl_bot))]

  cat(sprintf("  [%s decile ~ %s | strat by %s] all~%d  top10%%~%d  bot10%%~%d per decile\n",
              base_label, value_col, strat_label,
              round(mean(agg_all$n)), round(mean(agg_top$n)), round(mean(agg_bot$n))))
  agg
}

# Blue spectrum (GW T2D strata) and red spectrum (GW CAD strata).
# Vector order is dark, mid, light -> mapped to top/all/bottom strata levels.
blue_shades <- c("#004466", "#0072B2", "#7EC8E3")
red_shades  <- c("#8B3A00", "#D55E00", "#FFAB73")

stratum_linetypes <- c("solid", "dashed", "dotted")
stratum_shapes    <- c(16, 17, 15)

# ============================================================
# Combined-score aggregation: compare deciles of the GW PRS alone against
# deciles of an equal-weight combined score (GW PRS + cluster PGS).
# Both scores are standardized to mean 0 / SD 1 within the analysis sample
# before summing, so the combined score weights them equally. Both lines use
# the SAME individuals, re-ranked by each score (apples-to-apples top-decile).
# ============================================================
zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

build_combined_agg <- function(d, gw_col, pgs_col, value_col, gw_label, pgs_label,
                                threshold = NA_real_) {
  dd <- d[!is.na(get(gw_col)) & !is.na(get(pgs_col)) & !is.na(get(value_col))]
  dd[, gw_z := zscore(get(gw_col))]
  dd[, pgs_z := zscore(get(pgs_col))]
  dd[, combined_score := gw_z + pgs_z]

  combined_label <- sprintf("%s + %s", gw_label, pgs_label)

  dd[, dec_gw := dt_ntile(get(gw_col), 10L)]
  dd[, dec_cb := dt_ntile(combined_score, 10L)]

  # n_above = individuals exceeding the panel's clinical threshold within each decile
  agg_gw <- dd[, .(mean_val = mean(get(value_col)), n = .N,
                   n_above = sum(get(value_col) > threshold)), by = .(decile = dec_gw)]
  agg_gw[, stratum := gw_label]
  agg_cb <- dd[, .(mean_val = mean(get(value_col)), n = .N,
                   n_above = sum(get(value_col) > threshold)), by = .(decile = dec_cb)]
  agg_cb[, stratum := combined_label]

  agg <- rbindlist(list(agg_cb, agg_gw))
  # Combined first (darkest) -> GW PGS alone (lighter)
  agg[, stratum := factor(stratum, levels = c(combined_label, gw_label))]

  cat(sprintf("  [%s | y=%s] n/decile ~ %d (same %d individuals re-ranked by each score)\n",
              combined_label, value_col, round(mean(agg_gw$n)), nrow(dd)))
  agg
}

# ============================================================
# Top-tail aggregation: mean biomarker among everyone in the top X% of each
# score, for X = pcts (default 10..1). Cumulative (nested) tail groups. x is
# encoded so the largest % (10%) sits on the left and the smallest (1%) on the
# right. Same equal-weight combined score as build_combined_agg.
# ============================================================
build_tail_agg <- function(d, gw_col, pgs_col, value_col, gw_label, pgs_label,
                           threshold = NA_real_, pcts = 10:1) {
  dd <- d[!is.na(get(gw_col)) & !is.na(get(pgs_col)) & !is.na(get(value_col))]
  dd[, combined_score := zscore(get(gw_col)) + zscore(get(pgs_col))]

  combined_label <- sprintf("%s + %s", gw_label, pgs_label)

  # mean / n / n_above among the top p% of `score_col`
  tail_stats <- function(score_col, p) {
    cutoff <- quantile(dd[[score_col]], 1 - p / 100, na.rm = TRUE)
    grp <- dd[get(score_col) >= cutoff]
    data.table(pct = p,
               mean_val = mean(grp[[value_col]]),
               n = nrow(grp),
               n_above = sum(grp[[value_col]] > threshold))
  }

  agg_gw <- rbindlist(lapply(pcts, function(p) tail_stats(gw_col, p)))
  agg_gw[, stratum := gw_label]
  agg_cb <- rbindlist(lapply(pcts, function(p) tail_stats("combined_score", p)))
  agg_cb[, stratum := combined_label]

  agg <- rbindlist(list(agg_cb, agg_gw))
  # x position: 10% -> left (1), 1% -> right (max). Avoids a reversed axis.
  agg[, decile := (max(pcts) + 1L) - pct]
  agg[, stratum := factor(stratum, levels = c(combined_label, gw_label))]

  cat(sprintf("  [%s | y=%s | tail top %d%%..%d%%] %d individuals\n",
              combined_label, value_col, max(pcts), min(pcts), nrow(dd)))
  agg
}

# ============================================================
# Plot builder: turn an aggregated table (decile, mean_val, stratum) into a panel
# ============================================================
plot_agg <- function(agg, spec) {
  strata <- levels(agg$stratum)
  col_map <- setNames(spec$palette[seq_along(strata)], strata)
  lt_map  <- setNames(stratum_linetypes[seq_along(strata)], strata)
  sh_map  <- setNames(stratum_shapes[seq_along(strata)], strata)

  x_breaks <- if (!is.null(spec$x_breaks)) spec$x_breaks else 1:10
  x_labels <- if (!is.null(spec$x_labels)) spec$x_labels else waiver()

  p <- ggplot(agg, aes(x = decile, y = mean_val,
                       color = stratum, linetype = stratum, shape = stratum,
                       group = stratum))
  if (!is.null(spec$threshold) && !is.na(spec$threshold)) {
    p <- p + geom_hline(yintercept = spec$threshold, linetype = "dashed",
                        color = "grey50", linewidth = 0.5)
  }
  # Optional annotation: relative % difference in individuals exceeding the
  # clinical threshold (combined vs GW alone) at each x position.
  if (isTRUE(spec$annotate) && length(strata) == 2L) {
    wide <- dcast(agg, decile ~ stratum, value.var = "n_above")
    setnames(wide, c(strata[1], strata[2]), c("comb", "gw"))
    ytop <- agg[, .(y = max(mean_val)), by = decile]
    annot <- merge(wide, ytop, by = "decile")
    annot[, rel := 100 * (comb - gw) / gw]
    annot[, label := ifelse(gw > 0, sprintf("%+.0f%%", rel), "")]
    p <- p + geom_text(data = annot, aes(x = decile, y = y, label = label),
                       inherit.aes = FALSE, vjust = -0.9, size = 3, color = "grey20")
  }
  p <- p +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_color_manual(values = col_map, name = NULL) +
    scale_linetype_manual(values = lt_map, name = NULL) +
    scale_shape_manual(values = sh_map, name = NULL) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    labs(x = spec$x_label, y = spec$y_label, title = spec$title) +
    theme_minimal(base_size = 13) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text        = element_text(size = 10),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.direction = "vertical",
      legend.text      = element_text(size = 10),
      legend.key.width = unit(1.5, "cm")
    )

  # Optional shared y-axis across panels (coord_cartesian keeps all data)
  if (!is.null(spec$ylim)) p <- p + coord_cartesian(ylim = spec$ylim)
  p
}

# Strata panel (default/hard mode): GW PRS top10/all/bottom10 cluster-PGS strata
make_panel <- function(spec) {
  agg <- build_strata_agg(eur_val, spec$decile_col, spec$value_col, spec$strat_col,
                          spec$strat_q90, spec$strat_q10, spec$base_label, spec$strat_label)
  plot_agg(agg, spec)
}

# 5x2 cluster grid: mean `value_col` by score decile, GW PGS alone vs GW PGS +
# each META cluster cPGS (one panel per cluster K1..K10). The GW-alone reference
# line is identical across panels; only the combined cluster cPGS varies. Shared
# y-axis (comparable across clusters) and one shared bottom legend.
# Canonical META labels — source of truth: scripts/a1_analysis/06_ancestry_trait_barplots.R
cluster_grid_final <- function(gw_col, gw_label, value_col, y_label, threshold, palette) {
  cluster_labels <- c("Lpa", "Adiponectin", "Platelet", "SHBG",
                      "Blood Pressure-Stature", "Metabolic", "Triglycerides-HDL",
                      "ALP-LDL", "Obesity", "Glycemic")
  # Panel order = desired META display order (Glycemic -> Lpa); cluster_labels is
  # indexed by K number, so iterate over K numbers in the desired sequence.
  desired_k <- c(10, 9, 4, 2, 7, 8, 6, 3, 5, 1)
  specs <- lapply(desired_k, function(i) {
    list(type = "decile", gw_col = gw_col, pgs_col = sprintf("PRS_K%d", i),
         value_col = value_col, gw_label = gw_label, pgs_label = "cluster cPGS",
         palette = palette, threshold = threshold, title = cluster_labels[i],
         x_label = "Score Decile", y_label = y_label)
  })

  aggs <- lapply(specs, function(s)
    build_combined_agg(eur_val, s$gw_col, s$pgs_col, s$value_col,
                       s$gw_label, s$pgs_label, s$threshold))
  ylim <- range(unlist(lapply(aggs, function(a) a$mean_val)))
  specs <- lapply(specs, function(s) { s$ylim <- ylim; s })

  panels <- Map(plot_agg, aggs, specs)
  # One shared horizontal legend; target the bottom guide-box explicitly to avoid
  # the ggplot >=3.5 "multiple components" ambiguity.
  legend_src <- panels[[1]] + theme(legend.direction = "horizontal")
  legend <- get_plot_component(legend_src, "guide-box-bottom")
  panels_nolegend <- lapply(panels, function(p) p + theme(legend.position = "none"))
  grid <- plot_grid(plotlist = panels_nolegend, ncol = 2, align = "hv", axis = "tblr")
  plot_grid(grid, legend, ncol = 1, rel_heights = c(1, 0.05))
}

# Combined-score panel: GW PGS alone vs GW PGS + cluster PGS.
# type = "decile" -> deciles 1-10; type = "tail" -> cumulative top 10%..1%.
make_combined_panel <- function(spec) {
  if (!is.null(spec$type) && spec$type == "tail") {
    agg <- build_tail_agg(eur_val, spec$gw_col, spec$pgs_col, spec$value_col,
                          spec$gw_label, spec$pgs_label, spec$threshold)
  } else {
    agg <- build_combined_agg(eur_val, spec$gw_col, spec$pgs_col, spec$value_col,
                              spec$gw_label, spec$pgs_label, spec$threshold)
  }
  plot_agg(agg, spec)
}

# ============================================================
# Panel specifications (soft default = 2x2; hard = original 1x2)
# ============================================================
cat("\n--- Generating plots ---\n")

out_name <- "decile_line_cluster_gw_prs.png"

if (cad_cluster_grid) {
  # Mean PCE risk vs score decile, GW CAD PGS alone vs GW CAD PGS + each cluster cPGS
  final <- cluster_grid_final(gw_col = "GW_CAD_PRS", gw_label = "GW CAD PGS",
                              value_col = "PCE_ASCVD", y_label = "Mean PCE Risk (%)",
                              threshold = 7.5, palette = red_shades[c(1, 2)])
  fig_w <- 12; fig_h <- 22
  out_name <- "decile_line_cad_cluster_grid.png"
} else if (t2d_cluster_grid) {
  # Mean HbA1c vs score decile, GW T2D PGS alone vs GW T2D PGS + each cluster cPGS
  final <- cluster_grid_final(gw_col = "GW_T2D_PRS", gw_label = "GW T2D PGS",
                              value_col = "HbA1c_pct", y_label = "Mean HbA1c (%)",
                              threshold = 5.7, palette = blue_shades[c(1, 2)])
  fig_w <- 12; fig_h <- 22
  out_name <- "decile_line_t2d_cluster_grid.png"
} else if (combined_mode) {
  # Combined-score (META) clusters: K10 = Glycemic, K5 = Blood Pressure-Stature.
  # Two lines per panel: GW PGS alone vs equal-weight (GW PGS + cluster PGS).
  # Left column = diabetes/HbA1c, right column = CAD/PCE.
  #   Top row (A,B): x = score decile 1-10.
  #   Bottom row (C,D): x = cumulative top 10%..1% of the score.
  # Each x is annotated with the relative % difference in individuals exceeding
  # the clinical threshold (combined vs GW PGS alone).
  tail_pct_labels <- paste0(10:1, "%")
  specs <- list(
    list(type = "decile", gw_col = "GW_T2D_PRS", pgs_col = "PRS_K10", value_col = "HbA1c_pct",
         gw_label = "GW T2D PGS", pgs_label = "Glycemic cPGS",
         palette = blue_shades[c(1, 2)], threshold = 5.7,
         x_label = "Score Decile", y_label = "Mean HbA1c (%)"),
    list(type = "decile", gw_col = "GW_CAD_PRS", pgs_col = "PRS_K2", value_col = "PCE_ASCVD",
         gw_label = "GW CAD PGS", pgs_label = "Adiponectin cPGS",
         palette = red_shades[c(1, 2)], threshold = 7.5,
         x_label = "Score Decile", y_label = "Mean PCE Risk (%)"),
    list(type = "tail", gw_col = "GW_T2D_PRS", pgs_col = "PRS_K10", value_col = "HbA1c_pct",
         gw_label = "GW T2D PGS", pgs_label = "Glycemic cPGS",
         palette = blue_shades[c(1, 2)], threshold = 5.7,
         x_label = "Top % of score", x_breaks = 1:10, x_labels = tail_pct_labels,
         y_label = "Mean HbA1c (%)"),
    list(type = "tail", gw_col = "GW_CAD_PRS", pgs_col = "PRS_K2", value_col = "PCE_ASCVD",
         gw_label = "GW CAD PGS", pgs_label = "Adiponectin cPGS",
         palette = red_shades[c(1, 2)], threshold = 7.5,
         x_label = "Top % of score", x_breaks = 1:10, x_labels = tail_pct_labels,
         y_label = "Mean PCE Risk (%)")
  )
  panels <- lapply(specs, make_combined_panel)
  final <- plot_grid(plotlist = panels, ncol = 2, align = "hv", axis = "tblr",
                     labels = c("A", "B", "C", "D"), label_size = 14)
  fig_w <- 14; fig_h <- 12
  out_name <- "decile_line_combined_score.png"
} else if (hard_assignment) {
  # Hard-assignment clusters: K9 = Glucose 2, K4 = Lipid.
  # Original orientation: x = cluster PGS decile, strata = GW PRS top/bottom 10%.
  specs <- list(
    list(decile_col = "PRS_K9", value_col = "HbA1c_pct", strat_col = "GW_T2D_PRS",
         strat_q90 = t2d_q90, strat_q10 = t2d_q10, base_label = "Glucose 2",
         strat_label = "GW T2D PRS", palette = blue_shades, threshold = 5.7,
         x_label = "Glucose 2 PGS Decile", y_label = "Mean HbA1c (%)"),
    list(decile_col = "PRS_K4", value_col = "PCE_ASCVD", strat_col = "GW_CAD_PRS",
         strat_q90 = cad_q90, strat_q10 = cad_q10, base_label = "Lipid",
         strat_label = "GW CAD PRS", palette = red_shades, threshold = 7.5,
         x_label = "Lipid PGS Decile", y_label = "Mean PCE Risk (%)")
  )
  panels <- lapply(specs, make_panel)
  final <- plot_grid(plotlist = panels, ncol = 2, align = "h", axis = "tb",
                     labels = c("A", "B"), label_size = 14)
  fig_w <- 14; fig_h <- 7
} else {
  # Soft (META) clusters: K10 = Glycemic, K5 = Blood Pressure-Stature.
  # Flipped orientation: x = GW PRS decile, strata = cluster PGS top/bottom 10%.
  specs <- list(
    list(decile_col = "GW_T2D_PRS", value_col = "HbA1c_pct", strat_col = "PRS_K10",
         strat_q90 = k10_q90, strat_q10 = k10_q10, base_label = "GW T2D PRS",
         strat_label = "Glycemic PGS", palette = blue_shades, threshold = 5.7,
         x_label = "GW T2D PRS Decile", y_label = "Mean HbA1c (%)"),
    list(decile_col = "GW_CAD_PRS", value_col = "PCE_ASCVD", strat_col = "PRS_K5",
         strat_q90 = k5_q90, strat_q10 = k5_q10, base_label = "GW CAD PRS",
         strat_label = "Blood Pressure-Stature PGS", palette = red_shades, threshold = 7.5,
         x_label = "GW CAD PRS Decile", y_label = "Mean PCE Risk (%)"),
    list(decile_col = "GW_CAD_PRS", value_col = "LDL_corrected", strat_col = "PRS_K5",
         strat_q90 = k5_q90, strat_q10 = k5_q10, base_label = "GW CAD PRS",
         strat_label = "Blood Pressure-Stature PGS", palette = red_shades, threshold = 3.0,
         x_label = "GW CAD PRS Decile", y_label = "Mean LDL (mmol/L)"),
    list(decile_col = "GW_CAD_PRS", value_col = "BMI", strat_col = "PRS_K5",
         strat_q90 = k5_q90, strat_q10 = k5_q10, base_label = "GW CAD PRS",
         strat_label = "Blood Pressure-Stature PGS", palette = red_shades, threshold = 25,
         x_label = "GW CAD PRS Decile", y_label = expression("Mean BMI (kg/m"^2*")"))
  )
  panels <- lapply(specs, make_panel)
  final <- plot_grid(plotlist = panels, ncol = 2, align = "hv", axis = "tblr",
                     labels = c("A", "B", "C", "D"), label_size = 14)
  fig_w <- 14; fig_h <- 12
}

out_path <- file.path(results_dir, out_name)
ggsave(out_path, final, width = fig_w, height = fig_h, dpi = 300)
cat(sprintf("  Saved: %s\n", out_path))

cat("\n=== Done ===\n")
