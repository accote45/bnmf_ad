#!/usr/bin/env Rscript
#
# Render standalone horizontal "-logP" color-bar legends for the A2 heatmaps so
# they can be cropped into a manually-stitched figure. The per-cluster heatmaps
# draw with the legend suppressed (draw(..., show_heatmap_legend = FALSE)), so
# this regenerates just the legend bar for each analysis.
#
# Uses the SAME primitives as the heatmaps (ComplexHeatmap::Legend +
# circlize::colorRamp2) so the bars are pixel-faithful. The palette / range / at
# specs below are transcribed (with citations) from the three source scripts.
# IMPORTANT: if a source script's col_fun or heatmap_legend_param changes, update
# the matching spec here.
#
# Usage:
#   module load gcc/14.2.0 R/4.2.0
#   Rscript scripts/a2_analysis/a2_make_legends.R

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(data.table)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
out_dir  <- file.path(base_dir, "results/a2_analysis/legends")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# GTEx uses a data-derived range: max(ceiling(max(log10_p)), 1) over its results
# (a2_gtex_heatmaps.R:303-305). Read it from the cached results so the legend
# tracks the actual figures.
gtex_csv <- file.path(base_dir, "results/a2_analysis/gtex_ttest_results_nearest.csv")
gtex_log10p <- fread(gtex_csv)$log10_p
gtex_range  <- max(ceiling(max(gtex_log10p, na.rm = TRUE)), 1)
message(sprintf("GTEx data-derived color range: 0 to %d", gtex_range))

# The source scripts set the TS/catlas legend `at = seq(-log10(1), -log10(1e-4),
# length = 5)` = 0,1,2,3,4 (the permutation's resolvable p floor at n_perm=1e4).
# But the col_fun spans 0-10 and cells with p==0 clamp to the darkest color, so a
# 0-4 bar would misrepresent the figure. We render the FULL col_fun range (0-10,
# ticks every 2) so the legend matches the actual cell colors. (For a literal
# match to the scripts' `at`, use seq(0, 4, length = 5) instead.)
ticks_full_10 <- seq(0, 10, by = 2)

# Per-analysis legend specs (transcribed from the source heatmap scripts)
specs <- list(
  gtex = list(
    # a2_gtex_heatmaps.R:308-311 (palette) + :346 (legend param, default `at`)
    palette = c("#F0FAF0", "#B2E2B2", "#66C266", "#2E8B57", "#1B4D2E"),
    range   = gtex_range,
    at      = NULL
  ),
  tabula_sapiens = list(
    # a2_tabula_sapiens_heatmaps.R:320-323 (palette) + :362-367 (legend param)
    palette = c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"),
    range   = 10,
    at      = ticks_full_10
  ),
  catlas = list(
    # a2_catlas_pancreas_heatmaps.R:281-284 (palette) + :315-320 (legend param)
    # (catlas_esophagus uses the identical purple scale)
    palette = c("#F2F0F7", "#CBC9E2", "#9E9AC8", "#756BB1", "#54278F"),
    range   = 10,
    at      = ticks_full_10
  )
)

# Build one horizontal "-logP" legend, matching the heatmap construction exactly:
#   colorRamp2(seq(0, range, length.out = 10), colorRampPalette(pal)(10))
make_legend <- function(spec) {
  col_fun <- colorRamp2(
    seq(0, spec$range, length.out = 10),
    colorRampPalette(spec$palette)(10)
  )
  args <- list(
    col_fun        = col_fun,
    title          = "-logP",
    direction      = "horizontal",
    legend_width   = unit(4, "cm"),
    title_position = "topcenter"
  )
  if (!is.null(spec$at)) args$at <- spec$at
  do.call(Legend, args)
}

for (nm in names(specs)) {
  lgd <- make_legend(specs[[nm]])
  png_path <- file.path(out_dir, sprintf("%s_legend.png", nm))
  # Small canvas sized to the legend; transparent background for clean overlay.
  png(png_path, width = 2.4, height = 1.0, units = "in", res = 600, bg = "transparent")
  grid.newpage()
  draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
  dev.off()
  message(sprintf("  Saved %s", basename(png_path)))
}

# Also a stacked panel of all three for convenience.
combined_path <- file.path(out_dir, "a2_legends_combined.png")
png(combined_path, width = 2.6, height = 2.6, units = "in", res = 600, bg = "transparent")
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1)))
i <- 1
for (nm in names(specs)) {
  pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
  draw(make_legend(specs[[nm]]), x = unit(0.5, "npc"), y = unit(0.5, "npc"))
  popViewport()
  i <- i + 1
}
popViewport()
dev.off()
message(sprintf("  Saved %s", basename(combined_path)))

message("\nDone.")
