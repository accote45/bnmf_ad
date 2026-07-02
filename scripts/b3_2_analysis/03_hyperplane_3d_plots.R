#!/usr/bin/env Rscript
# 03_hyperplane_3d_plots.R
# 3D GAM surface plots: PGS percentile as a function of PCE risk and HbA1c
# per bNMF cluster. Renders a smooth surface colored by predicted PGS percentile.
#
# Usage:
#   Rscript scripts/b3_2_analysis/03_hyperplane_3d_plots.R --hard-assignment
#   Rscript scripts/b3_2_analysis/03_hyperplane_3d_plots.R --config config/b3_2_config.yaml --hard-assignment
#   Rscript scripts/b3_2_analysis/03_hyperplane_3d_plots.R --hard-assignment --theta 45 --phi 30

library(data.table)
library(yaml)
library(mgcv)
library(ragg)
library(fields)

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b3_2_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}
hard_assignment <- "--hard-assignment" %in% args

theta_val <- 35
if ("--theta" %in% args) {
  theta_val <- as.numeric(args[which(args == "--theta") + 1])
}
phi_val <- 25
if ("--phi" %in% args) {
  phi_val <- as.numeric(args[which(args == "--phi") + 1])
}

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
if (hard_assignment) {
  results_dir <- file.path(results_dir, "prs_hard_assignment")
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B3.2: 3D Hyperplane Plots ===\n")
cat(sprintf("  Config: %s\n", config_path))
cat(sprintf("  Output: %s\n", results_dir))
cat(sprintf("  Viewing angle: theta=%g, phi=%g\n", theta_val, phi_val))
if (hard_assignment) cat("  ** HARD ASSIGNMENT MODE **\n")

# ============================================================
# Load and merge data
# ============================================================
cat("\n--- Loading data ---\n")

prs_all <- fread(cfg$prs_file)
biomarker <- fread(cfg$biomarker_file,
                   select = c("FID", "IID", "PCE_ASCVD", "HbA1c_corrected"))

prs_cols <- grep("^PRS_K", colnames(prs_all), value = TRUE)
n_clusters <- length(prs_cols)
cat(sprintf("  PRS: %d individuals, %d clusters\n", nrow(prs_all), n_clusters))
cat(sprintf("  Biomarkers: %d individuals\n", nrow(biomarker)))

dat <- merge(prs_all, biomarker, by = c("FID", "IID"), all.x = TRUE)

eur_val <- dat[group == "eur_validation"]
cat(sprintf("  EUR validation: %d individuals\n", nrow(eur_val)))

eur_val[, HbA1c_pct := 0.09148 * HbA1c_corrected + 2.152]

if (!is.null(cfg$cluster_labels)) {
  cluster_labels <- unlist(cfg$cluster_labels)
} else {
  cluster_labels <- setNames(gsub("PRS_", "", prs_cols), gsub("PRS_", "", prs_cols))
}

# ============================================================
# Helper: 3D hyperplane plot
# ============================================================
draw_3d_hyperplane <- function(data, prs_col, k_id, cluster_name,
                               theta, phi, grid_n = 50,
                               cex.lab = 0.9, cex.axis = 0.7,
                               cex.main = 1.3, cex.dot = 1.2) {
  dt <- copy(data)
  dt[, pgs_pctl := frank(get(prs_col), ties.method = "average") / .N * 100]
  sub <- dt[!is.na(PCE_ASCVD) & !is.na(HbA1c_pct) & !is.na(pgs_pctl)]

  fit <- gam(pgs_pctl ~ s(PCE_ASCVD, HbA1c_pct, k = 20), data = sub)

  xg <- seq(quantile(sub$PCE_ASCVD, 0.01),
            quantile(sub$PCE_ASCVD, 0.99), length.out = grid_n)
  yg <- seq(quantile(sub$HbA1c_pct, 0.01),
            quantile(sub$HbA1c_pct, 0.99), length.out = grid_n)
  grid <- CJ(PCE_ASCVD = xg, HbA1c_pct = yg)
  grid[, pred := predict(fit, newdata = grid)]
  zmat <- matrix(grid$pred, nrow = grid_n, ncol = grid_n, byrow = FALSE)

  zfacet <- (zmat[-1, -1] + zmat[-grid_n, -1] +
             zmat[-1, -grid_n] + zmat[-grid_n, -grid_n]) / 4
  color_ramp <- designer.colors(256, c("#2166AC", "#F7F7F7", "#B2182B"))
  z_range <- range(zfacet, na.rm = TRUE)
  z_scaled <- (zfacet - z_range[1]) / max(z_range[2] - z_range[1], 1e-6)
  col_idx <- pmin(pmax(round(z_scaled * 255) + 1, 1), 256)
  col_mat <- matrix(color_ramp[col_idx], nrow = nrow(zfacet))

  zlim <- range(zmat, na.rm = TRUE)

  pmat <- persp(xg, yg, zmat,
                theta = theta, phi = phi,
                zlim = zlim,
                xlab = "PCE Risk (%)",
                ylab = "HbA1c (%)",
                zlab = "PGS Percentile",
                col = col_mat,
                shade = 0.3,
                border = NA,
                ticktype = "detailed",
                nticks = 5,
                cex.lab = cex.lab,
                cex.axis = cex.axis)

  x_ticks <- pretty(xg, n = 5)
  x_ticks <- x_ticks[x_ticks >= min(xg) & x_ticks <= max(xg)]
  y_ticks <- pretty(yg, n = 5)
  y_ticks <- y_ticks[y_ticks >= min(yg) & y_ticks <= max(yg)]
  z_ticks <- pretty(zlim, n = 5)
  z_ticks <- z_ticks[z_ticks >= zlim[1] & z_ticks <= zlim[2]]

  for (xt in x_ticks) {
    l <- trans3d(xt, range(yg), zlim[1], pmat)
    lines(l, col = "grey70", lwd = 0.5)
    l <- trans3d(xt, min(yg), c(zlim[1], zlim[2]), pmat)
    lines(l, col = "grey70", lwd = 0.5)
  }
  for (yt in y_ticks) {
    l <- trans3d(range(xg), yt, zlim[1], pmat)
    lines(l, col = "grey70", lwd = 0.5)
    l <- trans3d(min(xg), yt, c(zlim[1], zlim[2]), pmat)
    lines(l, col = "grey70", lwd = 0.5)
  }
  for (zt in z_ticks) {
    l <- trans3d(min(xg), range(yg), zt, pmat)
    lines(l, col = "grey70", lwd = 0.5)
    l <- trans3d(range(xg), max(yg), zt, pmat)
    lines(l, col = "grey70", lwd = 0.5)
  }

  corner_z <- zmat[grid_n, grid_n]
  pt <- trans3d(max(xg), max(yg), corner_z, pmat)
  points(pt, pch = 16, col = "black", cex = cex.dot)

  title(main = paste0(cluster_name, " (", k_id, ")"),
        cex.main = cex.main, font.main = 1)
}

make_3d_hyperplane <- function(data, prs_col, k_id, cluster_name,
                               out_path, theta, phi, grid_n = 50) {
  agg_png(out_path, width = 2400, height = 1800, res = 300)
  par(mar = c(1, 1, 3, 1), bg = "white")
  draw_3d_hyperplane(data, prs_col, k_id, cluster_name, theta, phi, grid_n)
  dev.off()
}

# ============================================================
# Generate plots for each cluster
# ============================================================
cat("\n--- Generating 3D hyperplane plots ---\n")

for (prs_col in prs_cols) {
  k_id <- gsub("PRS_", "", prs_col)
  cluster_name <- cluster_labels[[k_id]]
  safe_name <- gsub("[^A-Za-z0-9]", "", cluster_name)
  out_path <- file.path(results_dir,
                        sprintf("hyperplane_3d_%s_%s.png", k_id, safe_name))

  cat(sprintf("  %s (%s)...", k_id, cluster_name))
  make_3d_hyperplane(eur_val, prs_col, k_id, cluster_name,
                     out_path, theta_val, phi_val)
  cat(sprintf(" saved: %s\n", basename(out_path)))
}

# ============================================================
# Combined 2x2 figure (Lipid, Glucose 1, MetSyn-Inflam, Lipid-Liver)
# ============================================================
cat("\n--- Generating combined 2x2 hyperplane figure ---\n")

combined_panels <- list(
  list(k = "K4",  label = "a"),
  list(k = "K3",  label = "b"),
  list(k = "K7",  label = "c"),
  list(k = "K10", label = "d")
)

out_combined <- file.path(results_dir, "hyperplane_3d_combined_2x2.png")
agg_png(out_combined, width = 4000, height = 3200, res = 300)
par(mfrow = c(2, 2), mar = c(1, 0.5, 3, 0.5), oma = c(0, 0, 0, 0), bg = "white")

for (panel in combined_panels) {
  prs_col <- paste0("PRS_", panel$k)
  k_id <- panel$k
  cluster_name <- cluster_labels[[k_id]]
  draw_3d_hyperplane(eur_val, prs_col, k_id, cluster_name,
                     theta_val, phi_val,
                     cex.lab = 0.8, cex.axis = 0.6,
                     cex.main = 1.1, cex.dot = 1.0)
  mtext(panel$label, side = 3, adj = 0, line = 1.5, cex = 1.2, font = 2)
}

dev.off()
cat(sprintf("  saved: %s\n", basename(out_combined)))

cat("\n=== Done ===\n")
