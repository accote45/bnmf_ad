#!/usr/bin/env Rscript
#
# 09_cad_t2d_spectrum.R
# Show that META bNMF clusters lie on a T2D–CAD similarity spectrum.
#
# Four metrics per cluster mapped to scatter aesthetics:
#   x: % variants GWS for CAD    y: % variants GWS for T2D
#   color: h2_spectrum (log2 ratio of W-weighted mean z²)    size: cluster variant count
#
# Usage:
#   Rscript scripts/a1_analysis/09_cad_t2d_spectrum.R [--config path] [--ancestry META]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
})

# --- CLI args ----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/a1_config.yaml"
ancestry    <- "META"
if ("--config"   %in% args) config_path <- args[which(args == "--config") + 1]
if ("--ancestry" %in% args) ancestry    <- args[which(args == "--ancestry") + 1]

cfg         <- read_yaml(config_path)
base_dir    <- getwd()
results_dir <- cfg$results_dir
figures_dir <- file.path(results_dir, "figures")
anc_dir     <- file.path(results_dir, ancestry)
gws_p       <- 5e-8

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== T2D–CAD Spectrum: %s ===\n", ancestry))
cat(sprintf("  GWS threshold: %.0e\n", gws_p))

# --- Canonical cluster labels ------------------------------------------------

cluster_labels <- c(
  K1 = "Lpa",                            K2 = "Adiponectin",
  K3 = "Platelet",                       K4 = "SHBG",
  K5 = "Blood Pressure-Stature",         K6 = "Metabolic",
  K7 = "Triglycerides-HDL",              K8 = "ALP-LDL",
  K9 = "Obesity",                        K10 = "Glycemic"
)

# --- Step 1: Load W matrix ---------------------------------------------------

cat("\nStep 1: Loading W matrix...\n")
W <- fread(file.path(anc_dir, sprintf("W_matrix_%s.tsv", ancestry)))
var_ids <- W$VAR_ID
w_mat <- as.matrix(W[, -1, with = FALSE])
rownames(w_mat) <- var_ids
n_var <- nrow(w_mat)
n_k   <- ncol(w_mat)
cat(sprintf("  %d variants x %d clusters\n", n_var, n_k))

# --- Step 2: Reconstruct disease source membership ---------------------------

cat("\nStep 2: Querying reference GWAS for disease z-scores...\n")

ref_gwas <- cfg$ref_gwas[[ancestry]]
cad_keys <- grep("^CAD_", names(ref_gwas), value = TRUE)
t2d_keys <- grep("^T2D_", names(ref_gwas), value = TRUE)
cat(sprintf("  CAD sources: %d | T2D sources: %d\n", length(cad_keys), length(t2d_keys)))

query_ref_gwas <- function(keys, label) {
  z_list <- list()
  p_list <- list()
  for (i in seq_along(keys)) {
    fp <- file.path(base_dir, ref_gwas[[keys[i]]])
    if (!file.exists(fp)) {
      cat(sprintf("    [%s] MISSING: %s\n", label, basename(fp)))
      next
    }
    ss <- fread(fp, select = c("VAR_ID", "P_VALUE", "BETA", "SE"))
    ss <- ss[VAR_ID %in% var_ids & !is.na(BETA) & !is.na(SE) & SE > 0]
    ss[, abs_z := abs(BETA / SE)]
    ss[, abs_beta := abs(BETA)]
    z_list[[keys[i]]] <- ss[, .(VAR_ID, abs_z, abs_beta)]
    p_list[[keys[i]]] <- ss[, .(VAR_ID, P_VALUE)]
    if (i %% 5 == 0 || i == length(keys))
      cat(sprintf("    [%s %d/%d] %s: %d variants\n", label, i, length(keys), keys[i], nrow(ss)))
  }
  z_all <- rbindlist(z_list)
  p_all <- rbindlist(p_list)
  # For each variant, keep the row with max |z| (and its paired BETA)
  z_max <- z_all[z_all[, .I[which.max(abs_z)], by = VAR_ID]$V1]
  setnames(z_max, c("VAR_ID", "max_abs_z", "max_z_beta"))
  p_min <- p_all[, .(min_p = min(P_VALUE, na.rm = TRUE)), by = VAR_ID]
  merge(z_max, p_min, by = "VAR_ID")
}

# Cache intermediate results to avoid re-reading 36 GWAS files on reruns
cache_file <- file.path(anc_dir, sprintf("disease_stats_cache_v2_%s.tsv", ancestry))
if (file.exists(cache_file)) {
  cat("  Loading cached disease stats (v2)...\n")
  disease_dt <- fread(cache_file)
} else {
  cad_stats <- query_ref_gwas(cad_keys, "CAD")
  t2d_stats <- query_ref_gwas(t2d_keys, "T2D")

  setnames(cad_stats, c("max_abs_z", "min_p", "max_z_beta"),
           c("z_CAD", "p_CAD", "beta_CAD"))
  setnames(t2d_stats, c("max_abs_z", "min_p", "max_z_beta"),
           c("z_T2D", "p_T2D", "beta_T2D"))

  disease_dt <- merge(
    data.table(VAR_ID = var_ids),
    cad_stats, by = "VAR_ID", all.x = TRUE
  )
  disease_dt <- merge(disease_dt, t2d_stats, by = "VAR_ID", all.x = TRUE)

  disease_dt[is.na(z_CAD), z_CAD := 0]
  disease_dt[is.na(z_T2D), z_T2D := 0]
  disease_dt[is.na(p_CAD), p_CAD := 1]
  disease_dt[is.na(p_T2D), p_T2D := 1]
  disease_dt[is.na(beta_CAD), beta_CAD := 0]
  disease_dt[is.na(beta_T2D), beta_T2D := 0]

  fwrite(disease_dt, cache_file, sep = "\t")
  cat(sprintf("  Cached disease stats (v2): %s\n", cache_file))
}

disease_dt[, gws_CAD := p_CAD < gws_p]
disease_dt[, gws_T2D := p_T2D < gws_p]
disease_dt[, category := fcase(
  gws_CAD & gws_T2D,  "Both",
  gws_CAD & !gws_T2D, "CAD-only",
  !gws_CAD & gws_T2D, "T2D-only",
  default = "Neither"
)]

cat(sprintf("\n  Variant classification:\n"))
cat(sprintf("    CAD-only:  %d\n", sum(disease_dt$category == "CAD-only")))
cat(sprintf("    T2D-only:  %d\n", sum(disease_dt$category == "T2D-only")))
cat(sprintf("    Both:      %d\n", sum(disease_dt$category == "Both")))
cat(sprintf("    Neither:   %d\n", sum(disease_dt$category == "Neither")))

# --- Step 3: Spectrum position per cluster -----------------------------------

cat("\nStep 3: Computing spectrum positions...\n")

z_cad_vec    <- disease_dt[match(var_ids, VAR_ID), z_CAD]
z_t2d_vec    <- disease_dt[match(var_ids, VAR_ID), z_T2D]
beta_cad_vec <- disease_dt[match(var_ids, VAR_ID), beta_CAD]
beta_t2d_vec <- disease_dt[match(var_ids, VAR_ID), beta_T2D]

spectrum_dt <- data.table(cluster = colnames(w_mat))

# Original z-based scores (kept for comparison)
spectrum_dt[, `:=`(
  score_CAD = sapply(seq_len(n_k), function(k) {
    sum(w_mat[, k] * z_cad_vec) / sum(w_mat[, k])
  }),
  score_T2D = sapply(seq_len(n_k), function(k) {
    sum(w_mat[, k] * z_t2d_vec) / sum(w_mat[, k])
  })
)]
spectrum_dt[, log2_ratio := log2(pmax(score_CAD, 1e-3) / pmax(score_T2D, 1e-3))]

# h2-proxy scores: W-weighted mean z² (chi-squared), precision-adjusted
# Uses z = BETA/SE instead of raw BETA to prevent imprecise large-effect
# variants from dominating (e.g. HK1 variant in the K3 cluster had beta_CAD=0.89
# but p_CAD=0.06, contributing 75% of the old BETA²-based CAD score)
spectrum_dt[, `:=`(
  h2_proxy_CAD = sapply(seq_len(n_k), function(k) {
    sum(w_mat[, k] * z_cad_vec^2) / sum(w_mat[, k])
  }),
  h2_proxy_T2D = sapply(seq_len(n_k), function(k) {
    sum(w_mat[, k] * z_t2d_vec^2) / sum(w_mat[, k])
  })
)]
global_mean_z2_CAD <- mean(z_cad_vec^2)
global_mean_z2_T2D <- mean(z_t2d_vec^2)
spectrum_dt[, h2_spectrum := log2(
  pmax(h2_proxy_CAD / global_mean_z2_CAD, 1e-10) /
  pmax(h2_proxy_T2D / global_mean_z2_T2D, 1e-10)
)]
spectrum_dt[, label := cluster_labels[cluster]]

cat("  Cluster spectrum positions (z²-proxy):\n")
for (i in seq_len(nrow(spectrum_dt))) {
  r <- spectrum_dt[i]
  cat(sprintf("    %s (%s): z2_CAD=%.4e, z2_T2D=%.4e, h2_spectrum=%.3f (old log2=%.3f)\n",
              r$cluster, r$label, r$h2_proxy_CAD, r$h2_proxy_T2D, r$h2_spectrum, r$log2_ratio))
}

# --- Step 4: Donut proportions per cluster -----------------------------------

cat("\nStep 4: Computing donut proportions...\n")

max_cluster <- colnames(w_mat)[max.col(w_mat, ties.method = "first")]
disease_dt[, assigned_cluster := max_cluster[match(VAR_ID, var_ids)]]

donut_dt <- disease_dt[category != "Neither",
  .N, by = .(assigned_cluster, category)]
donut_dt[, total := sum(N), by = assigned_cluster]
donut_dt[, prop := N / total]

cat("  Variants per cluster (max-loading assignment, excl. Neither):\n")
cluster_totals <- donut_dt[, .(n = sum(N)), by = assigned_cluster][order(assigned_cluster)]
for (i in seq_len(nrow(cluster_totals))) {
  r <- cluster_totals[i]
  cat(sprintf("    %s: %d variants\n", r$assigned_cluster, r$n))
}

# --- Step 5: Save scores table -----------------------------------------------

out_table <- merge(spectrum_dt,
  dcast(donut_dt, assigned_cluster ~ category, value.var = "N", fill = 0),
  by.x = "cluster", by.y = "assigned_cluster", all.x = TRUE)
out_file <- file.path(anc_dir, sprintf("cad_t2d_spectrum_scores_%s.tsv", ancestry))
fwrite(out_table, out_file, sep = "\t")
cat(sprintf("\nSaved scores table: %s\n", out_file))

# --- Step 6: Generate scatterplot ---------------------------------------------

cat("\nStep 6: Generating spectrum scatterplot...\n")

# Per-cluster scatter metrics from donut counts
scatter_dt <- donut_dt[, .(
  n_CAD_only = sum(N[category == "CAD-only"]),
  n_T2D_only = sum(N[category == "T2D-only"]),
  n_Both     = sum(N[category == "Both"]),
  n_total    = sum(N)
), by = assigned_cluster]

scatter_dt[, pct_CAD := 100 * (n_CAD_only + n_Both) / n_total]
scatter_dt[, pct_T2D := 100 * (n_T2D_only + n_Both) / n_total]

# Fold-enrichment over background: corrects for T2D having higher base rate
global_CAD_rate <- sum(disease_dt$gws_CAD) / nrow(disease_dt)
global_T2D_rate <- sum(disease_dt$gws_T2D) / nrow(disease_dt)
scatter_dt[, enrich_CAD := (pct_CAD / 100) / global_CAD_rate]
scatter_dt[, enrich_T2D := (pct_T2D / 100) / global_T2D_rate]

cat(sprintf("  Global GWS rates: CAD=%.1f%%, T2D=%.1f%%\n",
            global_CAD_rate * 100, global_T2D_rate * 100))

# Merge with spectrum positions (W-weighted mean log2 ratio) and labels
scatter_dt <- merge(scatter_dt, spectrum_dt[, .(cluster, log2_ratio, h2_spectrum, label)],
                    by.x = "assigned_cluster", by.y = "cluster")

cat("  Scatter data (fold-enrichment):\n")
for (i in seq_len(nrow(scatter_dt))) {
  r <- scatter_dt[i]
  cat(sprintf("    %s (%s): CAD=%.2fx, T2D=%.2fx, log2=%.2f, n=%d\n",
              r$assigned_cluster, r$label, r$enrich_CAD, r$enrich_T2D,
              r$log2_ratio, r$n_total))
}

p <- ggplot(scatter_dt, aes(x = enrich_CAD, y = enrich_T2D)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_point(aes(fill = h2_spectrum, size = n_total),
             shape = 21, color = "black", stroke = 0.4, alpha = 0.85) +
  geom_text_repel(
    aes(label = label),
    size = 3.5, fontface = "bold",
    box.padding = 0.6, point.padding = 0.4,
    min.segment.length = 0.2, segment.color = "black",
    max.overlaps = Inf
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = expression(h^2 ~ "Spectrum")
  ) +
  scale_size_continuous(
    range = c(3, 14),
    name = "Cluster size (SNVs)"
  ) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(
    x = "CAD Similarity",
    y = "T2D Similarity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.line.x = element_line(arrow = arrow(length = unit(0.25, "cm"), ends = "last", type = "closed")),
    axis.line.y = element_line(arrow = arrow(length = unit(0.25, "cm"), ends = "last", type = "closed")),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_colorbar(order = 2),
    size = guide_legend(order = 1)
  )

ggsave(file.path(figures_dir, sprintf("cad_t2d_spectrum_%s.png", ancestry)),
       p, width = 9, height = 7, dpi = 300)
cat(sprintf("Saved: %s\n", file.path(figures_dir,
            sprintf("cad_t2d_spectrum_%s.png", ancestry))))

# --- Step 7: 1D Cartoon spectrum -----------------------------------------------

cat("\nStep 7: Generating 1D spectrum cartoon...\n")

cartoon_dt <- copy(spectrum_dt)
cartoon_dt <- cartoon_dt[order(h2_spectrum)]
cartoon_dt[, rank := seq_len(.N)]

# Alternate label position above/below
cartoon_dt[, y_label := fifelse(rank %% 2 == 1, 0.45, -0.45)]
cartoon_dt[, y_nudge := fifelse(rank %% 2 == 1, 0.12, -0.12)]

p_cartoon <- ggplot(cartoon_dt, aes(x = rank, y = 0)) +
  annotate("segment", x = 0.2, xend = n_k + 0.8, y = 0, yend = 0,
           arrow = arrow(ends = "both", length = unit(0.25, "cm"), type = "closed"),
           linewidth = 0.7, color = "grey40") +
  geom_point(aes(fill = h2_spectrum), shape = 21, size = 7,
             color = "black", stroke = 0.6) +
  geom_segment(aes(x = rank, xend = rank, y = y_nudge, yend = y_label - sign(y_label) * 0.06),
               linewidth = 0.3, color = "grey50") +
  geom_text(aes(y = y_label, label = label), size = 3.2, fontface = "bold") +
  annotate("text", x = 0.0, y = 0, label = "T2D", fontface = "bold",
           size = 5, color = "#2166AC", hjust = 1) +
  annotate("text", x = n_k + 1.0, y = 0, label = "CAD", fontface = "bold",
           size = 5, color = "#B2182B", hjust = 0) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "Spectrum Score"
  ) +
  scale_x_continuous(expand = expansion(add = 1.8)) +
  coord_cartesian(ylim = c(-0.7, 0.7)) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2.5, "cm"),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    plot.margin  = margin(15, 25, 10, 25),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(figures_dir, sprintf("cad_t2d_spectrum_1d_%s.png", ancestry)),
       p_cartoon, width = 12, height = 4, dpi = 300)
cat(sprintf("Saved: %s\n", file.path(figures_dir,
            sprintf("cad_t2d_spectrum_1d_%s.png", ancestry))))

cat("\nDone.\n")
