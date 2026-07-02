#!/usr/bin/env Rscript
# 01_decile_biomarker_risk.R
# Decile-stratified biomarker and cardiovascular risk plots per bNMF cluster
# and genome-wide PRS (sanity check).
# For each score: 2x2 grid of LDL, HbA1c, PCE, PREVENT by PGS decile,
# sex-stratified (linetype), with clinical thresholds.
#
# Usage:
#   Rscript scripts/b3_2_analysis/01_decile_biomarker_risk.R --config config/b3_2_config.yaml --hard-assignment

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

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
if (hard_assignment) {
  results_dir <- file.path(results_dir, "prs_hard_assignment")
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B3.2: Decile Biomarker & Risk Plots ===\n")
cat(sprintf("  Config: %s\n", config_path))
cat(sprintf("  Output: %s\n", results_dir))
if (hard_assignment) cat("  ** HARD ASSIGNMENT MODE **\n")

# ============================================================
# Load and merge data
# ============================================================
cat("\n--- Loading data ---\n")

prs_all <- fread(cfg$prs_file)
biomarker <- fread(cfg$biomarker_file)

prs_cols <- grep("^PRS_K", colnames(prs_all), value = TRUE)
n_clusters <- length(prs_cols)
cat(sprintf("  PRS: %d individuals, %d clusters\n", nrow(prs_all), n_clusters))
cat(sprintf("  Biomarkers: %d individuals\n", nrow(biomarker)))

dat <- merge(prs_all, biomarker, by = c("FID", "IID"), all.x = TRUE)

# Load genome-wide PRS for sanity check
gw_t2d <- fread("results/a0_analysis/prs_ct/prsice_output/T2D.all_score")
gw_cad <- fread("results/a0_analysis/prs_ct/prsice_output/CAD.all_score")

# Select best-fit thresholds (use Pt_0.01 as reasonable default)
# Find available threshold columns
t2d_thresh_cols <- grep("^Pt_", colnames(gw_t2d), value = TRUE)
cad_thresh_cols <- grep("^Pt_", colnames(gw_cad), value = TRUE)
cat(sprintf("  GW T2D thresholds: %s\n", paste(t2d_thresh_cols, collapse = ", ")))
cat(sprintf("  GW CAD thresholds: %s\n", paste(cad_thresh_cols, collapse = ", ")))

# Use Pt_1 (all SNPs) as the most inclusive threshold for sanity check
gw_t2d_col <- "Pt_1"
gw_cad_col <- "Pt_1"
if (!gw_t2d_col %in% colnames(gw_t2d)) gw_t2d_col <- t2d_thresh_cols[length(t2d_thresh_cols)]
if (!gw_cad_col %in% colnames(gw_cad)) gw_cad_col <- cad_thresh_cols[length(cad_thresh_cols)]

setnames(gw_t2d, gw_t2d_col, "GW_T2D_PRS")
setnames(gw_cad, gw_cad_col, "GW_CAD_PRS")

dat <- merge(dat, gw_t2d[, .(FID, IID, GW_T2D_PRS)], by = c("FID", "IID"), all.x = TRUE)
dat <- merge(dat, gw_cad[, .(FID, IID, GW_CAD_PRS)], by = c("FID", "IID"), all.x = TRUE)

# Filter to EUR validation
eur_val <- dat[group == "eur_validation"]
cat(sprintf("  EUR validation: %d individuals\n", nrow(eur_val)))

# Map sex: UKB 0=female, 1=male
eur_val[, sex_label := fifelse(sex == 1L, "Men", "Women")]
cat(sprintf("  Men: %d, Women: %d\n",
            sum(eur_val$sex_label == "Men", na.rm = TRUE),
            sum(eur_val$sex_label == "Women", na.rm = TRUE)))

# Unit conversions
eur_val[, LDL_mgdl := LDL_corrected * 38.67]
eur_val[, HbA1c_pct := 0.09148 * HbA1c_corrected + 2.152]

# Cluster labels from config
if (!is.null(cfg$cluster_labels)) {
  cluster_labels <- unlist(cfg$cluster_labels)
} else {
  cluster_labels <- setNames(gsub("PRS_", "", prs_cols), gsub("PRS_", "", prs_cols))
}

# ============================================================
# Define metrics
# ============================================================
metrics <- list(
  list(col = "LDL_mgdl",      label = "LDL (mg/dL)",     color = "orange",  threshold = 70, ylim_min = NA),
  list(col = "HbA1c_pct",     label = "HbA1c (%)",       color = "#2ca02c", threshold = 6.5, ylim_min = NA),
  list(col = "PCE_ASCVD",     label = "PCE Risk (%)",    color = "#800020", threshold = 5.0, ylim_min = 0),
  list(col = "PREVENT_ASCVD", label = "PREVENT Risk (%)", color = "#E74C3C", threshold = 5.0, ylim_min = 0)
)

# data.table ntile
dt_ntile <- function(x, n = 10L) {
  ranks <- frank(x, ties.method = "random", na.last = "keep")
  as.integer(ceiling(ranks / (sum(!is.na(x)) / n)))
}

# ============================================================
# Helper: generate 2x2 plot for a given PRS column
# ============================================================
make_2x2_plot <- function(data, prs_col_name, title_label) {
  data[, decile := dt_ntile(get(prs_col_name), 10L)]

  panel_plots <- list()

  for (m in metrics) {
    agg <- data[!is.na(get(m$col)) & !is.na(decile) & !is.na(sex_label),
                .(mean_val = mean(get(m$col), na.rm = TRUE)),
                by = .(decile, sex_label)]

    p <- ggplot(agg, aes(x = decile, y = mean_val,
                          linetype = sex_label, group = sex_label)) +
      geom_hline(yintercept = m$threshold, linetype = "dashed",
                 color = "grey50", linewidth = 0.5) +
      geom_point(size = 2.5, color = m$color) +
      geom_line(linewidth = 0.8, color = m$color) +
      scale_linetype_manual(values = c("Men" = "solid", "Women" = "dashed"),
                            name = "Sex") +
      scale_x_continuous(breaks = 1:10) +
      labs(x = "PGS Decile", y = m$label) +
      theme_minimal(base_size = 13) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )

    if (!is.na(m$ylim_min)) {
      y_max <- max(agg$mean_val, m$threshold, na.rm = TRUE) * 1.05
      p <- p + scale_y_continuous(limits = c(m$ylim_min, y_max))
    }

    panel_plots[[length(panel_plots) + 1]] <- p
  }

  combined <- plot_grid(
    panel_plots[[1]], panel_plots[[2]],
    panel_plots[[3]], panel_plots[[4]],
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 14
  )

  titled <- ggdraw() +
    draw_label(title_label, fontface = "bold", size = 16, x = 0.5, y = 0.5) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  final <- plot_grid(titled, combined, ncol = 1, rel_heights = c(0.05, 1))

  data[, decile := NULL]
  final
}

# ============================================================
# Generate plots per cluster
# ============================================================
cat("\n--- Generating cluster PRS plots ---\n")

for (prs_col in prs_cols) {
  k_id <- gsub("PRS_", "", prs_col)
  cluster_name <- cluster_labels[[k_id]]
  safe_name <- gsub("[^A-Za-z0-9]", "", cluster_name)

  cat(sprintf("  Cluster %s (%s)...\n", k_id, cluster_name))

  final <- make_2x2_plot(eur_val, prs_col, cluster_name)

  out_path <- file.path(results_dir,
                         sprintf("decile_biomarker_%s_%s.png", k_id, safe_name))
  ggsave(out_path, final, width = 12, height = 10, dpi = 300)
  cat(sprintf("    Saved: %s\n", basename(out_path)))
}

# ============================================================
# Sanity check: genome-wide PRS plots
# ============================================================
cat("\n--- Generating genome-wide PRS sanity check plots ---\n")

gw_prs <- list(
  list(col = "GW_T2D_PRS", name = "Genome-Wide T2D PRS", file = "decile_biomarker_GW_T2D.png"),
  list(col = "GW_CAD_PRS", name = "Genome-Wide CAD PRS", file = "decile_biomarker_GW_CAD.png")
)

for (gw in gw_prs) {
  cat(sprintf("  %s...\n", gw$name))

  final <- make_2x2_plot(eur_val, gw$col, gw$name)

  out_path <- file.path(results_dir, gw$file)
  ggsave(out_path, final, width = 12, height = 10, dpi = 300)
  cat(sprintf("    Saved: %s\n", gw$file))
}

# ============================================================
# Combined cluster PGS: T2D+CAD joint model weights
# ============================================================
cat("\n--- Generating combined cluster PGS plot ---\n")

pheno <- fread("results/a0_analysis/prs_ct/phenotypes_combined.txt")
covar <- fread(cfg$covariate_file)

# Merge phenotypes and covariates; drop overlapping columns already in eur_val
eur_val_model <- merge(eur_val, pheno[, .(FID, IID, T2D, CAD)], by = c("FID", "IID"), all.x = TRUE)
covar_new_cols <- setdiff(colnames(covar), colnames(eur_val_model))
eur_val_model <- merge(eur_val_model, covar[, c("FID", "IID", covar_new_cols), with = FALSE],
                       by = c("FID", "IID"), all.x = TRUE)

eur_val_model[, Synthetic := fifelse(T2D == 1 & CAD == 1, 1L,
                             fifelse(T2D == 0 | CAD == 0, 0L, NA_integer_))]

covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))
model_data <- eur_val_model[, c("Synthetic", prs_cols, covar_terms), with = FALSE]
model_data <- na.omit(model_data)

if ("Batch" %in% colnames(model_data)) {
  model_data[, Batch := as.factor(Batch)]
}

# Standardize PRS for model fitting
prs_means <- sapply(prs_cols, function(pc) mean(model_data[[pc]]))
prs_sds   <- sapply(prs_cols, function(pc) sd(model_data[[pc]]))
for (pc in prs_cols) {
  model_data[[pc]] <- scale(model_data[[pc]])[, 1]
}

rhs <- paste(c(prs_cols, covar_terms), collapse = " + ")
fit <- glm(as.formula(paste("Synthetic ~", rhs)), data = model_data, family = binomial)

prs_betas <- coef(fit)[prs_cols]
cat(sprintf("  Joint model betas (Synthetic):\n"))
for (i in seq_along(prs_cols)) {
  cat(sprintf("    %s: %.4f\n", prs_cols[i], prs_betas[i]))
}

# Compute composite PGS for ALL eur_val individuals using the model betas
# Standardize using the same means/sds from the model fitting data
prs_std <- as.matrix(eur_val[, ..prs_cols])
for (i in seq_along(prs_cols)) {
  prs_std[, i] <- (prs_std[, i] - prs_means[i]) / prs_sds[i]
}
eur_val[, Combined_PGS := as.numeric(prs_std %*% prs_betas)]

cat(sprintf("  Combined PGS: mean=%.4f, SD=%.4f\n",
            mean(eur_val$Combined_PGS, na.rm = TRUE),
            sd(eur_val$Combined_PGS, na.rm = TRUE)))

final <- make_2x2_plot(eur_val, "Combined_PGS", "Combined Cluster PGS (T2D+CAD Model)")

out_path <- file.path(results_dir, "decile_biomarker_Combined_Cluster_PGS.png")
ggsave(out_path, final, width = 12, height = 10, dpi = 300)
cat(sprintf("  Saved: %s\n", basename(out_path)))

# ============================================================
# Forest plot: top/bottom 20% PCE & HbA1c per cluster
# ============================================================
cat("\n--- Generating top/bottom 20% forest plot ---\n")

forest_list <- list()
for (prs_col in prs_cols) {
  k_id <- gsub("PRS_", "", prs_col)
  eur_val[, tmp_decile := dt_ntile(get(prs_col), 10L)]

  for (bio in list(list(col = "PCE_ASCVD", label = "PCE Risk (%)"),
                   list(col = "HbA1c_pct", label = "HbA1c (%)"))) {
    for (grp in list(list(decs = c(9L, 10L), label = "Top 20%"),
                     list(decs = c(1L, 2L),  label = "Bottom 20%"))) {
      subset <- eur_val[tmp_decile %in% grp$decs & !is.na(get(bio$col))]
      m <- mean(subset[[bio$col]], na.rm = TRUE)
      s <- sd(subset[[bio$col]], na.rm = TRUE)
      forest_list[[length(forest_list) + 1]] <- data.table(
        cluster = k_id,
        cluster_label = cluster_labels[[k_id]],
        biomarker = bio$label,
        group = grp$label,
        mean_val = m,
        sd_lower = m - s,
        sd_upper = m + s
      )
    }
  }
  eur_val[, tmp_decile := NULL]
}

forest_dt <- rbindlist(forest_list)
forest_dt[, color_key := paste0(biomarker, " — ", group)]
forest_dt[, cluster_label := factor(cluster_label,
  levels = rev(cluster_labels[paste0("K", 1:n_clusters)]))]

color_map <- c(
  "PCE Risk (%) — Top 20%"    = "#D55E00",
  "PCE Risk (%) — Bottom 20%" = "#F5A07A",
  "HbA1c (%) — Top 20%"       = "#0072B2",
  "HbA1c (%) — Bottom 20%"    = "#8ECAE6"
)

p_forest <- ggplot(forest_dt, aes(x = mean_val, y = cluster_label, color = color_key)) +
  geom_pointrange(aes(xmin = sd_lower, xmax = sd_upper),
                  position = position_dodge(width = 0.6),
                  size = 0.5, linewidth = 0.7) +
  facet_wrap(~ biomarker, scales = "free_x") +
  scale_color_manual(values = color_map, name = NULL) +
  labs(x = "Mean ± SD", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

forest_path <- file.path(results_dir, "forest_top_bottom_decile.png")
ggsave(forest_path, p_forest, width = 12, height = 6, dpi = 300)
cat(sprintf("  Saved: %s\n", basename(forest_path)))

cat("\n=== Done ===\n")
