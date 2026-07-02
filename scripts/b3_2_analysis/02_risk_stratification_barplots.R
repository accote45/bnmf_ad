#!/usr/bin/env Rscript
# 02_risk_stratification_barplots.R
# Grouped barplots: within each clinical risk category (Low/Borderline/Intermediate/High),
# show mean predicted ASCVD risk stratified by cluster PGS tertile (Low/Intermediate/High).
#
# Usage:
#   Rscript scripts/b3_2_analysis/02_risk_stratification_barplots.R --config config/b3_2_config.yaml --hard-assignment

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

cat("=== B3.2: Risk Stratification Barplots ===\n")

# ============================================================
# Load data
# ============================================================
cat("\n--- Loading data ---\n")

prs_all <- fread(cfg$prs_file)
biomarker <- fread(cfg$biomarker_file,
                   select = c("FID", "IID", "PCE_ASCVD", "PREVENT_ASCVD"))

dat <- merge(prs_all, biomarker, by = c("FID", "IID"), all.x = TRUE)
eur_val <- dat[group == "eur_validation"]
cat(sprintf("  EUR validation: %d individuals\n", nrow(eur_val)))

# Cluster labels
cluster_labels <- unlist(cfg$cluster_labels)

# Selected clusters
selected <- c("K3", "K4", "K7", "K10")
selected_labels <- cluster_labels[selected]
cat(sprintf("  Selected clusters: %s\n",
            paste(sprintf("%s (%s)", selected, selected_labels), collapse = ", ")))

# ============================================================
# Risk categories and PGS strata
# ============================================================
risk_breaks <- c(0, 5, 7.5, 20, Inf)
risk_labels <- c("Low\n(<5%)", "Borderline\n(5-7.4%)", "Intermediate\n(7.5-20%)", "High\n(>20%)")

pgs_colors <- c("Low PGS" = "#FFD700", "Intermediate PGS" = "#FF8C00", "High PGS" = "#CC0000")

# ============================================================
# Helper: make barplot for one risk score
# ============================================================
make_risk_barplots <- function(data, risk_col, risk_label, selected_clusters) {
  panel_plots <- list()

  for (k_id in selected_clusters) {
    prs_col <- paste0("PRS_", k_id)
    cluster_name <- cluster_labels[[k_id]]

    sub <- data[!is.na(get(risk_col)) & !is.na(get(prs_col))]

    # Risk categories
    sub[, risk_cat := cut(get(risk_col), breaks = risk_breaks,
                          labels = risk_labels, right = FALSE)]

    # PGS strata: bottom 20%, middle 60%, top 20%
    q20 <- quantile(sub[[prs_col]], 0.20, na.rm = TRUE)
    q80 <- quantile(sub[[prs_col]], 0.80, na.rm = TRUE)
    sub[, pgs_cat := fifelse(get(prs_col) <= q20, "Low PGS",
                    fifelse(get(prs_col) >= q80, "High PGS",
                            "Intermediate PGS"))]
    sub[, pgs_cat := factor(pgs_cat, levels = c("Low PGS", "Intermediate PGS", "High PGS"))]

    # Aggregate: mean ± SE
    agg <- sub[!is.na(risk_cat),
               .(mean_risk = mean(get(risk_col), na.rm = TRUE),
                 se_risk = sd(get(risk_col), na.rm = TRUE) / sqrt(.N),
                 n = .N),
               by = .(risk_cat, pgs_cat)]
    agg[, ci_lower := mean_risk - 1.96 * se_risk]
    agg[, ci_upper := mean_risk + 1.96 * se_risk]

    p <- ggplot(agg, aes(x = risk_cat, y = mean_risk, fill = pgs_cat)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                    position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.5) +
      scale_fill_manual(values = pgs_colors, name = "Cluster PGS") +
      labs(x = paste(risk_label, "Risk Category"),
           y = paste("Mean", risk_label, "Risk (%)"),
           title = cluster_name) +
      theme_minimal(base_size = 13) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )

    panel_plots[[k_id]] <- p
  }

  combined <- plot_grid(
    panel_plots[[1]], panel_plots[[2]],
    panel_plots[[3]], panel_plots[[4]],
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 14
  )

  combined
}

# ============================================================
# Generate PCE figure
# ============================================================
cat("\n--- Generating PCE risk stratification figure ---\n")

pce_plot <- make_risk_barplots(eur_val, "PCE_ASCVD", "PCE", selected)

out_pce <- file.path(results_dir, "risk_stratification_PCE.png")
ggsave(out_pce, pce_plot, width = 14, height = 10, dpi = 300)
cat(sprintf("  Saved: %s\n", out_pce))

# ============================================================
# Generate PREVENT figure
# ============================================================
cat("\n--- Generating PREVENT risk stratification figure ---\n")

prevent_plot <- make_risk_barplots(eur_val, "PREVENT_ASCVD", "PREVENT", selected)

out_prevent <- file.path(results_dir, "risk_stratification_PREVENT.png")
ggsave(out_prevent, prevent_plot, width = 14, height = 10, dpi = 300)
cat(sprintf("  Saved: %s\n", out_prevent))

# ============================================================
# Combined cluster PGS risk stratification
# ============================================================
cat("\n--- Computing combined cluster PGS (T2D+CAD joint model) ---\n")

prs_cols <- grep("^PRS_K", colnames(eur_val), value = TRUE)

pheno <- fread("results/a0_analysis/prs_ct/phenotypes_combined.txt")
covar <- fread(cfg$covariate_file)

eur_val_model <- merge(eur_val, pheno[, .(FID, IID, T2D, CAD)], by = c("FID", "IID"), all.x = TRUE)
covar_new_cols <- setdiff(colnames(covar), colnames(eur_val_model))
eur_val_model <- merge(eur_val_model, covar[, c("FID", "IID", covar_new_cols), with = FALSE],
                       by = c("FID", "IID"), all.x = TRUE)

eur_val_model[, Synthetic := fifelse(T2D == 1 & CAD == 1, 1L,
                             fifelse(T2D == 0 | CAD == 0, 0L, NA_integer_))]

covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))
model_data <- eur_val_model[, c("Synthetic", prs_cols, covar_terms), with = FALSE]
model_data <- na.omit(model_data)
if ("Batch" %in% colnames(model_data)) model_data[, Batch := as.factor(Batch)]

prs_means <- sapply(prs_cols, function(pc) mean(model_data[[pc]]))
prs_sds   <- sapply(prs_cols, function(pc) sd(model_data[[pc]]))
for (pc in prs_cols) model_data[[pc]] <- scale(model_data[[pc]])[, 1]

rhs <- paste(c(prs_cols, covar_terms), collapse = " + ")
fit <- glm(as.formula(paste("Synthetic ~", rhs)), data = model_data, family = binomial)
prs_betas <- coef(fit)[prs_cols]

prs_std <- as.matrix(eur_val[, ..prs_cols])
for (i in seq_along(prs_cols)) prs_std[, i] <- (prs_std[, i] - prs_means[i]) / prs_sds[i]
eur_val[, Combined_PGS := as.numeric(prs_std %*% prs_betas)]

cat(sprintf("  Combined PGS: mean=%.4f, SD=%.4f\n",
            mean(eur_val$Combined_PGS, na.rm = TRUE),
            sd(eur_val$Combined_PGS, na.rm = TRUE)))

# PGS strata
q20 <- quantile(eur_val$Combined_PGS, 0.20, na.rm = TRUE)
q80 <- quantile(eur_val$Combined_PGS, 0.80, na.rm = TRUE)
eur_val[, pgs_cat := fifelse(Combined_PGS <= q20, "Low PGS",
                    fifelse(Combined_PGS >= q80, "High PGS",
                            "Intermediate PGS"))]
eur_val[, pgs_cat := factor(pgs_cat, levels = c("Low PGS", "Intermediate PGS", "High PGS"))]

# Helper: single-panel barplot for one risk score
make_combined_barplot <- function(data, risk_col, risk_label) {
  sub <- data[!is.na(get(risk_col)) & !is.na(pgs_cat)]
  sub[, risk_cat := cut(get(risk_col), breaks = risk_breaks,
                        labels = risk_labels, right = FALSE)]

  agg <- sub[!is.na(risk_cat),
             .(mean_risk = mean(get(risk_col), na.rm = TRUE),
               se_risk = sd(get(risk_col), na.rm = TRUE) / sqrt(.N),
               n = .N),
             by = .(risk_cat, pgs_cat)]
  agg[, ci_lower := mean_risk - 1.96 * se_risk]
  agg[, ci_upper := mean_risk + 1.96 * se_risk]

  ggplot(agg, aes(x = risk_cat, y = mean_risk, fill = pgs_cat)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.5) +
    scale_fill_manual(values = pgs_colors, name = "Combined Cluster PGS") +
    labs(x = paste(risk_label, "Risk Category"),
         y = paste("Mean", risk_label, "Risk (%)"),
         title = paste("Combined Cluster PGS:", risk_label, "Risk Stratification")) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

cat("\n--- Generating combined PGS risk stratification figures ---\n")

pce_combined <- make_combined_barplot(eur_val, "PCE_ASCVD", "PCE")
out_pce_c <- file.path(results_dir, "risk_stratification_PCE_combined.png")
ggsave(out_pce_c, pce_combined, width = 8, height = 6, dpi = 300)
cat(sprintf("  Saved: %s\n", out_pce_c))

prevent_combined <- make_combined_barplot(eur_val, "PREVENT_ASCVD", "PREVENT")
out_prevent_c <- file.path(results_dir, "risk_stratification_PREVENT_combined.png")
ggsave(out_prevent_c, prevent_combined, width = 8, height = 6, dpi = 300)
cat(sprintf("  Saved: %s\n", out_prevent_c))

cat("\n=== Done ===\n")
