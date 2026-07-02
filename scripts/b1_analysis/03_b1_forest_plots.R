#!/usr/bin/env Rscript
# 03_b1_forest_plots.R
# Publication-quality forest plots comparing bNMF cluster PRS vs genome-wide PRS.
# EUR validation + cross-ancestry (AFR/EAS/SAS) portability.
#
# Usage: Rscript scripts/b1_analysis/03_b1_forest_plots.R

library(data.table)
library(ggplot2)
library(yaml)
library(pROC)

args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/b1_config.yaml"
if ("--config" %in% args) {
  config_path <- args[which(args == "--config") + 1]
}
hard_assignment <- "--hard-assignment" %in% args

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
if (hard_assignment) {
  results_dir <- file.path(results_dir, "prs_hard_assignment")
  cat("  ** HARD ASSIGNMENT MODE **\n")
  cat(sprintf("  Output dir: %s\n", results_dir))
}

cat("=== B1 Step 3: Forest Plots ===\n")

# --- Nagelkerke R-squared ---
nagelkerke_r2 <- function(full_model, null_model, n) {
  ll_full <- as.numeric(logLik(full_model))
  ll_null <- as.numeric(logLik(null_model))
  cox_snell <- 1 - exp((2 / n) * (ll_null - ll_full))
  max_r2 <- 1 - exp((2 / n) * ll_null)
  cox_snell / max_r2
}

# ============================================================
# Step 1: Load data
# ============================================================
cat("\n--- Loading data ---\n")

# Cluster PRS association results (already computed)
cluster_results <- fread(file.path(results_dir, "association_results_all.csv"))
cluster_results <- cluster_results[model_type == "individual"]
cat(sprintf("  Cluster results: %d rows\n", nrow(cluster_results)))

# Genome-wide PRS (bNMF-matched; built by 01b_compute_genomewide_prs.R).
# Columns: GW_T2D, GW_CAD (trait-specific) + GW_T2D_combined, GW_CAD_combined
# (supplementary, shared variant universe). Carries the sample `group`.
gw_all <- fread(file.path(results_dir, "prs/genomewide_prs_all.tsv"))
cat(sprintf("  GW PRS: %d individuals, %d groups\n",
            nrow(gw_all), length(unique(gw_all$group))))

# Phenotypes and covariates
pheno <- fread(cfg$phenotypes$phenotype_file)
covar <- fread(cfg$phenotypes$covariate_file)

# Merge GW PRS (carries group) with phenotypes + covariates
dat <- merge(gw_all, pheno, by = c("FID", "IID"))
dat <- merge(dat, covar, by = c("FID", "IID"))

# Create synthetic outcome
dat[, Synthetic := fifelse(T2D == 1 & CAD == 1, 1L,
                   fifelse(T2D == 0 | CAD == 0, 0L, NA_integer_))]

cat(sprintf("  Merged dataset: %d individuals\n", nrow(dat)))

covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))

# ============================================================
# Step 2 & 3: Genome-wide PRS regressions (fixed 5e-8 — no threshold scan)
# ============================================================
# The genome-wide PRS now uses the bNMF variant pipeline (single 5e-8 set), so
# there is no p-value-threshold sweep. Main comparator = trait-specific GW PRS;
# the two combined-variant scores are run too and flagged for supplementary use.
cat("\n--- Running genome-wide PRS regressions ---\n")

gw_prs_configs <- list(
  list(name = "GW T2D",                     col = "GW_T2D",          supp = FALSE),
  list(name = "GW CAD",                     col = "GW_CAD",          supp = FALSE),
  list(name = "GW T2D (combined variants)", col = "GW_T2D_combined", supp = TRUE),
  list(name = "GW CAD (combined variants)", col = "GW_CAD_combined", supp = TRUE)
)

outcomes <- c("T2D", "CAD", "Synthetic")
val_groups <- c("eur_validation", "afr_validation", "eas_validation", "sas_validation")

gw_results_list <- list()

for (gw in gw_prs_configs) {
  for (grp in val_groups) {
    grp_data <- dat[group == grp]

    for (outcome in outcomes) {
      prs_col <- gw$col
      model_data <- grp_data[, c(outcome, prs_col, covar_terms), with = FALSE]
      model_data <- na.omit(model_data)

      n_cases <- sum(model_data[[outcome]] == 1)
      n_controls <- sum(model_data[[outcome]] == 0)

      if (n_cases < 10 || n_controls < 10) {
        cat(sprintf("  SKIP %s ~ %s in %s: cases=%d, controls=%d\n",
                    outcome, gw$name, grp, n_cases, n_controls))
        next
      }

      model_data[[prs_col]] <- scale(model_data[[prs_col]])[, 1]

      rhs_full <- paste(c(prs_col, covar_terms), collapse = " + ")
      rhs_null <- paste(covar_terms, collapse = " + ")

      fit_full <- tryCatch(
        glm(as.formula(paste(outcome, "~", rhs_full)),
            data = model_data, family = binomial),
        error = function(e) NULL
      )
      if (is.null(fit_full)) next

      fit_null <- glm(as.formula(paste(outcome, "~", rhs_null)),
                      data = model_data, family = binomial)

      coef_summ <- summary(fit_full)$coefficients
      prs_row <- coef_summ[prs_col, ]

      or <- exp(prs_row["Estimate"])
      ci_lower <- exp(prs_row["Estimate"] - 1.96 * prs_row["Std. Error"])
      ci_upper <- exp(prs_row["Estimate"] + 1.96 * prs_row["Std. Error"])
      p_val <- prs_row["Pr(>|z|)"]
      r2 <- nagelkerke_r2(fit_full, fit_null, nrow(model_data))

      cat(sprintf("  %s ~ %s [%s]: OR=%.3f [%.3f-%.3f], P=%.2e\n",
                  outcome, gw$name, grp, or, ci_lower, ci_upper, p_val))

      gw_results_list[[length(gw_results_list) + 1]] <- data.table(
        group = grp, outcome = outcome, model_type = "individual",
        predictor = gw$name, is_supplementary = isTRUE(gw$supp),
        OR = round(or, 4), CI_lower = round(ci_lower, 4),
        CI_upper = round(ci_upper, 4),
        p_value = p_val, nagelkerke_r2 = round(r2, 6),
        n_cases = n_cases, n_controls = n_controls,
        n_total = n_cases + n_controls
      )
    }
  }
}

gw_results <- rbindlist(gw_results_list)
fwrite(gw_results, file.path(results_dir, "genome_wide_prs_results.csv"))
cat(sprintf("\n  Genome-wide PRS results: %d rows\n", nrow(gw_results)))

# ============================================================
# Step 4: Combine cluster + genome-wide results for plotting
# ============================================================
cat("\n--- Preparing forest plot data ---\n")

# Rename cluster predictors for display
cluster_plot <- cluster_results[group %in% val_groups,
  .(group, outcome, predictor, OR, CI_lower, CI_upper, p_value)]
cluster_plot[, predictor := gsub("PRS_", "", predictor)]

# Main figures use only the trait-specific GW PRS; the combined-variant scores
# are reserved for the supplementary forests built at the end of this script.
gw_plot <- gw_results[is_supplementary == FALSE,
  .(group, outcome, predictor, OR, CI_lower, CI_upper, p_value)]

plot_data <- rbindlist(list(cluster_plot, gw_plot))

# Predictor ordering: GW T2D (leftmost), clusters in desired display order, GW CAD (rightmost)
k_preds <- sort(unique(grep("^K\\d+$", plot_data$predictor, value = TRUE)))
n_clusters <- length(k_preds)

# Semantic cluster labels from config (falls back to generic K1..Kn)
if (!is.null(cfg$cluster_labels)) {
  cluster_labels <- unlist(cfg$cluster_labels)
} else {
  cluster_labels <- setNames(k_preds, k_preds)
}

# Desired META cluster display order: Glycemic, Obesity, SHBG, Adiponectin,
# Triglycerides-HDL, ALP-LDL, Metabolic, Platelet, Blood Pressure-Stature, Lpa
desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
k_ordered <- c(desired_k_order[desired_k_order %in% k_preds],
               sort(setdiff(k_preds, desired_k_order)))

# Remap K predictors to semantic names
for (k in names(cluster_labels)) {
  plot_data[predictor == k, predictor := cluster_labels[[k]]]
}
# GW T2D forced leftmost, clusters in desired order, GW CAD forced rightmost
pred_levels <- c("GW T2D", unname(cluster_labels[k_ordered]), "GW CAD")
plot_data[, predictor := factor(predictor, levels = pred_levels)]

# Grey shading behind the two genome-wide predictors (now leftmost + rightmost)
gw_band <- data.frame(xmin = c(0.5, n_clusters + 1.5),
                      xmax = c(1.5, n_clusters + 2.5),
                      ymin = -Inf, ymax = Inf)

# Outcome labels
outcome_labels <- c("T2D" = "Type 2 Diabetes", "CAD" = "Coronary Artery Disease",
                     "Synthetic" = "T2D and CAD")
plot_data[, outcome_label := outcome_labels[outcome]]
plot_data[, outcome_label := factor(outcome_label, levels = outcome_labels)]

# Ancestry labels for non-EUR
ancestry_labels <- c("afr_validation" = "AFR", "eas_validation" = "EAS",
                      "sas_validation" = "SAS", "eur_validation" = "EUR")
plot_data[, ancestry := ancestry_labels[group]]

# Significance annotation
plot_data[, sig := fifelse(p_value < 0.001, "***",
                  fifelse(p_value < 0.01, "**",
                  fifelse(p_value < 0.05, "*", "")))]

# ============================================================
# Step 5: EUR validation forest plot
# ============================================================
cat("\n--- Creating EUR validation forest plot ---\n")

outcome_colors <- c("Type 2 Diabetes" = "#E74C3C",
                     "Coronary Artery Disease" = "#3498DB",
                     "T2D and CAD" = "#9B59B6")

eur_data <- plot_data[group == "eur_validation"]

p_eur <- ggplot(eur_data, aes(x = predictor, y = OR, color = outcome_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.6), width = 0.25, linewidth = 0.7) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  geom_rect(data = gw_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, alpha = 0.08, fill = "grey50") +
  labs(
    x = "PRS Predictor",
    y = "OR/SD",
    title = "EUR Validation: Cluster PRS vs Genome-Wide PRS"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "forest_eur_validation.png"),
       p_eur, width = 12, height = 7, dpi = 300)
cat("  Saved: forest_eur_validation.png\n")

# ============================================================
# Step 6: Non-EUR combined forest plot
# ============================================================
cat("\n--- Creating non-EUR combined forest plot ---\n")

noneur_data <- plot_data[group %in% c("afr_validation", "eas_validation", "sas_validation")]
noneur_data[, ancestry := factor(ancestry, levels = c("AFR", "EAS", "SAS"))]

p_noneur <- ggplot(noneur_data, aes(x = predictor, y = OR, color = outcome_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.6), width = 0.25, linewidth = 0.6) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  geom_rect(data = gw_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, alpha = 0.08, fill = "grey50") +
  facet_wrap(~ancestry, ncol = 1, scales = "free_y") +
  labs(
    x = "PRS Predictor",
    y = "OR/SD",
    title = "Cross-Ancestry Portability: EUR-Derived Cluster PRS vs Genome-Wide PRS"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "forest_noneur_validation.png"),
       p_noneur, width = 12, height = 12, dpi = 300)
cat("  Saved: forest_noneur_validation.png\n")

# ============================================================
# Step 6b: Combined validation forest plot (EUR + non-EUR, one figure)
# ============================================================
# Same data/geoms/theme as the EUR and non-EUR plots above; all four validation
# ancestries are stacked as facets (EUR first), each with a free y-axis. Uses the
# faceted (non-EUR) point/linewidth styling for consistency across panels.
cat("\n--- Creating combined (EUR + non-EUR) validation forest plot ---\n")

combined_data <- plot_data[group %in% c("eur_validation", "afr_validation",
                                        "eas_validation", "sas_validation")]
combined_data[, ancestry := factor(ancestry, levels = c("EUR", "AFR", "EAS", "SAS"))]

p_combined <- ggplot(combined_data, aes(x = predictor, y = OR, color = outcome_label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.6), width = 0.25, linewidth = 0.6) +
  scale_color_manual(values = outcome_colors, name = "Outcome") +
  geom_rect(data = gw_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, alpha = 0.08, fill = "grey50") +
  facet_wrap(~ancestry, ncol = 1, scales = "free_y") +
  labs(
    x = "PRS Predictor",
    y = "OR/SD",
    title = "Cluster PRS vs Genome-Wide PRS across Ancestries"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "forest_combined_validation.png"),
       p_combined, width = 12, height = 16, dpi = 300)
cat("  Saved: forest_combined_validation.png\n")

# ============================================================
# Step 7: Cluster PRS correlation heatmap (EUR only)
# ============================================================
cat("\n--- Creating cluster PRS correlation heatmap (EUR) ---\n")

prs_all <- fread(file.path(results_dir, "prs/cluster_prs_all.tsv"))
eur_prs <- prs_all[group %in% c("eur_train", "eur_validation")]
prs_cols <- grep("^PRS_K", colnames(eur_prs), value = TRUE)

cor_mat <- cor(eur_prs[, ..prs_cols], use = "complete.obs", method = "pearson")
cat("  Pearson correlation matrix:\n")
print(round(cor_mat, 3))

# Load W matrix for cluster size counts
w_mat <- fread(cfg$a1_results$w_matrix)
kcols <- setdiff(colnames(w_mat), "VAR_ID")
n_k <- length(kcols)
# Heatmap axis order follows the desired META cluster display order (Glycemic -> Lpa);
# desired_k_order is defined in the forest-ordering step above.
k_levels <- c(desired_k_order[desired_k_order %in% kcols],
              sort(setdiff(kcols, desired_k_order)))

plot_rows <- list()

if (hard_assignment) {
  # Hard assignment: full Pearson r matrix, cluster sizes on diagonal
  cluster_assignment <- kcols[max.col(as.matrix(w_mat[, ..kcols]), ties.method = "first")]
  cluster_sizes <- table(factor(cluster_assignment, levels = kcols))
  cat("\n  Cluster sizes (hard assignment):\n")
  print(cluster_sizes)

  for (i in seq_along(k_levels)) {
    for (j in seq_along(k_levels)) {
      ki <- k_levels[i]
      kj <- k_levels[j]
      if (i == j) {
        plot_rows[[length(plot_rows) + 1]] <- data.table(
          Var1 = ki, Var2 = kj, r = NA_real_,
          label = sprintf("n=%d", cluster_sizes[[ki]]),
          region = "diagonal"
        )
      } else if (i < j) {
        rval <- cor_mat[paste0("PRS_", ki), paste0("PRS_", kj)]
        plot_rows[[length(plot_rows) + 1]] <- data.table(
          Var1 = ki, Var2 = kj, r = rval,
          label = sprintf("%.2f", rval),
          region = "upper"
        )
      } else {
        plot_rows[[length(plot_rows) + 1]] <- data.table(
          Var1 = ki, Var2 = kj, r = NA_real_,
          label = "",
          region = "lower"
        )
      }
    }
  }
} else {
  # Soft assignment: upper (top-right) triangle + diagonal show the Pearson
  # correlation coefficient; the lower (bottom-left) triangle is left blank.
  for (i in seq_along(k_levels)) {
    for (j in seq_along(k_levels)) {
      ki <- k_levels[i]
      kj <- k_levels[j]
      if (i < j) {
        # bottom-left triangle: blank (no fill, no label)
        plot_rows[[length(plot_rows) + 1]] <- data.table(
          Var1 = ki, Var2 = kj, r = NA_real_, label = "", region = "blank"
        )
      } else {
        rval <- if (i == j) 1 else cor_mat[paste0("PRS_", ki), paste0("PRS_", kj)]
        plot_rows[[length(plot_rows) + 1]] <- data.table(
          Var1 = ki, Var2 = kj, r = rval,
          label = sprintf("%.2f", rval),
          region = "shown"
        )
      }
    }
  }
}

hm_dt <- rbindlist(plot_rows)

# Apply semantic cluster labels to heatmap axes
k_display <- unname(cluster_labels[k_levels])
for (k in names(cluster_labels)) {
  hm_dt[Var1 == k, Var1 := cluster_labels[[k]]]
  hm_dt[Var2 == k, Var2 := cluster_labels[[k]]]
}
hm_dt[, Var1 := factor(Var1, levels = k_display)]
hm_dt[, Var2 := factor(Var2, levels = rev(k_display))]

if (hard_assignment) {
  p_heatmap <- ggplot(hm_dt, aes(x = Var1, y = Var2)) +
    geom_tile(data = hm_dt[region == "upper"],
              aes(fill = r), color = "white", linewidth = 0.5) +
    geom_tile(data = hm_dt[region == "diagonal"],
              fill = "grey90", color = "white", linewidth = 0.5) +
    geom_tile(data = hm_dt[region == "lower"],
              fill = NA, color = NA) +
    geom_text(data = hm_dt[region != "lower"],
              aes(label = label), size = 3.8, color = "black") +
    scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r",
                         na.value = "grey90") +
    labs(x = NULL, y = NULL,
         title = "Hard-Assignment Cluster PRS: Pearson Correlation") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed()
} else {
  p_heatmap <- ggplot(hm_dt, aes(x = Var1, y = Var2)) +
    geom_tile(data = hm_dt[region != "blank"],
              aes(fill = r), color = "white", linewidth = 0.5) +
    geom_text(data = hm_dt[region != "blank"],
              aes(label = label), size = 3.8, color = "black") +
    scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    labs(x = NULL, y = NULL,
         title = "Cluster PRS: Pearson Correlation") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed()
}

ggsave(file.path(results_dir, "heatmap_cluster_prs_corr.png"),
       p_heatmap, width = 10, height = 8, dpi = 300)
cat("  Saved: heatmap_cluster_prs_corr.png\n")

# ============================================================
# Step 8 (supplementary): forests including the combined-variant GW PRS
# ============================================================
cat("\n--- Creating supplementary forests (incl. combined-variant GW PRS) ---\n")

# Cluster PRS + ALL four GW predictors (trait-specific + combined-variant).
gw_plot_all <- gw_results[, .(group, outcome, predictor, OR, CI_lower, CI_upper, p_value)]
plot_supp <- rbindlist(list(cluster_plot, gw_plot_all))
for (k in names(cluster_labels)) plot_supp[predictor == k, predictor := cluster_labels[[k]]]

gw_levels_supp <- c("GW T2D", "GW CAD",
                    "GW T2D (combined variants)", "GW CAD (combined variants)")
pred_levels_supp <- c(unname(cluster_labels[k_ordered]), gw_levels_supp)
plot_supp[, predictor := factor(predictor, levels = pred_levels_supp)]
plot_supp <- plot_supp[!is.na(predictor)]

plot_supp[, outcome_label := factor(outcome_labels[outcome], levels = outcome_labels)]
plot_supp[, ancestry := ancestry_labels[group]]

# Grey band behind the four GW predictors (contiguous block on the right).
gw_band_supp <- data.frame(xmin = n_clusters + 0.5, xmax = n_clusters + 4.5,
                           ymin = -Inf, ymax = Inf)

supp_forest <- function(df, facet, title, fname, w, h) {
  p <- ggplot(df, aes(x = predictor, y = OR, color = outcome_label)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_rect(data = gw_band_supp,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, alpha = 0.08, fill = "grey50") +
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  position = position_dodge(width = 0.6), width = 0.25, linewidth = 0.6) +
    scale_color_manual(values = outcome_colors, name = "Outcome") +
    labs(x = "PRS Predictor", y = "OR/SD", title = title) +
    theme_minimal(base_size = 14) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 11),
          legend.position = "bottom", legend.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
  if (facet) {
    p <- p + facet_wrap(~ancestry, ncol = 1, scales = "free_y") +
      theme(strip.text = element_text(face = "bold", size = 14))
  }
  ggsave(file.path(results_dir, fname), p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s\n", fname))
}

supp_forest(plot_supp[group == "eur_validation"], FALSE,
  "EUR Validation (supp.): Cluster vs Trait-specific & Combined-variant GW PRS",
  "forest_eur_validation_supp.png", 14, 7)

supp_noneur <- plot_supp[group %in% c("afr_validation", "eas_validation", "sas_validation")]
supp_noneur[, ancestry := factor(ancestry, levels = c("AFR", "EAS", "SAS"))]
supp_forest(supp_noneur, TRUE,
  "Cross-Ancestry (supp.): Cluster vs Trait-specific & Combined-variant GW PRS",
  "forest_noneur_validation_supp.png", 14, 12)

# ============================================================
# Summary
# ============================================================
cat("\n=== Summary ===\n")
cat("\n  EUR validation OR per SD:\n")
print(eur_data[, .(predictor, outcome, OR, CI_lower, CI_upper, p_value, sig)])
cat("\n  Non-EUR validation OR per SD:\n")
print(noneur_data[, .(ancestry, predictor, outcome, OR, CI_lower, CI_upper, p_value, sig)])

cat("\n=== Done ===\n")
