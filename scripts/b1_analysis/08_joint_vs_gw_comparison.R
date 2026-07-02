#!/usr/bin/env Rscript
# 08_joint_vs_gw_comparison.R
# Compare joint cluster PRS model vs genome-wide PRS for T2D/CAD/Synthetic.
# Models: Joint Cluster PRS, GW T2D PRS, GW CAD PRS, Combined (cluster + GW).
# Outputs grouped bar chart (R² and AUC) and CSV summary table.
#
# Usage:
#   Rscript scripts/b1_analysis/08_joint_vs_gw_comparison.R --config config/b1_config.yaml --hard-assignment

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
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== B1 Step 8: Joint Cluster PRS vs Genome-Wide PRS Comparison ===\n")
cat(sprintf("  Config: %s\n", config_path))
cat(sprintf("  Output: %s\n", results_dir))
if (hard_assignment) cat("  ** HARD ASSIGNMENT MODE **\n")

nagelkerke_r2 <- function(full_model, null_model, n) {
  ll_full <- as.numeric(logLik(full_model))
  ll_null <- as.numeric(logLik(null_model))
  cox_snell <- 1 - exp((2 / n) * (ll_null - ll_full))
  max_r2 <- 1 - exp((2 / n) * ll_null)
  cox_snell / max_r2
}

# ============================================================
# Load data
# ============================================================
cat("\n--- Loading data ---\n")

prs_all <- fread(file.path(results_dir, "prs/cluster_prs_all.tsv"))
pheno   <- fread(cfg$phenotypes$phenotype_file)
covar   <- fread(cfg$phenotypes$covariate_file)

prs_cols <- grep("^PRS_K", colnames(prs_all), value = TRUE)
cat(sprintf("  Cluster PRS: %d individuals, %d clusters\n", nrow(prs_all), length(prs_cols)))

merged <- merge(prs_all, pheno, by = c("FID", "IID"), all.x = TRUE)
merged <- merge(merged, covar, by = c("FID", "IID"), all.x = TRUE)

# Load genome-wide PRS
gw_t2d <- fread("results/a0_analysis/prs_ct/prsice_output/T2D.all_score")
gw_cad <- fread("results/a0_analysis/prs_ct/prsice_output/CAD.all_score")

t2d_thresh_cols <- grep("^Pt_", colnames(gw_t2d), value = TRUE)
cad_thresh_cols <- grep("^Pt_", colnames(gw_cad), value = TRUE)

# Use the best-fit threshold selected on eur_train by 08a_select_gw_threshold.R
# (Nagelkerke R2). Falls back to Pt_1 if the selection file is absent or a
# trait/threshold is missing. Selection CSV lives in the base results_dir.
select_best_threshold <- function(trait_name, avail_cols) {
  sel_path <- file.path(cfg$results_dir, "gw_threshold_selection.csv")
  fallback <- if ("Pt_1" %in% avail_cols) "Pt_1" else avail_cols[length(avail_cols)]
  if (!file.exists(sel_path)) {
    cat(sprintf("  [warn] %s not found; using fallback threshold for %s\n", sel_path, trait_name))
    return(fallback)
  }
  sel <- fread(sel_path)
  # base-R subsetting to avoid data.table NSE/column-name collisions
  best <- sel$threshold[sel$trait == trait_name & sel$best == TRUE]
  if (length(best) == 1 && best %in% avail_cols) best else fallback
}

gw_t2d_col <- select_best_threshold("T2D", t2d_thresh_cols)
gw_cad_col <- select_best_threshold("CAD", cad_thresh_cols)
cat(sprintf("  GW T2D threshold (best-fit on eur_train): %s\n", gw_t2d_col))
cat(sprintf("  GW CAD threshold (best-fit on eur_train): %s\n", gw_cad_col))

setnames(gw_t2d, gw_t2d_col, "GW_T2D_PRS")
setnames(gw_cad, gw_cad_col, "GW_CAD_PRS")

merged <- merge(merged, gw_t2d[, .(FID, IID, GW_T2D_PRS)], by = c("FID", "IID"), all.x = TRUE)
merged <- merge(merged, gw_cad[, .(FID, IID, GW_CAD_PRS)], by = c("FID", "IID"), all.x = TRUE)

# Synthetic outcome
merged[, Synthetic := fifelse(T2D == 1 & CAD == 1, 1L,
                     fifelse(T2D == 0 | CAD == 0, 0L, NA_integer_))]

# EUR validation only
eur_val <- merged[group == "eur_validation"]
cat(sprintf("  EUR validation: %d individuals\n", nrow(eur_val)))
cat(sprintf("  T2D cases: %d, CAD cases: %d, Synthetic cases: %d\n",
            sum(eur_val$T2D == 1, na.rm = TRUE),
            sum(eur_val$CAD == 1, na.rm = TRUE),
            sum(eur_val$Synthetic == 1, na.rm = TRUE)))

if ("Batch" %in% colnames(eur_val)) {
  eur_val[, Batch := as.factor(Batch)]
}

# ============================================================
# Model fitting
# ============================================================
cat("\n--- Fitting models ---\n")

covar_terms <- c("age", "age2", "sex", "Batch", paste0("PC", 1:10))
outcomes <- c("T2D", "CAD", "Synthetic")
outcome_labels <- c(T2D = "Type 2 Diabetes", CAD = "Coronary Artery Disease",
                    Synthetic = "T2D and CAD")

models <- list(
  list(name = "Joint Cluster PRS", predictors = prs_cols),
  list(name = "GW T2D PRS",        predictors = "GW_T2D_PRS"),
  list(name = "GW CAD PRS",        predictors = "GW_CAD_PRS"),
  list(name = "Combined",          predictors = c(prs_cols, "GW_T2D_PRS", "GW_CAD_PRS"))
)

results_list <- list()

for (outcome in outcomes) {
  cat(sprintf("\n  Outcome: %s\n", outcome))

  all_prs <- c(prs_cols, "GW_T2D_PRS", "GW_CAD_PRS")
  model_data <- eur_val[, c(outcome, all_prs, covar_terms), with = FALSE]
  model_data <- na.omit(model_data)

  n <- nrow(model_data)
  n_cases <- sum(model_data[[outcome]] == 1)
  n_controls <- sum(model_data[[outcome]] == 0)
  cat(sprintf("    N=%d (cases=%d, controls=%d)\n", n, n_cases, n_controls))

  # Standardize all PRS columns
  for (pc in all_prs) {
    model_data[[pc]] <- scale(model_data[[pc]])[, 1]
  }

  # Null model
  formula_null <- as.formula(paste(outcome, "~", paste(covar_terms, collapse = " + ")))
  fit_null <- glm(formula_null, data = model_data, family = binomial)

  for (mod in models) {
    rhs <- paste(c(mod$predictors, covar_terms), collapse = " + ")
    formula_full <- as.formula(paste(outcome, "~", rhs))

    fit_full <- tryCatch(
      glm(formula_full, data = model_data, family = binomial),
      error = function(e) { cat(sprintf("    ERROR %s: %s\n", mod$name, e$message)); NULL }
    )
    if (is.null(fit_full)) next

    r2 <- nagelkerke_r2(fit_full, fit_null, n)
    pred_probs <- predict(fit_full, type = "response")
    auc_val <- tryCatch(
      as.numeric(auc(roc(model_data[[outcome]], pred_probs, quiet = TRUE))),
      error = function(e) NA_real_
    )

    cat(sprintf("    %s: R²=%.5f, AUC=%.4f\n", mod$name, r2, auc_val))

    results_list[[length(results_list) + 1]] <- data.table(
      outcome = outcome,
      outcome_label = outcome_labels[[outcome]],
      model = mod$name,
      nagelkerke_r2 = round(r2, 6),
      AUC = round(auc_val, 4),
      n_cases = n_cases,
      n_controls = n_controls,
      n_total = n
    )
  }
}

results <- rbindlist(results_list)

# ============================================================
# Write CSV
# ============================================================
csv_path <- file.path(results_dir, "joint_vs_gw_comparison.csv")
fwrite(results, csv_path)
cat(sprintf("\n--- Saved CSV: %s ---\n", csv_path))

# ============================================================
# Bar chart
# ============================================================
cat("\n--- Generating comparison bar chart ---\n")

model_colors <- c(
  "Joint Cluster PRS" = "#1F78B4",
  "GW T2D PRS"        = "#E74C3C",
  "GW CAD PRS"        = "#3498DB",
  "Combined"           = "#9B59B6"
)

results[, model := factor(model, levels = names(model_colors))]
results[, outcome_label := factor(outcome_label,
  levels = c("Type 2 Diabetes", "Coronary Artery Disease", "T2D and CAD"))]

plot_long <- melt(results,
                  id.vars = c("outcome", "outcome_label", "model"),
                  measure.vars = c("nagelkerke_r2", "AUC"),
                  variable.name = "metric", value.name = "value")

plot_long[, metric_label := fifelse(metric == "nagelkerke_r2",
                                    "Nagelkerke R²", "AUC")]
plot_long[, metric_label := factor(metric_label, levels = c("Nagelkerke R²", "AUC"))]

p <- ggplot(plot_long, aes(x = outcome_label, y = value, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", value)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = model_colors, name = "Model") +
  labs(x = NULL, y = NULL,
       title = "Joint Cluster PRS vs Genome-Wide PRS: Model Comparison",
       subtitle = "EUR Validation") +
  theme_minimal(base_size = 13) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(size = 10, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 13)
  )

# Expand y-axis to make room for text labels
p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

out_png <- file.path(results_dir, "joint_vs_gw_comparison.png")
ggsave(out_png, p, width = 12, height = 6, dpi = 300)
cat(sprintf("  Saved: %s\n", out_png))

# Print summary
cat("\n=== Summary ===\n")
print(results[, .(outcome_label, model, nagelkerke_r2, AUC)], nrow = Inf)

cat("\n=== Done ===\n")
