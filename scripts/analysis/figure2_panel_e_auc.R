# figure2_panel_e_auc.R
# Panel E: AUC dot plot comparing individual cluster PRS, GRS, and combined
# models, evaluated on a held-out validation set.
#
# Models are fit on a 50/50 training split and AUC is computed on validation.
# Each model includes standard covariates (age, age^2, sex, PC1-PC16).
#
# For each of EUR and META GWAS ancestries:
#   1. outcome ~ K1 + covariates
#   2. outcome ~ K2 + covariates
#   3. outcome ~ GRS + covariates
#   4. outcome ~ K1 + K2 + GRS + covariates
#
# Assumes figure_utils.R and figure2_utils.R have been sourced.

library(pROC)


#' Build model specifications for AUC comparison
#'
#' For a given GWAS ancestry, creates specs for:
#'   - Each individual cluster PRS (K1, K2, ...)
#'   - GRS total
#'   - Combined (all cluster PRS + GRS)
#'
#' @param prs_meta  tibble with PRS column metadata
#' @param gwas_anc  GWAS ancestry string (e.g., "EUR", "META")
#' @return List of lists, each with: model_label, prs_cols, gwas_ancestry
build_model_specs <- function(prs_meta, gwas_anc) {
  anc_meta <- prs_meta %>% filter(gwas_ancestry == gwas_anc)
  cluster_meta <- anc_meta %>% filter_out(is_grs) %>% arrange(prs_col)
  grs_meta     <- anc_meta %>% filter(is_grs)

  specs <- list()

  # Individual cluster models
  for (i in seq_len(nrow(cluster_meta))) {
    specs[[length(specs) + 1]] <- list(
      model_label   = paste0("K", i),
      prs_cols      = cluster_meta$prs_col[i],
      gwas_ancestry = gwas_anc
    )
  }

  # GRS total model
  if (nrow(grs_meta) > 0) {
    specs[[length(specs) + 1]] <- list(
      model_label   = "GRS",
      prs_cols      = grs_meta$prs_col[1],
      gwas_ancestry = gwas_anc
    )
  }

  # Combined model (all cluster PRS + GRS)
  all_cols <- c(cluster_meta$prs_col, grs_meta$prs_col)
  if (length(all_cols) > 1) {
    specs[[length(specs) + 1]] <- list(
      model_label   = "Combined",
      prs_cols      = all_cols,
      gwas_ancestry = gwas_anc
    )
  }

  specs
}


#' Fit a logistic model on training data and compute AUC on validation data
#'
#' @param train     Training tibble
#' @param val       Validation tibble
#' @param outcome   Name of the binary outcome column
#' @param prs_cols  Character vector of PRS column(s) to include
#' @param ci_level  Confidence interval level (default: 0.95)
#' @return One-row tibble with auc, auc_lo, auc_hi, n_train, n_val
compute_model_auc <- function(train, val, outcome, prs_cols, ci_level = 0.95) {

  pc_terms <- paste0("PC", 1:16)
  rhs <- paste(c(prs_cols, "age", "I(age^2)", "sex_at_birth", pc_terms),
               collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  keep_cols <- c(outcome, prs_cols, "age", "sex_at_birth", pc_terms)
  train_sub <- train %>% drop_na(all_of(keep_cols))
  val_sub   <- val   %>% drop_na(all_of(keep_cols))

  fit <- tryCatch(
    glm(fml, data = train_sub, family = binomial),
    error = function(e) NULL
  )

  if (is.null(fit) || !fit$converged) {
    return(tibble(
      auc = NA_real_, auc_lo = NA_real_, auc_hi = NA_real_,
      n_train = nrow(train_sub), n_val = nrow(val_sub)
    ))
  }

  # Predict on validation set
  pred_prob <- predict(fit, newdata = val_sub, type = "response")

  # Compute AUC with DeLong CI
  roc_obj <- tryCatch(
    roc(val_sub[[outcome]], pred_prob, quiet = TRUE),
    error = function(e) NULL
  )

  if (is.null(roc_obj)) {
    return(tibble(
      auc = NA_real_, auc_lo = NA_real_, auc_hi = NA_real_,
      n_train = nrow(train_sub), n_val = nrow(val_sub)
    ))
  }

  ci <- ci.auc(roc_obj, conf.level = ci_level, method = "delong")

  tibble(
    auc     = as.numeric(auc(roc_obj)),
    auc_lo  = ci[1],
    auc_hi  = ci[3],
    n_train = nrow(train_sub),
    n_val   = nrow(val_sub)
  )
}


#' Create Panel E: AUC dot plot with lineranges
#'
#' Splits data 50/50, fits models on training, computes AUC on validation.
#' Models: individual K1, K2, GRS, and Combined (K1+K2+GRS) for EUR and META.
#'
#' @param dat       tibble with standardized PRS columns and outcomes
#' @param prs_meta  tibble with PRS column metadata
#' @param outcomes  Named character vector (display name = column name)
#' @param seed      Random seed for train/val split (default: 47)
#' @param show_ci   Whether to show 95% CI error bars (default: TRUE)
#' @param ci_level  Confidence interval level (default: 0.95)
#' @return ggplot object
plot_panel_e_auc <- function(dat, prs_meta,
    outcomes = c("Type 2 Diabetes" = "t2d",
                 "CAD" = "cad",
                 "Confluence" = "confluence"),
    seed = 47,
    show_ci = TRUE,
    ci_level = 0.95) {

  # Split data
  splits <- split_train_val(dat, seed = seed)
  train  <- splits$train
  val    <- splits$val

  # Build model specs for EUR and META
  gwas_ancestries <- c("EUR", "META")
  all_specs <- list()
  for (anc in gwas_ancestries) {
    anc_in_meta <- prs_meta %>% filter(gwas_ancestry == anc)
    if (nrow(anc_in_meta) > 0) {
      all_specs <- c(all_specs, build_model_specs(prs_meta, anc))
    }
  }

  # Run all models across outcomes
  cat("Computing AUC for Panel E...\n")
  auc_results <- map_dfr(names(outcomes), function(outcome_label) {
    outcome_col <- outcomes[[outcome_label]]
    cat(sprintf("  Fitting models for %s...\n", outcome_label))

    map_dfr(all_specs, function(spec) {
      res <- compute_model_auc(train, val, outcome_col, spec$prs_cols,
                               ci_level = ci_level)
      res %>% mutate(
        outcome_label = outcome_label,
        model_label   = spec$model_label,
        gwas_ancestry = spec$gwas_ancestry
      )
    })
  })

  # Drop failed models
  auc_results <- auc_results %>% filter_out(is.na(auc))

  # Factor ordering
  model_order <- c("K1", "K2", "GRS", "Combined")
  model_order <- intersect(model_order, unique(auc_results$model_label))
  auc_results <- auc_results %>%
    mutate(
      model_label   = factor(model_label, levels = model_order),
      outcome_label = factor(outcome_label, levels = names(outcomes))
    )

  # Ancestry colors (consistent with all other panels)
  anc_colors <- c("EUR" = "#4A90D9", "AFR" = "#D94A4A", "META" = "#50B878")

  p <- ggplot(auc_results,
              aes(x = model_label, y = auc, color = gwas_ancestry)) +
    geom_point(
      position = position_dodge(width = 0.5),
      size = 3
    ) +
    geom_linerange(
      aes(ymin = auc_lo, ymax = auc_hi),
      position = position_dodge(width = 0.5),
      linewidth = 0.8
    ) +
    facet_wrap(~ outcome_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = anc_colors, name = "GWAS\nAncestry") +
    theme_big_text() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      strip.text  = element_text(face = "bold")
    ) +
    labs(
      x = "Model",
      y = "AUC (validation set)"
    )
  
  p
}
