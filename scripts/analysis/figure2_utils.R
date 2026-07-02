# figure2_utils.R
# Shared utilities for Figure 2 (PRS analysis):
#   data loading, PRS metadata parsing, and statistical helpers.

# Assumes figure_utils.R has been sourced (for theme_big_text).

library(tidyverse)


# ============================================================================
# Data Loading
# ============================================================================

#' Load and merge phenotype + PRS data, build PRS metadata
#'
#' @param pheno_path Path to phenotype CSV
#' @param prs_path   Path to PRS TSV (output of compute_prs.R)
#' @return List with:
#'   - dat: tibble with merged phenotype + standardized PRS columns
#'   - prs_meta: tibble describing each PRS column
load_prs_data <- function(pheno_path, prs_path) {

  cat("Loading phenotype data...\n")
  pheno <- read_csv(pheno_path, show_col_types = FALSE)
  cat(sprintf("  Phenotype: %d individuals x %d columns\n", nrow(pheno), ncol(pheno)))

  cat("Loading PRS data...\n")
  prs <- read_tsv(prs_path, show_col_types = FALSE)
  cat(sprintf("  PRS: %d individuals x %d columns\n", nrow(prs), ncol(prs)))

  # Inner join on person_id
  dat <- inner_join(pheno, prs, by = "person_id")
  cat(sprintf("  After merge: %d individuals\n", nrow(dat)))

  # Drop rows missing required covariates
  dat <- dat %>%
    filter(!is.na(age), !is.na(sex_at_birth), sex_at_birth %in% c("Male", "Female"))
  cat(sprintf("  After covariate QC: %d individuals\n", nrow(dat)))

  # Derive confluence outcome: 1 if T2D or CAD, else 0
  dat <- dat %>%
    mutate(confluence = as.integer(t2d == 1L | cad == 1L))

  cat(sprintf("  T2D prevalence:        %.2f%% (%d cases)\n",
    100 * mean(dat$t2d), sum(dat$t2d)))
  cat(sprintf("  CAD prevalence:        %.2f%% (%d cases)\n",
    100 * mean(dat$cad), sum(dat$cad)))
  cat(sprintf("  Confluence prevalence: %.2f%% (%d cases)\n",
    100 * mean(dat$confluence), sum(dat$confluence)))

  # --- Build PRS metadata ---
  prs_cols <- setdiff(names(prs), "person_id")
  prs_meta <- map_dfr(prs_cols, function(col) {
    # Parse column name: prs_k{N}_{anc} or grs_total_{anc}
    is_grs <- str_detect(col, "^grs_total_")
    anc <- toupper(str_extract(col, "[^_]+$"))
    if (is_grs) {
      label <- paste(anc, "GRS")
    } else {
      # Extract cluster number from prs_k{N}
      k_num <- str_replace(col, "^prs_k(\\d+)_.*", "\\1")
      label <- paste(anc, paste0("K", k_num))
    }
    tibble(prs_col = col, gwas_ancestry = anc, label = label, is_grs = is_grs)
  })

  cat(sprintf("  PRS columns: %d (%d cluster, %d GRS total)\n",
    nrow(prs_meta), sum(!prs_meta$is_grs), sum(prs_meta$is_grs)))

  # --- Drop PRS columns that are entirely NA (cluster doesn't exist) ---
  has_data <- map_lgl(prs_meta$prs_col, ~ !all(is.na(dat[[.x]])))
  n_dropped <- sum(!has_data)
  if (n_dropped > 0) {
    cat(sprintf("  Dropping %d PRS columns with no data: %s\n",
      n_dropped, paste(prs_meta$label[!has_data], collapse = ", ")))
    prs_meta <- prs_meta[has_data, ]
  }

  # --- Standardize PRS columns ---
  dat <- dat %>%
    mutate(across(all_of(prs_meta$prs_col), ~ {
      if (all(is.na(.x))) .x
      else (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)
    }))

  list(dat = dat, prs_meta = prs_meta)
}


# ============================================================================
# Statistical Helpers
# ============================================================================

#' Fit logistic regression for one PRS column and extract OR + CI
#'
#' @param dat       tibble with outcome, PRS, and covariates
#' @param outcome   Name of the binary outcome column (e.g., "t2d")
#' @param prs_col   Name of the (standardized) PRS column
#' @return One-row tibble with OR, CI_lo, CI_hi, p_value, n_cases, n_total
fit_prs_or <- function(dat, outcome, prs_col) {

  # Build formula: outcome ~ prs + age + age^2 + sex + 16 PCs
  pc_terms <- paste0("PC", 1:16)
  rhs <- paste(c(prs_col, "age", "I(age^2)", "sex_at_birth", pc_terms), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  # Subset to complete cases for this PRS column
  keep_cols <- c(outcome, prs_col, "age", "sex_at_birth", pc_terms)
  sub_dat <- dat %>% drop_na(all_of(keep_cols))

  fit <- tryCatch(
    glm(fml, data = sub_dat, family = binomial),
    error = function(e) NULL
  )

  if (is.null(fit) || !fit$converged) {
    return(tibble(
      OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_,
      p_value = NA_real_, n_cases = sum(sub_dat[[outcome]]),
      n_total = nrow(sub_dat), converged = FALSE
    ))
  }

  # Extract coefficient for the PRS term (first non-intercept)
  coef_idx <- 2  # PRS is the first predictor after intercept
  beta <- coef(fit)[coef_idx]
  se <- summary(fit)$coefficients[coef_idx, "Std. Error"]
  p_val <- summary(fit)$coefficients[coef_idx, "Pr(>|z|)"]

  # Wald CI (fast, appropriate for large samples)
  ci_lo <- beta - 1.96 * se
  ci_hi <- beta + 1.96 * se

  tibble(
    OR      = exp(beta),
    CI_lo   = exp(ci_lo),
    CI_hi   = exp(ci_hi),
    p_value = p_val,
    n_cases = sum(sub_dat[[outcome]]),
    n_total = nrow(sub_dat),
    converged = TRUE
  )
}


#' Compute Nagelkerke pseudo-R-squared from a glm object
#'
#' @param model A fitted glm object (family = binomial)
#' @return Scalar Nagelkerke R-squared
nagelkerke_r2 <- function(model) {
  n <- nobs(model)
  # Use deviance values already stored in the glm object to avoid
  # update(model, . ~ 1) scoping issues inside lapply closures.
  # glm stores: deviance = -2*logLik(full), null.deviance = -2*logLik(intercept-only)
  r2_cs  <- 1 - exp((model$deviance - model$null.deviance) / n)
  r2_max <- 1 - exp(-model$null.deviance / n)
  r2_cs / r2_max
}


#' Convert observed-scale R-squared to liability-scale (Lee et al. 2012)
#'
#' @param r2_obs Observed Nagelkerke R-squared (incremental)
#' @param K      Population prevalence
#' @param P      Sample case proportion
#' @return Liability-scale R-squared
liability_r2 <- function(r2_obs, K, P) {
  t_thresh <- qnorm(1 - K)
  z <- dnorm(t_thresh)
  r2_obs * K^2 * (1 - K)^2 / (z^2 * P * (1 - P))
}


#' Run logistic regressions for all PRS x outcome combinations
#'
#' @param dat       tibble from load_prs_data()
#' @param prs_meta  tibble from load_prs_data()
#' @param outcomes  Named character vector of outcome columns
#'                  (e.g., c("Type 2 Diabetes" = "t2d", "CAD" = "cad", ...))
#' @return tibble with columns: outcome_label, prs_col, label,
#'         gwas_ancestry, is_grs, OR, CI_lo, CI_hi, p_value, n_cases, n_total
run_all_regressions <- function(dat, prs_meta,
    outcomes = c("Type 2 Diabetes" = "t2d",
                 "CAD" = "cad",
                 "Confluence" = "confluence")) {

  map_dfr(names(outcomes), function(outcome_label) {
    outcome_col <- outcomes[[outcome_label]]
    cat(sprintf("  Fitting regressions for %s...\n", outcome_label))

    map_dfr(seq_len(nrow(prs_meta)), function(i) {
      res <- fit_prs_or(dat, outcome_col, prs_meta$prs_col[i])
      res %>% mutate(
        outcome_label = outcome_label,
        prs_col       = prs_meta$prs_col[i],
        label         = prs_meta$label[i],
        gwas_ancestry = prs_meta$gwas_ancestry[i],
        is_grs        = prs_meta$is_grs[i]
      )
    })
  })
}


#' Compute liability-scale R-squared for all PRS x outcome combinations
#'
#' @param dat       tibble from load_prs_data()
#' @param prs_meta  tibble from load_prs_data()
#' @param outcomes  Named character vector of outcome columns
#' @return tibble with columns: outcome_label, prs_col, label,
#'         gwas_ancestry, is_grs, r2_obs, r2_liability, prevalence
run_all_liability_r2 <- function(dat, prs_meta,
    outcomes = c("Type 2 Diabetes" = "t2d",
                 "CAD" = "cad",
                 "Confluence" = "confluence")) {

  pc_terms <- paste0("PC", 1:16)
  covar_rhs <- paste(c("age", "I(age^2)", "sex_at_birth", pc_terms), collapse = " + ")

  map_dfr(names(outcomes), function(outcome_label) {
    outcome_col <- outcomes[[outcome_label]]
    cat(sprintf("  Computing liability R-squared for %s...\n", outcome_label))

    # Fit null model once per outcome (covariates only)
    null_fml <- as.formula(paste(outcome_col, "~", covar_rhs))
    keep_cols_null <- c(outcome_col, "age", "sex_at_birth", pc_terms)
    dat_complete_null <- dat %>% drop_na(all_of(keep_cols_null))
    null_fit <- glm(null_fml, data = dat_complete_null, family = binomial)
    r2_null <- nagelkerke_r2(null_fit)

    P <- mean(dat_complete_null[[outcome_col]], na.rm = TRUE)
    K <- P  # use observed prevalence as population prevalence

    map_dfr(seq_len(nrow(prs_meta)), function(i) {
      prs_col_name <- prs_meta$prs_col[i]

      # Full model: covariates + PRS
      full_rhs <- paste(c(prs_col_name, "age", "I(age^2)", "sex_at_birth", pc_terms),
                        collapse = " + ")
      full_fml <- as.formula(paste(outcome_col, "~", full_rhs))

      keep_cols <- c(outcome_col, prs_col_name, "age", "sex_at_birth", pc_terms)
      sub_dat <- dat %>% drop_na(all_of(keep_cols))

      full_fit <- tryCatch(
        glm(full_fml, data = sub_dat, family = binomial),
        error = function(e) NULL
      )

      if (is.null(full_fit) || !full_fit$converged) {
        return(tibble(
          outcome_label = outcome_label, prs_col = prs_col_name,
          label = prs_meta$label[i], gwas_ancestry = prs_meta$gwas_ancestry[i],
          is_grs = prs_meta$is_grs[i],
          r2_obs = NA_real_, r2_liability = NA_real_, prevalence = P
        ))
      }

      # Refit null on same subset for fair comparison
      null_fit_sub <- glm(null_fml, data = sub_dat, family = binomial)
      r2_null_sub <- nagelkerke_r2(null_fit_sub)
      r2_full <- nagelkerke_r2(full_fit)

      r2_incremental <- r2_full - r2_null_sub
      r2_liab <- liability_r2(max(r2_incremental, 0), K, P)

      tibble(
        outcome_label = outcome_label, prs_col = prs_col_name,
        label = prs_meta$label[i], gwas_ancestry = prs_meta$gwas_ancestry[i],
        is_grs = prs_meta$is_grs[i],
        r2_obs = r2_incremental, r2_liability = r2_liab, prevalence = P
      )
    })
  })
}


#' Split data into training and validation sets (stratified 50/50)
#'
#' Stratifies on the joint distribution of t2d and cad to preserve
#' outcome prevalence in both splits.
#'
#' @param dat       tibble from load_prs_data()
#' @param seed      Random seed (default: 47)
#' @param train_frac Fraction allocated to training (default: 0.5)
#' @return List with:
#'   - train: tibble (training set)
#'   - val:   tibble (validation set)
split_train_val <- function(dat, seed = 47, train_frac = 0.5) {
  set.seed(seed)

  dat <- dat %>% mutate(.strat = paste(t2d, cad, sep = "_"))

  train_idx <- dat %>%
    mutate(.row = row_number()) %>%
    group_by(.strat) %>%
    slice_sample(prop = train_frac) %>%
    ungroup() %>%
    pull(.row)

  train <- dat %>% slice(train_idx) %>% select(-.strat)
  val   <- dat %>% slice(-train_idx) %>% select(-.strat)

  cat(sprintf("Train/val split (seed=%d): %d train, %d val\n",
    seed, nrow(train), nrow(val)))
  cat(sprintf("  Train T2D prevalence: %.2f%%\n", 100 * mean(train$t2d)))
  cat(sprintf("  Val   T2D prevalence: %.2f%%\n", 100 * mean(val$t2d)))
  cat(sprintf("  Train CAD prevalence: %.2f%%\n", 100 * mean(train$cad)))
  cat(sprintf("  Val   CAD prevalence: %.2f%%\n", 100 * mean(val$cad)))

  list(train = train, val = val)
}
