#!/usr/bin/env Rscript
# 04b_liability_h2.R — Convert LDSC observed-scale SNP heritabilities to the
# liability scale (Lee et al. 2011, AJHG) for the a0 EUR traits.
#
# Liability-scale transformation:
#   h2_liab = h2_obs * K^2 (1-K)^2 / ( P (1-P) z^2 )
#   z = dnorm(qnorm(1 - K))           # height of the std-normal density at the
#                                     # liability threshold t = qnorm(1 - K)
# where
#   K = population prevalence       (config: pop_prevalence; UKB-derived here)
#   P = sample (case) prevalence    = n_cases / (n_cases + n_controls)
# The transformation is linear in h2_obs, so the SE scales by the same factor.
#
# IMPORTANT — P must match the N convention used when each trait was munged:
#   * traits munged on TOTAL N  -> P = n_cases / n_total   (the case here)
#   * traits munged on EFFECTIVE N would instead use P = 0.5
# In this pipeline every trait is munged on total N (Suzuki T2D is forced to
# total N via the `n_total` override in 00_preprocess_sumstats.py), so P is
# always n_cases / n_total.
#
# Observed-scale h2 is read straight from the LDSC --rg logs:
#   * each comparison trait  -> "Heritability of phenotype 2/2" block of rg_<trait>.log
#   * the reference trait     -> "Heritability of phenotype 1"   block (identical
#                                across all logs; read from the first one)
#
# Implemented in base R (+ yaml) so it runs under plain Rscript without the
# rocky9 tidyverse/GLIBCXX environment recipe.
#
# Usage:
#   Rscript scripts/a0_analysis/04b_liability_h2.R \
#       --config config/a0_config.yaml \
#       --results-dir results/a0_analysis
#
# Outputs:
#   results/a0_analysis/liability_h2.csv          (one row per trait)
#   augments results/a0_analysis/rg_results.csv   (adds liability-scale columns)

suppressMessages(library(yaml))

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
config_path <- get_arg("--config", "config/a0_config.yaml")
results_dir <- get_arg("--results-dir", "results/a0_analysis")

cfg        <- yaml::read_yaml(config_path)
traits_cfg <- cfg$traits

# ---------------------------------------------------------------------------
# Parse observed-scale h2 (+ SE) from an LDSC --rg log.
# `phenotype` is 1 (reference) or 2 (comparison). Grabs the first
# "Total Observed scale h2: <est> (<se>)" line at/after the matching header —
# robust to the trailing "Ratio: x (y)" vs "Ratio < 0 ..." line variants.
# ---------------------------------------------------------------------------
parse_obs_h2 <- function(log_path, phenotype) {
  txt    <- readLines(log_path, warn = FALSE)
  header <- if (phenotype == 1) "Heritability of phenotype 1" else "Heritability of phenotype 2/2"
  h_idx  <- grep(header, txt, fixed = TRUE)
  if (length(h_idx) == 0L) stop(sprintf("'%s' not found in %s", header, log_path))
  h2_lines <- grep("Total Observed scale h2:", txt, fixed = TRUE)
  h2_idx   <- h2_lines[h2_lines >= h_idx[1L]][1L]
  m <- regmatches(
    txt[h2_idx],
    regexec("Total Observed scale h2:\\s*([-0-9.eE]+)\\s*\\(([-0-9.eE]+)\\)", txt[h2_idx])
  )[[1L]]
  list(h2_obs = as.numeric(m[2L]), h2_obs_se = as.numeric(m[3L]))
}

# ---------------------------------------------------------------------------
# Identify reference vs comparison traits
# ---------------------------------------------------------------------------
keys      <- names(traits_cfg)
is_ref    <- vapply(traits_cfg, function(x) isTRUE(x$is_reference), logical(1))
ref_key   <- keys[is_ref][1L]
comp_keys <- keys[!is_ref]

# Reference h2: read from the first comparison's log (phenotype 1 block)
ref_log <- file.path(results_dir, sprintf("rg_%s.log", comp_keys[1L]))
ref_h2  <- parse_obs_h2(ref_log, phenotype = 1)

# ---------------------------------------------------------------------------
# Assemble per-trait table (reference first, then comparisons in config order)
# ---------------------------------------------------------------------------
ordered_keys <- c(ref_key, comp_keys)
rows <- lapply(ordered_keys, function(key) {
  info <- traits_cfg[[key]]
  obs  <- if (key == ref_key) ref_h2 else
            parse_obs_h2(file.path(results_dir, sprintf("rg_%s.log", key)), phenotype = 2)
  n_cases    <- as.numeric(info$n_cases)
  n_controls <- as.numeric(info$n_controls)
  K          <- as.numeric(info$pop_prevalence)
  P          <- n_cases / (n_cases + n_controls)
  z          <- dnorm(qnorm(1 - K))
  conv       <- K^2 * (1 - K)^2 / (P * (1 - P) * z^2)
  data.frame(
    trait        = key,
    label        = info$label,
    is_reference = key == ref_key,
    n_cases      = n_cases,
    n_controls   = n_controls,
    n_total      = n_cases + n_controls,
    samp_prev_P  = P,
    pop_prev_K   = K,
    z            = z,
    conv_factor  = conv,
    h2_obs       = obs$h2_obs,
    h2_obs_se    = obs$h2_obs_se,
    h2_liab      = obs$h2_obs * conv,
    h2_liab_se   = obs$h2_obs_se * conv,
    stringsAsFactors = FALSE
  )
})
tab <- do.call(rbind, rows)

# Round numeric columns for the written CSV / console echo
round_df <- function(df, digits) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], signif, digits = digits)
  df
}

out_csv <- file.path(results_dir, "liability_h2.csv")
write.csv(round_df(tab, 6), out_csv, row.names = FALSE)
message("Wrote: ", out_csv)
print(round_df(tab, 4), row.names = FALSE)

# ---------------------------------------------------------------------------
# Augment rg_results.csv with liability-scale h2 (reference = trait1,
# comparison = trait2), matched on label.
# ---------------------------------------------------------------------------
rg_path <- file.path(results_dir, "rg_results.csv")
if (file.exists(rg_path)) {
  rg  <- read.csv(rg_path, stringsAsFactors = FALSE, check.names = FALSE)
  lut <- tab[, c("label", "h2_liab", "h2_liab_se", "samp_prev_P", "pop_prev_K")]

  t1 <- lut; names(t1) <- c("label", "h2_liab_trait1", "h2_liab_se_trait1", "P_trait1", "K_trait1")
  t2 <- lut; names(t2) <- c("label", "h2_liab_trait2", "h2_liab_se_trait2", "P_trait2", "K_trait2")

  orig_cols <- names(rg)
  rg <- merge(rg, t1, by.x = "trait1", by.y = "label", all.x = TRUE, sort = FALSE)
  rg <- merge(rg, t2, by.x = "trait2", by.y = "label", all.x = TRUE, sort = FALSE)
  # merge() hoists the join keys to the front; restore original column order
  # (trait1, trait2, ...) then append the new liability-scale columns.
  new_cols <- setdiff(names(rg), orig_cols)
  rg <- rg[, c(orig_cols, new_cols)]

  write.csv(rg, rg_path, row.names = FALSE)
  message("Augmented: ", rg_path, " (added liability-scale h2 columns)")
} else {
  message("NOTE: ", rg_path, " not found; skipped augmentation.")
}
