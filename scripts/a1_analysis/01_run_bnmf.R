#!/usr/bin/env Rscript
# 01_run_bnmf.R
# Config-driven bNMF orchestrator for a single ancestry.
# Snakemake invokes this once per ancestry; YAML config supplies all paths.
#
# Usage:
#   Rscript scripts/a1_analysis/01_run_bnmf.R \
#     --config config/a1_config.yaml --ancestry EUR

library(data.table)
library(yaml)

# --- Parse CLI args ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
ancestry    <- NULL
if ("--config" %in% args)   config_path <- args[which(args == "--config") + 1]
if ("--ancestry" %in% args) ancestry    <- args[which(args == "--ancestry") + 1]

# Optional: force a specific K instead of the modal one (sensitivity analysis).
# CLI takes precedence over cfg$bnmf$select_k (set below once cfg is loaded).
select_k <- NULL
if ("--select-k" %in% args) select_k <- as.integer(args[which(args == "--select-k") + 1])

if (is.null(config_path) || is.null(ancestry)) {
  stop("Usage: Rscript 01_run_bnmf.R --config <config.yaml> --ancestry <ANC>")
}

# --- Load config ---
cfg <- read_yaml(config_path)

valid_ancestries <- cfg$ancestries
if (!ancestry %in% valid_ancestries) {
  stop(sprintf("Invalid ancestry: %s. Must be one of: %s",
               ancestry, paste(valid_ancestries, collapse = ", ")))
}

# --- Resolve project root ---
project_root <- getwd()
if (!file.exists(file.path(project_root, "scripts", "a1_analysis", "bnmf_algorithm.R"))) {
  script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(script_file)) {
    project_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
  }
}

# --- Source analysis scripts ---
source(file.path(project_root, "scripts", "a1_analysis", "bnmf_algorithm.R"))
source(file.path(project_root, "scripts", "a1_analysis", "prep_bnmf.R"))
source(file.path(project_root, "scripts", "a1_analysis", "format_results.R"))

# --- Extract parameters (coerce to numeric for YAML string values like 5e-8) ---
n_reps      <- as.integer(cfg$bnmf$nreps)
K_init      <- as.integer(cfg$bnmf$K)
phi_value   <- as.numeric(cfg$bnmf$phi)
maf_thresh  <- as.numeric(cfg$bnmf$maf_threshold)
p_threshold <- as.numeric(cfg$bnmf$p_threshold[[ancestry]])

ref_gwas_paths <- unlist(cfg$ref_gwas[[ancestry]])
trait_gwas     <- cfg$trait_gwas[[ancestry]]
trait_names    <- names(trait_gwas)
trait_paths    <- unlist(trait_gwas)

ref_panel_prefix <- cfg$ld_clump$ref_panel[[ancestry]]
clump_r2 <- as.numeric(cfg$ld_clump$r2)
clump_kb <- as.integer(cfg$ld_clump$kb)
# Ancestries that get a second LD clump across the pooled union of SNPs
union_clump_ancestries <- unlist(cfg$ld_clump$union_clump_ancestries)
union_clump <- ancestry %in% union_clump_ancestries

# Trait missingness parameters (Smith et al. 2024)
miss_cfg <- cfg$trait_missingness
use_trait_missingness <- !is.null(miss_cfg)
if (use_trait_missingness) {
  missingness_threshold <- as.numeric(ifelse(is.null(miss_cfg$missingness_threshold), 0.20, miss_cfg$missingness_threshold))
  proxy_r2_threshold    <- as.numeric(ifelse(is.null(miss_cfg$proxy_r2_threshold), 0.80, miss_cfg$proxy_r2_threshold))
  proxy_ld_window_kb    <- as.integer(ifelse(is.null(miss_cfg$proxy_ld_window_kb), 1000, miss_cfg$proxy_ld_window_kb))
  # trait_correlation_threshold can be a per-ancestry map or a scalar
  tct_raw <- miss_cfg$trait_correlation_threshold
  if (is.list(tct_raw)) {
    trait_corr_threshold <- as.numeric(tct_raw[[ancestry]])
    if (is.null(trait_corr_threshold) || is.na(trait_corr_threshold)) {
      trait_corr_threshold <- 0.85
    }
  } else {
    trait_corr_threshold <- as.numeric(ifelse(is.null(tct_raw), 0.85, tct_raw))
  }
  impute_method         <- ifelse(is.null(miss_cfg$impute_missing), "median", miss_cfg$impute_missing)
}

# Config fallback for select_k (CLI --select-k already parsed above wins)
if (is.null(select_k) && !is.null(cfg$bnmf$select_k)) select_k <- as.integer(cfg$bnmf$select_k)

results_dir <- cfg$results_dir
# When forcing a K, write to a separate <ancestry>_k<N> dir so the default
# (modal-K) results are left intact for side-by-side comparison.
output_subdir <- if (is.null(select_k)) ancestry else sprintf("%s_k%d", ancestry, select_k)
output_dir  <- file.path(project_root, results_dir, output_subdir)
if (!is.null(select_k)) cat(sprintf("select_k = %d -> writing to %s/\n", select_k, output_subdir))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Resolve relative paths
resolve_path <- function(p) {
  if (!startsWith(p, "/")) file.path(project_root, p) else p
}
ref_gwas_paths <- sapply(ref_gwas_paths, resolve_path, USE.NAMES = FALSE)
# Name-preserving copy (CAD_*/T2D_* keys) for allele alignment (Step 3a)
ref_gwas_named <- setNames(ref_gwas_paths, names(unlist(cfg$ref_gwas[[ancestry]])))
trait_paths    <- sapply(trait_paths, resolve_path, USE.NAMES = TRUE)
if (!is.null(ref_panel_prefix)) ref_panel_prefix <- resolve_path(ref_panel_prefix)

# HapMap3 restriction
hapmap3_path <- NULL
if (!is.null(cfg$hapmap3_snp_file)) {
  hapmap3_path <- resolve_path(cfg$hapmap3_snp_file)
  if (!file.exists(hapmap3_path)) {
    cat(sprintf("WARNING: HapMap3 file not found at %s, proceeding without restriction.\n", hapmap3_path))
    hapmap3_path <- NULL
  }
}

# ===== PIPELINE =====
cat(sprintf("\n%s\n=== A1 bNMF: %s (nreps=%d, K=%d, phi=%.1f, p=%.1e) ===\n%s\n",
            strrep("=", 60), ancestry, n_reps, K_init, phi_value, p_threshold, strrep("=", 60)))

# 1b. Resolve fallback traits (null entries use META/EUR GWAS)
fallback_order <- cfg$trait_gwas$fallback_order
if (!is.null(fallback_order)) {
  null_traits <- trait_names[sapply(trait_gwas, is.null)]
  if (length(null_traits) > 0) {
    cat(sprintf("\n--- Step 1b: Resolving %d fallback traits ---\n", length(null_traits)))
    for (tr in null_traits) {
      resolved <- FALSE
      for (fb_anc in fallback_order) {
        fb_path <- cfg$trait_gwas[[fb_anc]][[tr]]
        if (!is.null(fb_path)) {
          fb_path_full <- resolve_path(fb_path)
          if (file.exists(fb_path_full)) {
            trait_paths[tr] <- fb_path_full
            cat(sprintf("  %s: using %s GWAS (%s)\n", tr, fb_anc, basename(fb_path)))
            resolved <- TRUE
            break
          }
        }
      }
      if (!resolved) {
        cat(sprintf("  %s: no fallback found, will be skipped\n", tr))
        trait_paths[tr] <- NA_character_
      }
    }
    # Remove traits that couldn't be resolved
    valid <- !is.na(trait_paths)
    trait_paths <- trait_paths[valid]
    trait_names <- names(trait_paths)
  }
}

# 1. Validate reference GWAS
cat("\n--- Step 1: Validating reference GWAS ---\n")
for (rf in ref_gwas_paths) {
  if (!file.exists(rf)) {
    stop(sprintf("Reference GWAS not found: %s", rf))
  }
  validation <- validate_gwas_format(rf)
  if (!validation$valid) {
    stop(sprintf("Validation failed for %s: %s", basename(rf),
                 paste(validation$messages, collapse = "; ")))
  }
  cat(sprintf("  OK: %s\n", basename(rf)))
}

# 2. QC and filter variants (union from reference GWAS)
cat("\n--- Step 2: Variant QC and filtering ---\n")
if (!is.null(ref_panel_prefix) && !file.exists(paste0(ref_panel_prefix, ".1.bed"))) {
  cat(sprintf("WARNING: LD reference panel not found at %s\n", ref_panel_prefix))
  cat("LD clumping will be SKIPPED.\n")
  ref_panel_prefix <- NULL
}

qc_result <- qc_variants_multi(ref_gwas_paths, p_threshold = p_threshold,
                                maf_threshold = maf_thresh,
                                ref_panel_prefix = ref_panel_prefix,
                                clump_r2 = clump_r2, clump_kb = clump_kb,
                                union_clump = union_clump,
                                # For union-clump ancestries, HapMap3 is applied
                                # at the union level (before the union clump)
                                # inside qc_variants_multi; see block below.
                                union_hapmap3_file = if (union_clump) hapmap3_path else NULL,
                                exclude_regions = cfg$exclude_regions)
filtered <- qc_result$data

# HapMap3 restriction (union-clump ancestries already restricted in the union step)
if (!is.null(hapmap3_path) && !union_clump) {
  hapmap3 <- fread(hapmap3_path, select = "VAR_ID")
  n_before <- nrow(filtered)
  filtered <- filtered[filtered$VAR_ID %in% hapmap3$VAR_ID, ]
  cat(sprintf("HapMap3 restriction: %d -> %d variants\n", n_before, nrow(filtered)))
}

if (nrow(filtered) < 10) {
  stop(sprintf("Only %d variants passed QC (need >= 10). Check p-threshold or data.", nrow(filtered)))
}

# Save QC artifacts
fwrite(qc_result$qc_report, file.path(output_dir, sprintf("qc_report_%s.tsv", ancestry)), sep = "\t")
fwrite(filtered, file.path(output_dir, sprintf("filtered_variants_%s.tsv", ancestry)), sep = "\t")
cat(sprintf("Filtered variants: %d\n", nrow(filtered)))

# 3. Build Z-score matrix from ALL trait GWAS
cat("\n--- Step 3: Building Z-score matrix ---\n")

# Check which trait files exist
missing_traits <- trait_names[!file.exists(trait_paths)]
if (length(missing_traits) > 0) {
  cat(sprintf("WARNING: Missing %d trait files, skipping: %s\n",
              length(missing_traits), paste(missing_traits, collapse = ", ")))
  valid_idx <- file.exists(trait_paths)
  trait_paths <- trait_paths[valid_idx]
  trait_names <- trait_names[valid_idx]
}

trait_list <- as.list(trait_paths)
names(trait_list) <- trait_names

# --- Step 2b: Calculate trait missingness (Smith et al. 2024) ---
proxy_map <- NULL
if (use_trait_missingness) {
  cat(sprintf("\n--- Step 2b: Calculating trait missingness (threshold=%.0f%%) ---\n",
              missingness_threshold * 100))
  missingness_dt <- calculate_trait_missingness(filtered$VAR_ID, trait_list)

  # Identify variants needing proxy search (missingness > threshold)
  target_var_ids <- missingness_dt[missingness > missingness_threshold, VAR_ID]
  cat(sprintf("  Variants with missingness > %.0f%%: %d / %d\n",
              missingness_threshold * 100, length(target_var_ids), nrow(missingness_dt)))

  # Save missingness report
  fwrite(missingness_dt, file.path(output_dir, sprintf("trait_missingness_%s.tsv", ancestry)),
         sep = "\t")

  # --- Step 2c: Find and select proxy variants ---
  if (length(target_var_ids) > 0 && !is.null(ref_panel_prefix)) {
    cat(sprintf("\n--- Step 2c: Finding proxy variants (r2>=%.2f, %dkb window) ---\n",
                proxy_r2_threshold, proxy_ld_window_kb))
    proxy_candidates <- find_proxy_variants(target_var_ids, ref_panel_prefix,
                                             r2_threshold = proxy_r2_threshold,
                                             ld_window_kb = proxy_ld_window_kb)

    if (nrow(proxy_candidates) > 0) {
      proxy_map <- select_best_proxy(proxy_candidates, filtered$VAR_ID, trait_list,
                                      missingness_threshold = missingness_threshold)

      if (nrow(proxy_map) > 0) {
        cat(sprintf("  Selected proxies for %d / %d high-missingness variants\n",
                    nrow(proxy_map), length(target_var_ids)))
        fwrite(proxy_map, file.path(output_dir, sprintf("proxy_map_%s.tsv", ancestry)),
               sep = "\t")
      }
    }
  } else if (length(target_var_ids) > 0) {
    cat("  WARNING: No reference panel available, skipping proxy search.\n")
  }
}

# Build Z-matrix (with proxy support if available)
matrices <- build_z_matrix(filtered, trait_list,
                           proxy_map = proxy_map,
                           impute_method = if (use_trait_missingness) impute_method else "zero")

cat(sprintf("Z-matrix dimensions: %d variants x %d traits\n",
            nrow(matrices$z_matrix), ncol(matrices$z_matrix)))

# --- Step 3a: Allele alignment to disease-risk direction (optional) ---
# Orient each variant row to its CAD/T2D risk-raising allele before the
# non-negative expansion, to prevent mirror clusters (e.g. glucose +/-).
# Must run before Step 3b (correlation filter) and expand_to_nonneg, both
# of which are sign-sensitive.
aa_cfg <- cfg$allele_alignment
if (!is.null(aa_cfg) && isTRUE(aa_cfg$enabled)) {
  rule <- if (is.null(aa_cfg$conflict_rule)) "strongest" else aa_cfg$conflict_rule
  cat(sprintf("\n--- Step 3a: Aligning variants to disease-risk allele (rule=%s) ---\n", rule))
  aln <- align_z_to_disease_risk(matrices$z_matrix, filtered, ref_gwas_named,
                                 conflict_rule = rule)
  matrices$z_matrix <- aln$z_matrix
  fwrite(aln$report,
         file.path(output_dir, sprintf("allele_alignment_report_%s.tsv", ancestry)),
         sep = "\t")
}

# --- Step 3b: Post-matrix trait filtering (significance + correlation) ---
if (use_trait_missingness) {
  cat(sprintf("\n--- Step 3b: Filtering traits (Bonferroni + correlation r>%.2f) ---\n",
              trait_corr_threshold))
  filter_result <- filter_correlated_traits(matrices$z_matrix, trait_list,
                                             filtered$VAR_ID,
                                             correlation_threshold = trait_corr_threshold)
  matrices$z_matrix <- filter_result$z_matrix

  if (length(filter_result$removed_traits) > 0) {
    cat(sprintf("  Removed %d traits: %s\n",
                length(filter_result$removed_traits),
                paste(filter_result$removed_traits, collapse = ", ")))
    # Update n_matrix and trait_names to match
    matrices$n_matrix <- matrices$n_matrix[, colnames(matrices$z_matrix), drop = FALSE]
    trait_names <- colnames(matrices$z_matrix)
  }

  # Save correlation matrix for diagnostics
  cor_dt <- as.data.table(filter_result$correlation_matrix, keep.rownames = "trait")
  fwrite(cor_dt, file.path(output_dir, sprintf("trait_correlations_%s.tsv", ancestry)),
         sep = "\t")

  cat(sprintf("Z-matrix after filtering: %d variants x %d traits\n",
              nrow(matrices$z_matrix), ncol(matrices$z_matrix)))
}

# 4. Expand to non-negative
cat("\n--- Step 4: Expanding to non-negative matrix ---\n")
nonneg_matrix <- expand_to_nonneg(matrices$z_matrix)
cat(sprintf("Non-negative matrix: %d variants x %d columns\n",
            nrow(nonneg_matrix), ncol(nonneg_matrix)))

fwrite(data.table(VAR_ID = rownames(nonneg_matrix), nonneg_matrix),
       file.path(output_dir, sprintf("prepared_matrix_%s.tsv", ancestry)), sep = "\t")

# 5. Run bNMF
cat(sprintf("\n--- Step 5: Running bNMF (nreps=%d, K=%d, phi=%.1f) ---\n",
            n_reps, K_init, phi_value))
results <- run_bnmf(nonneg_matrix, n_reps = n_reps, K = K_init, K0 = K_init,
                    seed = 42, phi = phi_value)

# 6. Summarize and format
cat("\n--- Step 6: Formatting results ---\n")
summary_result <- summarize_bnmf(
  results_list = results,
  trait_names  = trait_names,
  variant_ids  = rownames(nonneg_matrix),
  output_dir   = output_dir,
  ancestry     = ancestry,
  select_k     = select_k
)

# 7. Heatmaps
cat("\n--- Step 7: Generating heatmaps ---\n")
plot_heatmaps(summary_result$W, summary_result$H, output_dir, ancestry = ancestry)

# 8. Summary report
cat("\n--- Step 8: Writing summary report ---\n")
generate_summary_report(summary_result, output_dir, ancestry = ancestry)

cat(sprintf("\n=== Completed %s. Results in: %s ===\n", ancestry, output_dir))
