#!/usr/bin/env Rscript
# ============================================================================
# compute_prs.R
# Compute per-individual PRS from bNMF W_matrix_{ANC} weights using plink2.
#
# For each ancestry, this script:
#   1. Loads W_matrix_{ancestry}.tsv + filtered_variants_{ancestry}.tsv
#   2. Lifts variant coordinates from GRCh37 → GRCh38 (AoU build)
#      (cached to disk — only runs liftOver once)
#   3. Creates a single multi-column plink2 score file (all clusters at once)
#   4. Runs plink2 --score with --score-col-nums across chromosomes 1-22
#   5. Merges results into a single output file
#
# Usage (run from AoU terminal with defaults):
#   Rscript scripts/compute_prs.R
#
# Or override paths:
#   Rscript scripts/compute_prs.R \
#     --geno-dir /custom/path/to/plink \
#     --results-dir /custom/path/to/genotypes \
#     --ancestry EUR
#
# Expected workspace layout (defaults):
#   ~/workspaces/duplicateofstatinresponse/
#   ├── genotype/plink/chr{1..22}.bed/.bim/.fam   (AoU PLINK1 genotypes)
#   └── multiancestry_polygenic/
#       ├── genotypes/{EUR,AFR,META}/              (W_matrix_{ANC} + filtered_variants_{ANC})
#       └── results/prs/{score_files,output,tmp}/  (outputs)
#   ├── genome_useful/{liftOver,hg19ToHg38.over.chain.gz}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  ws <- file.path(Sys.getenv("HOME"), "workspaces", "duplicateofstatinresponse")
  mp <- file.path(ws, "multiancestry_polygenic")

  params <- list(
    geno_dir      = file.path(ws, "genotype", "plink"),
    geno_prefix   = "chr",
    results_dir   = file.path(mp, "genotypes"),
    out_dir       = file.path(mp, "results", "prs"),
    liftover_path = file.path(ws, "genome_useful", "liftOver"),
    chain_path    = file.path(ws, "genome_useful", "hg19ToHg38.over.chain.gz"),
    ancestry      = "all",
    plink2_path   = "plink2",
    apply_weight_cutoff = TRUE
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--geno-dir"      = { params$geno_dir      <- args[i + 1]; i <- i + 2 },
      "--geno-prefix"   = { params$geno_prefix    <- args[i + 1]; i <- i + 2 },
      "--results-dir"   = { params$results_dir    <- args[i + 1]; i <- i + 2 },
      "--out-dir"       = { params$out_dir        <- args[i + 1]; i <- i + 2 },
      "--liftover-path" = { params$liftover_path  <- args[i + 1]; i <- i + 2 },
      "--chain-path"    = { params$chain_path     <- args[i + 1]; i <- i + 2 },
      "--ancestry"      = { params$ancestry       <- args[i + 1]; i <- i + 2 },
      "--plink2-path"   = { params$plink2_path    <- args[i + 1]; i <- i + 2 },
      "--no-weight-cutoff" = { params$apply_weight_cutoff <- FALSE; i <- i + 1 },
      "--help"          = {
        cat("Usage: Rscript compute_prs.R [options]\n\n")
        cat("All arguments are optional (workspace defaults are pre-configured).\n\n")
        cat("  --geno-dir       Directory with per-chr bed/bim/fam files\n")
        cat("                   (default: ~/workspaces/duplicateofstatinresponse/genotype/plink)\n")
        cat("  --geno-prefix    File prefix before chr number (default: 'chr')\n")
        cat("  --results-dir    Directory with {EUR,AFR,META}/ subfolders containing W_matrix_{ANC} + filtered_variants_{ANC}\n")
        cat("                   (default: ~/workspaces/duplicateofstatinresponse/multiancestry_polygenic/genotypes)\n")
        cat("  --out-dir        Output directory for PRS files\n")
        cat("                   (default: ~/workspaces/duplicateofstatinresponse/multiancestry_polygenic/results/prs)\n")
        cat("  --liftover-path  Path to liftOver binary\n")
        cat("  --chain-path     Path to hg19ToHg38 chain file\n")
        cat("  --ancestry       EUR, AFR, META, or 'all' (default: all)\n")
        cat("  --plink2-path    Path to plink2 binary (default: plink2)\n")
        cat("  --no-weight-cutoff  Disable bNMF weight cutoff (use all variants)\n")
        quit(status = 0)
      },
      { cat(sprintf("Unknown argument: %s\n", args[i])); quit(status = 1) }
    )
  }

  params
}

params <- parse_args(args)

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
cat("============================================\n")
cat("PRS Computation from bNMF Weights\n")
cat("============================================\n\n")

cat(sprintf("Genotype directory:  %s\n", params$geno_dir))
cat(sprintf("Genotype prefix:     %s\n", params$geno_prefix))
cat(sprintf("Results directory:   %s\n", params$results_dir))
cat(sprintf("Output directory:    %s\n", params$out_dir))
cat(sprintf("liftOver path:       %s\n", params$liftover_path))
cat(sprintf("Chain file:          %s\n", params$chain_path))
cat(sprintf("Ancestry:            %s\n", params$ancestry))
cat(sprintf("plink2 path:         %s\n", params$plink2_path))
cat(sprintf("Weight cutoff:       %s\n\n",
    if (params$apply_weight_cutoff) "enabled (two-line method)" else "disabled"))

stopifnot(
  "Genotype directory does not exist" = dir.exists(params$geno_dir),
  "Results directory does not exist"  = dir.exists(params$results_dir),
  "liftOver binary not found"         = file.exists(params$liftover_path),
  "Chain file not found"              = file.exists(params$chain_path)
)

dir.create(file.path(params$out_dir, "score_files"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(params$out_dir, "output"),       recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(params$out_dir, "tmp"),           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(params$out_dir, "diagnostics"),   recursive = TRUE, showWarnings = FALSE)

# Determine ancestries to process
if (tolower(params$ancestry) == "all") {
  ancestries <- c("EUR", "AFR", "META")
} else {
  ancestries <- toupper(params$ancestry)
}

# ---------------------------------------------------------------------------
# Helper: liftOver GRCh37 → GRCh38 (with disk cache)
# ---------------------------------------------------------------------------
liftover_variants <- function(dt, liftover_path, chain_path, cache_path, tmp_dir) {

  # Check for cached result (invalidate if missing VAR_ID_38_alt column)
  if (file.exists(cache_path)) {
    cached <- fread(cache_path, nrows = 1)
    if ("VAR_ID_38_alt" %in% names(cached)) {
      cat("  Loading cached liftover from disk...\n")
      cached <- fread(cache_path)
      cat(sprintf("  Loaded %d lifted variants from cache\n", nrow(cached)))
      return(cached)
    } else {
      cat("  Stale cache detected (missing VAR_ID_38_alt), regenerating...\n")
      file.remove(cache_path)
    }
  }

  cat("  Running liftOver GRCh37 → GRCh38 (first run — result will be cached)...\n")

  # Parse VAR_ID into CHR, POS, A1, A2
  dt[, c("CHR", "POS", "A1", "A2") := {
    parts <- tstrsplit(VAR_ID, "_", fixed = TRUE)
    list(parts[[1]], as.integer(parts[[2]]), parts[[3]], parts[[4]])
  }]

  # Write BED file (0-based start, 1-based end)
  bed_in  <- file.path(tmp_dir, "input.bed")
  bed_out <- file.path(tmp_dir, "output.bed")
  bed_unmapped <- file.path(tmp_dir, "unmapped.bed")

  bed_dt <- dt[, .(
    chrom = paste0("chr", CHR),
    start = POS - 1L,
    end   = POS,
    name  = VAR_ID
  )]
  fwrite(bed_dt, bed_in, sep = "\t", col.names = FALSE)

  # Run liftOver
  cmd <- sprintf("%s %s %s %s %s",
    shQuote(liftover_path), shQuote(bed_in), shQuote(chain_path),
    shQuote(bed_out), shQuote(bed_unmapped)
  )
  exit_code <- system(cmd, ignore.stdout = TRUE)
  if (exit_code != 0) {
    stop("liftOver failed with exit code ", exit_code)
  }

  # Read lifted coordinates
  lifted <- fread(bed_out, header = FALSE,
    col.names = c("chrom_38", "start_38", "end_38", "VAR_ID"))
  lifted[, POS_38 := end_38]
  lifted[, CHR_38 := sub("^chr", "", chrom_38)]

  n_unmapped <- nrow(dt) - nrow(lifted)
  cat(sprintf("  Lifted: %d variants, unmapped: %d\n", nrow(lifted), n_unmapped))

  # Merge back
  dt <- merge(dt, lifted[, .(VAR_ID, CHR_38, POS_38)], by = "VAR_ID", all.x = FALSE)

  # Construct GRCh38 variant IDs matching AoU .bim format: chr{CHR}:{POS}:{REF}:{ALT}
  # AoU uses REF:ALT ordering (not alphabetical), so we create both orderings.
  # plink2 will match whichever exists in the .bim and skip the other.
  dt[, VAR_ID_38     := paste0("chr", CHR_38, ":", POS_38, ":", A1, ":", A2)]
  dt[, VAR_ID_38_alt := paste0("chr", CHR_38, ":", POS_38, ":", A2, ":", A1)]

  # Cache to disk
  fwrite(dt, cache_path, sep = "\t")
  cat(sprintf("  Cached lifted variants to %s\n", cache_path))

  # Clean up temp files
  file.remove(bed_in, bed_out)
  if (file.exists(bed_unmapped)) file.remove(bed_unmapped)

  dt
}

# ---------------------------------------------------------------------------
# Helper: Write multi-column score file (all clusters in one file)
# ---------------------------------------------------------------------------
write_multi_score_file <- function(dt, cluster_cols, out_path) {
  # Keep rows where at least one score column has non-zero weight.
  # NOTE: use abs() so variants with negative BETA are not excluded.
  # Without abs(), rowSums > 0 drops all negative-BETA variants since

  # W_k >= 0 (bNMF constraint) makes every W_k*BETA term <= 0 when BETA < 0.
  keep <- rowSums(abs(dt[, cluster_cols, with = FALSE])) > 0
  dt_keep <- dt[keep]
  n_variants <- nrow(dt_keep)

  # Build rows with primary ID ordering (A1:A2)
  primary <- dt_keep[, c("VAR_ID_38", "Effect_Allele", cluster_cols), with = FALSE]
  setnames(primary, c("VAR_ID_38", "Effect_Allele"), c("ID", "A1"))

  # Build rows with alternate ID ordering (A2:A1) — same effect allele & weights
  alt <- dt_keep[, c("VAR_ID_38_alt", "Effect_Allele", cluster_cols), with = FALSE]
  setnames(alt, c("VAR_ID_38_alt", "Effect_Allele"), c("ID", "A1"))

  # Stack both orderings. plink2 will match whichever exists in .bim,
  # skip the other — no double-counting since only one ordering exists per variant.
  score_dt <- rbind(primary, alt)

  fwrite(score_dt, out_path, sep = "\t")
  cat(sprintf("  Score file: %s (%d unique variants x %d clusters, %d rows with both allele orderings)\n",
    basename(out_path), n_variants, length(cluster_cols), nrow(score_dt)))
  invisible(n_variants)
}

# ---------------------------------------------------------------------------
# Helper: Compute bNMF weight cutoff using two-regression-line method
#
#   Pools all non-zero raw W_k values across clusters, sorts them
#   descending, fits a line to the top 1% ("signal") and the bottom 80%
#   ("noise"), then finds the crossover where distance to the signal line
#   becomes shorter than distance to the noise line.
#
#   Reference: Udler et al. 2018 / Kim et al. 2023 — "To define a set of
#   strongest-weighted variants in each cluster and maximize the signal to
#   noise ratio of weights, we developed a method to determine a cluster
#   weight cutoff."
# ---------------------------------------------------------------------------
compute_weight_cutoff <- function(w_dt, cluster_cols,
                                  top_pct = 0.01, bot_pct = 0.80) {

  # Pool all raw W_k values across all clusters, exclude zeros
  all_weights <- unlist(w_dt[, cluster_cols, with = FALSE], use.names = FALSE)
  all_weights <- all_weights[all_weights > 0]
  n_total <- length(all_weights)

  empty_result <- list(
    cutoff = 0, n_total = n_total, n_above = n_total, n_below = 0,
    sorted_weights = sort(all_weights, decreasing = TRUE),
    signal_fit = NULL, noise_fit = NULL, crossover_rank = NA_integer_
  )

  # Need enough points to fit two meaningful lines
  if (n_total < 10) {
    warning("Too few non-zero weights (", n_total, "); skipping cutoff")
    return(empty_result)
  }

  # Sort descending: rank 1 = largest weight
  sorted_w <- sort(all_weights, decreasing = TRUE)
  ranks    <- seq_along(sorted_w)

  # Fit signal line to top 1% (at least 2 points for a valid lm)
  n_top   <- max(2L, ceiling(top_pct * n_total))
  top_idx <- seq_len(n_top)
  signal_fit <- lm(sorted_w[top_idx] ~ ranks[top_idx])

  # Fit noise line to bottom 80%
  bot_start <- max(1L, floor((1 - bot_pct) * n_total) + 1L)
  bot_idx   <- seq(bot_start, n_total)
  noise_fit <- lm(sorted_w[bot_idx] ~ ranks[bot_idx])

  # Perpendicular distance from point (x, y) to line y = a + b*x
  #   Line in general form: b*x - y + a = 0
  #   Distance = |b*x - y + a| / sqrt(b^2 + 1)
  perp_dist <- function(x, y, fit) {
    a <- coef(fit)[[1]]
    b <- coef(fit)[[2]]
    abs(b * x - y + a) / sqrt(b^2 + 1)
  }

  d_signal <- perp_dist(ranks, sorted_w, signal_fit)
  d_noise  <- perp_dist(ranks, sorted_w, noise_fit)

  # Scan from lowest weight upward: find the last rank (highest index)
  # where the point is still closer to the signal line than the noise line
  closer_to_signal <- which(d_signal < d_noise)

  if (length(closer_to_signal) == 0) {
    warning("No crossover found; all points closer to noise line. Setting cutoff = 0")
    empty_result$signal_fit <- signal_fit
    empty_result$noise_fit  <- noise_fit
    return(empty_result)
  }

  crossover_rank <- max(closer_to_signal)
  cutoff <- sorted_w[crossover_rank]

  list(
    cutoff         = cutoff,
    n_total        = n_total,
    n_above        = sum(sorted_w >= cutoff),
    n_below        = sum(sorted_w <  cutoff),
    sorted_weights = sorted_w,
    signal_fit     = signal_fit,
    noise_fit      = noise_fit,
    crossover_rank = crossover_rank
  )
}


# ---------------------------------------------------------------------------
# Helper: Save diagnostic plot of the weight cutoff
#   Uses base R graphics so we don't add a ggplot2 dependency.
# ---------------------------------------------------------------------------
plot_weight_cutoff <- function(cutoff_result, ancestry, out_dir) {

  sorted_w   <- cutoff_result$sorted_weights
  ranks      <- seq_along(sorted_w)
  cutoff     <- cutoff_result$cutoff
  cross_rank <- cutoff_result$crossover_rank
  signal_fit <- cutoff_result$signal_fit
  noise_fit  <- cutoff_result$noise_fit

  out_path <- file.path(out_dir, sprintf("weight_cutoff_%s.png", tolower(ancestry)))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  png(out_path, width = 800, height = 600, res = 150)

  plot(ranks, sorted_w, pch = 16, cex = 0.5,
       xlab = "Rank (descending weight)", ylab = "Raw W_k value",
       main = sprintf("bNMF Weight Cutoff - %s (cutoff = %.4f)", ancestry, cutoff))

  # Overlay fitted lines
  if (!is.null(signal_fit)) {
    abline(signal_fit, col = "blue", lwd = 2, lty = 2)
  }
  if (!is.null(noise_fit)) {
    abline(noise_fit, col = "red", lwd = 2, lty = 2)
  }

  # Cutoff horizontal line + marker
  abline(h = cutoff, col = "darkgreen", lwd = 2, lty = 3)
  if (!is.na(cross_rank)) {
    points(cross_rank, sorted_w[cross_rank],
           pch = 4, cex = 2, col = "darkgreen", lwd = 2)
  }

  legend("topright",
    legend = c("Signal line (top 1%)", "Noise line (bottom 80%)",
               sprintf("Cutoff = %.4f", cutoff)),
    col = c("blue", "red", "darkgreen"),
    lty = c(2, 2, 3), lwd = 2, cex = 0.6
  )

  dev.off()

  cat(sprintf("  Cutoff diagnostic plot saved: %s\n", out_path))
  invisible(out_path)
}


# ---------------------------------------------------------------------------
# Helper: Run plink2 --score with --score-col-nums across chromosomes
#   Scores ALL clusters in a single plink2 call per chromosome.
# ---------------------------------------------------------------------------
run_plink2_score <- function(score_file, n_clusters, geno_dir, geno_prefix,
                             out_prefix, plink2_path) {
  sscore_files <- c()

  # --score-col-nums: columns 3 through 3+n_clusters-1 (1-indexed)
  col_nums <- paste(seq(3, 3 + n_clusters - 1), collapse = ",")

  for (chr in 1:22) {
    bfile <- file.path(geno_dir, paste0(geno_prefix, chr))

    # Check if bed file exists (PLINK1 format)
    if (!file.exists(paste0(bfile, ".bed"))) {
      cat(sprintf("    WARNING: %s.bed not found, skipping chr%d\n", bfile, chr))
      next
    }

    chr_out <- paste0(out_prefix, "_chr", chr)

    cmd <- sprintf(
      "%s --bfile %s --score %s 1 2 header cols=+scoresums ignore-dup-ids --score-col-nums %s --out %s 2>&1",
      shQuote(plink2_path),
      shQuote(bfile),
      shQuote(score_file),
      col_nums,
      shQuote(chr_out)
    )

    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    sscore_path <- paste0(chr_out, ".sscore")
    if (file.exists(sscore_path)) {
      sscore_files <- c(sscore_files, sscore_path)
    }
  }

  cat(sprintf("    Completed %d chromosomes\n", length(sscore_files)))
  sscore_files
}

# ---------------------------------------------------------------------------
# Helper: Merge per-chromosome sscore files (multi-column aware)
# ---------------------------------------------------------------------------
merge_sscore_files <- function(sscore_files) {
  if (length(sscore_files) == 0) {
    warning("No sscore files to merge")
    return(NULL)
  }

  # Read first file to get structure
  combined <- fread(sscore_files[1])

  if (length(sscore_files) > 1) {
    # Identify all score sum columns (SCORE1_SUM, SCORE2_SUM, ...)
    score_cols <- grep("SCORE.*SUM", names(combined), value = TRUE)
    id_col <- ifelse("#IID" %in% names(combined), "#IID",
              ifelse("IID" %in% names(combined), "IID", names(combined)[1]))

    for (f in sscore_files[-1]) {
      dt <- fread(f)
      # Sum each score column across chromosomes
      for (sc in score_cols) {
        if (sc %in% names(dt)) {
          combined[, (sc) := get(sc) + dt[[sc]][match(combined[[id_col]], dt[[id_col]])]]
        }
      }
    }
  }

  combined
}

# ---------------------------------------------------------------------------
# Main: Process each ancestry
# ---------------------------------------------------------------------------
all_prs <- list()

for (anc in ancestries) {
  cat(sprintf("\n--- Processing ancestry: %s ---\n", anc))

  # File paths
  w_path <- file.path(params$results_dir, anc, sprintf("W_matrix_%s.tsv", anc))
  fv_path <- file.path(params$results_dir, anc,
    paste0("filtered_variants_", anc, ".tsv"))

  if (!file.exists(w_path)) {
    cat(sprintf("  WARNING: W_matrix not found at %s, skipping.\n", w_path))
    next
  }
  if (!file.exists(fv_path)) {
    cat(sprintf("  WARNING: filtered_variants not found at %s, skipping.\n", fv_path))
    next
  }

  # Load data
  w_mat <- fread(w_path)
  fv <- fread(fv_path)
  cat(sprintf("  W_matrix: %d variants x %d clusters\n",
    nrow(w_mat), ncol(w_mat) - 1))
  cat(sprintf("  Filtered variants: %d variants\n", nrow(fv)))

  # Get cluster column names (K1, K2, ...)
  cluster_cols <- setdiff(names(w_mat), "VAR_ID")
  cat(sprintf("  Clusters: %s\n", paste(cluster_cols, collapse = ", ")))

  # Join W_matrix with filtered_variants to get Effect_Allele and BETA
  dt <- merge(w_mat, fv[, .(VAR_ID, Effect_Allele, BETA)], by = "VAR_ID")
  cat(sprintf("  After joining: %d variants with effect alleles\n", nrow(dt)))

  # Liftover to GRCh38 (cached) — do this before weighting so the cache
  # only stores coordinate mapping, not computed scores.
  # NOTE: must run BEFORE weight cutoff because the cache stores the full
  # data.table; applying cutoff first would be overwritten by a cache hit.
  cache_path <- file.path(params$out_dir, "score_files",
    sprintf("lifted_%s.tsv", tolower(anc)))
  dt <- liftover_variants(dt, params$liftover_path, params$chain_path,
    cache_path, file.path(params$out_dir, "tmp"))

  # --- bNMF weight cutoff: zero out low-weight (noise) variants per cluster ---
  # Applied to raw W_k values BEFORE BETA multiplication so that noisy variants
  # do not contribute to cluster-specific PRS.  grs_total is unaffected (uses
  # all variants with BETA regardless of cluster weight).
  if (params$apply_weight_cutoff) {
    cutoff_result <- compute_weight_cutoff(dt, cluster_cols)
    cutoff_val <- cutoff_result$cutoff

    if (cutoff_val > 0) {
      cat(sprintf("  Weight cutoff: %.6f (retaining %d / %d non-zero weights)\n",
        cutoff_val, cutoff_result$n_above, cutoff_result$n_total))

      # Zero out W_k values below cutoff (per variant-cluster pair)
      for (k in cluster_cols) {
        n_before <- sum(dt[[k]] > 0)
        dt[get(k) < cutoff_val & get(k) > 0, (k) := 0]
        n_after <- sum(dt[[k]] > 0)
        cat(sprintf("    %s: %d -> %d variants (removed %d)\n",
          k, n_before, n_after, n_before - n_after))
      }

      # Save diagnostic plot
      plot_weight_cutoff(cutoff_result, anc,
        file.path(params$out_dir, "diagnostics"))
    } else {
      cat("  Weight cutoff: not applied (cutoff = 0 or insufficient data)\n")
    }
  }

  # Weight each cluster column by BETA: pPS_k = Σ(genotype × BETA × W_k)
  for (k in cluster_cols) {
    dt[, (k) := get(k) * BETA]
  }

  # Add total GRS column (unweighted BETA across all variants)
  dt[, grs_total := BETA]
  all_score_cols <- c(cluster_cols, "grs_total")
  cat(sprintf("  Score columns (BETA-weighted): %s\n", paste(all_score_cols, collapse = ", ")))

  # Write ONE multi-column score file for all clusters + total GRS
  score_path <- file.path(params$out_dir, "score_files",
    sprintf("score_%s.tsv", tolower(anc)))
  n_variants <- write_multi_score_file(dt, all_score_cols, score_path)

  if (n_variants == 0) {
    cat("  No variants with non-zero weights, skipping ancestry.\n")
    next
  }

  # Run plink2 (one call per chromosome, scoring all clusters + total GRS)
  out_prefix <- file.path(params$out_dir, "output",
    sprintf("prs_%s", tolower(anc)))
  cat(sprintf("  Running plink2 --score across 22 chromosomes (%d score columns per call)...\n",
    length(all_score_cols)))
  sscore_files <- run_plink2_score(
    score_path, length(all_score_cols),
    params$geno_dir, params$geno_prefix,
    out_prefix, params$plink2_path
  )

  # Merge chromosome results
  merged <- merge_sscore_files(sscore_files)
  if (!is.null(merged)) {
    id_col <- ifelse("#IID" %in% names(merged), "#IID",
              ifelse("IID" %in% names(merged), "IID", names(merged)[1]))
    score_cols <- grep("SCORE.*SUM", names(merged), value = TRUE)

    # Map SCORE{i}_SUM → prs_k{i}_{ancestry} + grs_total_{ancestry}
    prs_names <- c(
      sprintf("prs_%s_%s", tolower(cluster_cols), tolower(anc)),
      sprintf("grs_total_%s", tolower(anc))
    )

    if (length(score_cols) == length(prs_names)) {
      prs_dt <- merged[, c(id_col, score_cols), with = FALSE]
      setnames(prs_dt, c("person_id", prs_names))
      all_prs[[anc]] <- prs_dt
      cat(sprintf("  PRS computed for %d individuals (%s)\n",
        nrow(prs_dt), paste(prs_names, collapse = ", ")))
    } else {
      warning(sprintf("Expected %d score columns but got %d",
        length(prs_names), length(score_cols)))
    }

    # Clean up per-chromosome files
    for (f in sscore_files) {
      file.remove(f)
      log_f <- sub("\\.sscore$", ".log", f)
      if (file.exists(log_f)) file.remove(log_f)
    }
  }
}

# ---------------------------------------------------------------------------
# Merge all PRS into single output
# ---------------------------------------------------------------------------
cat("\n--- Merging all PRS scores ---\n")

if (length(all_prs) == 0) {
  stop("No PRS scores computed. Check input files and genotype data.")
}

# Start with the first ancestry, then left-join the rest
final_prs <- all_prs[[1]]
if (length(all_prs) > 1) {
  for (i in 2:length(all_prs)) {
    final_prs <- merge(final_prs, all_prs[[i]], by = "person_id", all = TRUE)
  }
}

# Write output
out_file <- file.path(params$out_dir, "output", "prs_all_clusters.tsv")
fwrite(final_prs, out_file, sep = "\t")

cat(sprintf("\nFinal output: %s\n", out_file))
cat(sprintf("  Individuals: %d\n", nrow(final_prs)))
cat(sprintf("  Columns: %s\n", paste(names(final_prs), collapse = ", ")))
cat("\nHead of output:\n")
print(head(final_prs))

cat("\n============================================\n")
cat("PRS computation complete!\n")
cat("============================================\n")
