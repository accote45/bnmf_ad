#!/usr/bin/env Rscript
# prs_scoring_helpers.R
# Shared plink2 --score scorer for the B1 PRS builders. Extracted from
# 01_compute_cluster_prs.R so the cluster PRS (01_) and the genome-wide PRS
# (01b_) score identically — same per-chromosome plink2 call, the same .sscore
# reuse cache, and the same per-group merge-by-summing. No analysis logic lives
# here beyond scoring; callers build the per-chromosome score files.
#
# Source it from a builder script:
#   source("scripts/b1_analysis/prs_scoring_helpers.R")

library(data.table)

#' Per-variant strongest-signal BETA across a list of GWAS files.
#'
#' For each VAR_ID (restricted to `var_ids`), returns the BETA / Effect_Allele /
#' P_VALUE from the GWAS file where the variant is most significant (minimum
#' P_VALUE). This mirrors how bNMF resolves a single effect per variant when
#' several GWAS contribute to the same trait's variant universe (allele
#' alignment conflict_rule = "strongest"). Harmonized files share the
#' CHR_POS_REF_ALT VAR_ID convention with Effect_Allele = ALT, so BETAs are
#' already co-oriented and can be compared/looked up by VAR_ID directly.
#'
#' @param var_ids Character vector of VAR_IDs to keep.
#' @param file_list List/character vector of harmonized GWAS file paths.
#' @return data.table with VAR_ID, P_VALUE, BETA, Effect_Allele (one row per VAR_ID).
strongest_beta <- function(var_ids, file_list) {
  var_ids <- unique(var_ids)
  acc <- vector("list", length(file_list))
  for (i in seq_along(file_list)) {
    d <- fread(file_list[[i]],
               select = c("VAR_ID", "P_VALUE", "BETA", "Effect_Allele"))
    acc[[i]] <- d[VAR_ID %in% var_ids]
  }
  combined <- rbindlist(acc)
  if (nrow(combined) == 0) return(combined)
  # Keep the most significant record per VAR_ID (NA P sorts last).
  setorder(combined, VAR_ID, P_VALUE, na.last = TRUE)
  combined[, .SD[1L], by = VAR_ID]
}

# --- Score cache helper ---------------------------------------------------
# A per-chromosome .sscore can be reused iff the inputs that determine it are
# unchanged: the score file (encodes A1 + weights), the sample keep file, and
# the genotype .bed. Key on md5 of the (small) score + keep files plus
# size+mtime of the (multi-GB) .bed — md5-ing the .bed would be too slow.
score_cache_key <- function(score_path, keep_path, bed_path) {
  bed_info <- file.info(bed_path)
  paste(tools::md5sum(score_path),
        tools::md5sum(keep_path),
        bed_info$size,
        as.integer(bed_info$mtime),
        sep = "|")
}

#' Score per ancestry group across chromosomes and merge by summing.
#'
#' @param chr_score_paths Named list, chr (as character) -> per-chr score file.
#'   Each score file must have columns: ID (variant id matching the .bim),
#'   A1 (effect allele), then one or more weight columns.
#' @param n_score_cols Number of weight columns in each score file.
#' @param out_names Character vector naming the summed score columns, in the
#'   same order as the weight columns (length == n_score_cols).
#' @param keep_files Named list, group name -> plink keep file (FID IID).
#' @param geno_prefix Genotype prefix; per-chr bed read as
#'   <dirname(geno_prefix)>/chr<chr>.bed (matches 01_compute_cluster_prs.R).
#' @param plink2_bin Path to the plink2 binary.
#' @param prs_dir Output directory for the .sscore + .cachekey files.
#' @return data.table with columns FID, IID, <out_names>, group (groups stacked).
score_by_group <- function(chr_score_paths, n_score_cols, out_names,
                           keep_files, geno_prefix, plink2_bin, prs_dir) {
  stopifnot(length(out_names) == n_score_cols)
  score_col_nums <- paste(seq(3, 2 + n_score_cols), collapse = ",")
  all_prs_list <- list()

  for (group_name in names(keep_files)) {
    keep_path <- keep_files[[group_name]]
    if (!file.exists(keep_path)) {
      cat(sprintf("  SKIP %s: keep file not found\n", group_name))
      next
    }

    cat(sprintf("\n  === %s ===\n", group_name))

    chr_sscore_files <- c()

    for (chr in names(chr_score_paths)) {
      chr_score_path <- chr_score_paths[[chr]]
      chr_geno      <- sprintf("%s/chr%s", dirname(geno_prefix), chr)
      out_prefix    <- file.path(prs_dir, sprintf("%s_chr%s", group_name, chr))
      sscore_path   <- paste0(out_prefix, ".sscore")
      cachekey_path <- paste0(out_prefix, ".sscore.cachekey")
      bed_path      <- paste0(chr_geno, ".bed")

      # Reuse a cached .sscore when the score file, keep file, and genotypes are
      # all unchanged since it was produced; otherwise (re)run plink2 and refresh
      # the cache key.
      cur_key <- score_cache_key(chr_score_path, keep_path, bed_path)
      cached  <- file.exists(sscore_path) && file.exists(cachekey_path) &&
                 identical(readLines(cachekey_path, warn = FALSE)[1], cur_key)

      if (cached) {
        cat(sprintf("    chr%s: cached (skipping plink2)\n", chr))
      } else {
        cmd <- sprintf(
          "OMP_NUM_THREADS=1 %s --bfile %s --keep %s --score %s 1 2 header cols=+scoresums ignore-dup-ids --score-col-nums %s --threads 1 --out %s 2>&1",
          plink2_bin, chr_geno, keep_path, chr_score_path, score_col_nums, out_prefix
        )
        ret <- system(cmd, intern = TRUE)
        if (file.exists(sscore_path)) writeLines(cur_key, cachekey_path)
      }

      if (file.exists(sscore_path)) {
        chr_sscore_files <- c(chr_sscore_files, sscore_path)
      } else {
        cat(sprintf("    WARNING: chr%s sscore not found\n", chr))
      }
    }

    cat(sprintf("    %d chromosome sscore files\n", length(chr_sscore_files)))

    if (length(chr_sscore_files) == 0) next

    # Merge per-chromosome scores by summing, joined by FID/IID
    group_dt <- NULL
    for (f in chr_sscore_files) {
      dt <- fread(f)
      setnames(dt, "#FID", "FID", skip_absent = TRUE)
      score_sum_cols <- grep("^SCORE.*_SUM$", colnames(dt), value = TRUE)

      chr_dt <- dt[, c("FID", "IID", score_sum_cols), with = FALSE]

      if (is.null(group_dt)) {
        group_dt <- chr_dt
      } else {
        group_dt <- merge(group_dt, chr_dt, by = c("FID", "IID"), all = TRUE,
                          suffixes = c("", ".new"))
        for (sc in score_sum_cols) {
          new_col <- paste0(sc, ".new")
          if (new_col %in% colnames(group_dt)) {
            group_dt[[sc]] <- fifelse(is.na(group_dt[[sc]]), 0, group_dt[[sc]]) +
                              fifelse(is.na(group_dt[[new_col]]), 0, group_dt[[new_col]])
            group_dt[[new_col]] <- NULL
          }
        }
      }
    }

    # Rename score columns (in plink's SCORE1_SUM, SCORE2_SUM, ... order) to
    # the caller's output names.
    score_sum_cols <- grep("^SCORE.*_SUM$", colnames(group_dt), value = TRUE)
    setnames(group_dt, score_sum_cols, out_names)

    group_dt[, group := group_name]
    all_prs_list[[group_name]] <- group_dt

    cat(sprintf("    %d individuals scored\n", nrow(group_dt)))
  }

  rbindlist(all_prs_list)
}

#' Map variants to per-chromosome imputed .bim files by base-pair position and
#' write one score file per chromosome.
#'
#' Mirrors the CHR:POS matching in 01_compute_cluster_prs.R: tries
#' <dirname(geno_prefix)>/chr<chr>.bim first, then the {chr}-substituted prefix.
#'
#' @param vars data.table with columns v_chr, v_pos, Effect_Allele, and the
#'   weight columns named in score_cols.
#' @param score_cols Character vector of weight column names to write.
#' @param geno_prefix Genotype prefix (may contain {chr}).
#' @param prs_dir Output directory.
#' @param tag Filename tag, e.g. "score" -> score_chr<chr>.tsv.
#' @return Named list, chr (as character) -> score file path.
write_chr_score_files <- function(vars, score_cols, geno_prefix, prs_dir,
                                  tag = "score") {
  chr_score_paths <- list()
  total_matched <- 0

  for (chr in 1:22) {
    chr_vars <- vars[v_chr == chr]
    if (nrow(chr_vars) == 0) next

    bim_path <- sprintf("%s/chr%d.bim", dirname(geno_prefix), chr)
    if (!file.exists(bim_path)) {
      bim_path <- paste0(gsub("\\{chr\\}", chr, geno_prefix), ".bim")
    }
    bim <- fread(bim_path,
                 col.names = c("CHR", "SNP", "CM", "POS", "A1_bim", "A2_bim"))

    chr_vars[, pos_key := v_pos]
    bim[, pos_key := POS]
    matched <- merge(chr_vars, bim[, .(pos_key, SNP)], by = "pos_key")

    if (nrow(matched) == 0) {
      cat(sprintf("  chr%d: 0/%d matched\n", chr, nrow(chr_vars)))
      next
    }

    total_matched <- total_matched + nrow(matched)

    score_dt <- matched[, c("SNP", "Effect_Allele", score_cols), with = FALSE]
    setnames(score_dt, c("SNP", "Effect_Allele"), c("ID", "A1"))

    chr_score_path <- file.path(prs_dir, sprintf("%s_chr%d.tsv", tag, chr))
    fwrite(score_dt, chr_score_path, sep = "\t")
    chr_score_paths[[as.character(chr)]] <- chr_score_path

    cat(sprintf("  chr%d: %d/%d matched\n", chr, nrow(matched), nrow(chr_vars)))
  }

  cat(sprintf("\nTotal matched: %d variants across %d chromosomes\n",
              total_matched, length(chr_score_paths)))
  chr_score_paths
}
