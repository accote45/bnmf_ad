# prep_bnmf.R
# GWAS preprocessing for bNMF analysis.
# Applies comprehensive QC: p-value, MAF, SNP-only, strand-ambiguous,
# multiallelic, MHC exclusion, and LD clumping (plink1.9).

library(data.table)


#' Validate GWAS summary statistics file format
#'
#' Checks that a file has the expected columns, correct types, and valid VAR_IDs.
#'
#' @param file_path Path to GWAS sumstats file (.txt, .txt.gz, .tsv.gz)
#' @param required_cols Character vector of required column names
#' @return List with $valid (logical), $messages (character vector of issues)
validate_gwas_format <- function(file_path,
                                  required_cols = c("VAR_ID", "P_VALUE", "BETA", "SE"),
                                  optional_cols = c("N", "MAF")) {
  messages <- character()
  warnings <- character()
  valid <- TRUE

  # Check file exists
  if (!file.exists(file_path)) {
    return(list(valid = FALSE, messages = paste("File not found:", file_path)))
  }

  # Check file is not empty
  info <- file.info(file_path)
  if (info$size == 0) {
    return(list(valid = FALSE, messages = paste("File is empty:", file_path)))
  }

  # Read header
  dt <- tryCatch(
    fread(file_path, nrows = 10),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) {
    return(list(valid = FALSE, messages = paste("Could not read file:", file_path)))
  }

  # Check required columns
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    messages <- c(messages, paste("Missing columns:", paste(missing_cols, collapse = ", ")))
    valid <- FALSE
  }

  # Check VAR_ID format (CHR_POS_REF_ALT; alleles may use D/I for indels)
  if ("VAR_ID" %in% names(dt)) {
    var_id_pattern <- "^[0-9]{1,2}_[0-9]+_[ATCGDI]+_[ATCGDI]+$"
    sample_ids <- head(dt$VAR_ID, 10)
    bad_ids <- sample_ids[!grepl(var_id_pattern, sample_ids)]
    if (length(bad_ids) > 0) {
      messages <- c(messages, paste("VAR_ID format invalid (expected CHR_POS_REF_ALT):",
                                    paste(head(bad_ids, 3), collapse = ", ")))
      valid <- FALSE
    }
  }

  # Check numeric columns (skip all-NA columns — fread reads them as logical)
  numeric_cols <- intersect(c("P_VALUE", "BETA", "SE", "N", "MAF"), names(dt))
  for (col in numeric_cols) {
    if (all(is.na(dt[[col]]))) next
    if (!is.numeric(dt[[col]])) {
      if (col %in% optional_cols) {
        warnings <- c(warnings, paste("Column", col, "is not numeric"))
      } else {
        messages <- c(messages, paste("Column", col, "is not numeric"))
        valid <- FALSE
      }
    }
  }

  # Check for all-NA required columns (hard fail)
  for (col in intersect(required_cols, names(dt))) {
    if (all(is.na(dt[[col]]))) {
      messages <- c(messages, paste("Column", col, "is entirely NA"))
      valid <- FALSE
    }
  }

  # Check for all-NA optional columns (warn only, no fail)
  for (col in intersect(optional_cols, names(dt))) {
    if (all(is.na(dt[[col]]))) {
      warnings <- c(warnings, paste("Column", col,
                                     "is entirely NA (sample-size normalization will be skipped)"))
    }
  }

  if (length(warnings) > 0) {
    cat(sprintf("  WARNINGS: %s\n", paste(warnings, collapse = "; ")))
  }

  if (valid && length(messages) == 0) messages <- "All checks passed"
  list(valid = valid, messages = messages)
}


#' Read GWAS sumstats file with auto-detection of separator
#'
#' @param file_path Path to GWAS file
#' @return data.table
read_gwas <- function(file_path) {
  fread(file_path, data.table = TRUE)
}


#' Run LD clumping via plink1.9
#'
#' Maps GWAS variants to reference panel RSIDs by CHR:POS (handles missing RSIDs),
#' then runs plink --clump per chromosome.
#'
#' @param dt data.table with CHR, POS, P_VALUE columns (parsed from VAR_ID)
#' @param ref_panel_prefix Path prefix for per-chromosome plink files
#'   (e.g., "reference/1kg_eur/1000G.EUR.QC" -> expects .1.bed, .2.bed, ...)
#' @param clump_r2 r-squared threshold for clumping (default 0.05)
#' @param clump_kb Clumping window in kb (default 500)
#' @param plink_bin Path to plink binary (default: uses module-loaded plink)
#' @return Character vector of VAR_IDs that survived clumping
# Ensure plink binary is available, trying module load and common HPC paths
# Returns resolved plink path, or NULL if not found
ensure_plink <- function(plink_bin = "plink") {
  if (Sys.which(plink_bin) != "") return(plink_bin)

  cat(sprintf("    plink not found in PATH, attempting: module load plink/1.90b6.21\n"))
  system("module load plink/1.90b6.21 2>/dev/null", intern = TRUE)
  modulecmd <- system("which modulecmd 2>/dev/null", intern = TRUE)
  if (length(modulecmd) > 0 && nchar(modulecmd[1]) > 0) {
    env_cmds <- system("modulecmd bash load plink/1.90b6.21 2>/dev/null", intern = TRUE)
    for (cmd in env_cmds) try(eval(parse(text = cmd)), silent = TRUE)
  }
  if (Sys.which(plink_bin) != "") return(plink_bin)

  # Fallback: try common HPC paths
  common_paths <- c(
    "/hpc/packages/minerva-centos7/plink/1.90b6.21/plink",
    "/opt/hpc/packages/minerva-centos7/plink/1.90b6.21/plink"
  )
  for (p in common_paths) {
    if (file.exists(p)) {
      cat(sprintf("    Using plink: %s\n", p))
      return(p)
    }
  }

  cat("    WARNING: plink not found.\n")
  return(NULL)
}

# Read BIM files from reference panel and return chr_pos -> RSID mapping
# Reusable by ld_clump() and find_proxy_variants()
read_ref_bim <- function(ref_panel_prefix) {
  rbindlist(lapply(1:22, function(chr) {
    bim_file <- paste0(ref_panel_prefix, ".", chr, ".bim")
    if (!file.exists(bim_file)) return(data.table())
    bim <- fread(bim_file, header = FALSE,
                 col.names = c("CHR", "RSID", "CM", "POS", "A1", "A2"))
    bim[, chr_pos := paste(CHR, POS, sep = ":")]
    bim
  }))
}

ld_clump <- function(dt, ref_panel_prefix, clump_r2 = 0.05, clump_kb = 500,
                     plink_bin = "plink") {
  plink_bin <- ensure_plink(plink_bin)
  if (is.null(plink_bin)) {
    cat("    WARNING: plink not found. LD clumping will be skipped.\n")
    return(dt$VAR_ID)
  }

  tmp_dir <- tempdir()
  on.exit({
    unlink(file.path(tmp_dir, "clump_input.tsv"))
    unlink(list.files(tmp_dir, pattern = "^clump_chr", full.names = TRUE))
  }, add = TRUE)

  # Build a CHR:POS -> RSID mapping from reference BIM files, and
  # a CHR:POS -> VAR_ID mapping from GWAS data
  gwas_map <- data.table(chr_pos = paste(dt$CHR, dt$POS, sep = ":"),
                         VAR_ID = dt$VAR_ID, P = dt$P_VALUE)

  # Read all reference BIM files to get RSIDs for position-based matching
  ref_snps <- read_ref_bim(ref_panel_prefix)[, .(chr_pos, RSID)]

  # Match GWAS variants to reference RSIDs by position
  matched <- merge(gwas_map, ref_snps, by = "chr_pos", all.x = FALSE)
  matched <- matched[!duplicated(chr_pos)]  # one RSID per position
  cat(sprintf("    LD clump: %d/%d GWAS variants matched to reference by position\n",
              nrow(matched), nrow(gwas_map)))

  if (nrow(matched) == 0) {
    cat("    WARNING: No variants matched reference panel. Skipping LD clumping.\n")
    return(dt$VAR_ID)  # return all if no matching possible
  }

  # Write sumstats for plink (SNP = reference RSID, P = GWAS p-value)
  sumstats_file <- file.path(tmp_dir, "clump_input.tsv")
  fwrite(matched[, .(SNP = RSID, P = P)], sumstats_file, sep = "\t")

  # Run plink --clump per chromosome
  clumped_rsids <- character()
  for (chr in 1:22) {
    bfile <- paste0(ref_panel_prefix, ".", chr)
    if (!file.exists(paste0(bfile, ".bed"))) next

    out_prefix <- file.path(tmp_dir, sprintf("clump_chr%d", chr))
    cmd <- sprintf(
      "%s --bfile %s --clump %s --clump-p1 1 --clump-p2 1 --clump-r2 %g --clump-kb %d --out %s --noweb 2>&1",
      plink_bin, bfile, sumstats_file, clump_r2, clump_kb, out_prefix
    )
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    clumped_file <- paste0(out_prefix, ".clumped")
    if (file.exists(clumped_file)) {
      clumped <- tryCatch(fread(clumped_file), error = function(e) NULL)
      if (!is.null(clumped) && nrow(clumped) > 0) {
        clumped_rsids <- c(clumped_rsids, trimws(clumped$SNP))
      }
    }

    # Clean up per-chromosome temp files
    unlink(list.files(tmp_dir, pattern = sprintf("clump_chr%d", chr), full.names = TRUE))
  }
  unlink(sumstats_file)

  # Map clumped RSIDs back to VAR_IDs
  clumped_varids <- matched[RSID %in% clumped_rsids, VAR_ID]

  cat(sprintf("    LD clump output: %d independent SNPs (r2 < %g, %dkb window)\n",
              length(clumped_varids), clump_r2, clump_kb))
  clumped_varids
}


#' Comprehensive variant QC for a single GWAS file
#'
#' Applies sequential QC filters with per-step reporting:
#'   1. P-value threshold
#'   2. MAF threshold
#'   3. SNPs only (single-nucleotide alleles)
#'   4. Strand-ambiguous removal
#'   5. Multiallelic removal (duplicate CHR_POS)
#'   6. MHC region exclusion (chr6:25-35Mb, GRCh37)
#'   7. LD clumping (plink1.9, ancestry-matched reference)
#'
#' @param gwas_file Path to harmonized GWAS sumstats
#' @param p_threshold P-value significance threshold
#' @param maf_threshold Minimum MAF
#' @param ref_panel_prefix Path prefix for per-chromosome plink reference files
#' @param clump_r2 r-squared threshold for LD clumping
#' @param clump_kb Clumping window in kb
#' @param plink_bin Path to plink binary
#' @return List with $data (filtered data.table) and $qc_report (data.frame of step counts)
qc_variants <- function(gwas_file, p_threshold = 5e-8, maf_threshold = 0.001,
                        ref_panel_prefix = NULL, clump_r2 = 0.05, clump_kb = 500,
                        plink_bin = "plink", hapmap3_file = NULL) {

  gwas_label <- sub("\\.gz$", "", basename(gwas_file))
  qc_log <- list()  # step_name -> variants_remaining

  cat(sprintf("\n  QC for: %s\n", gwas_label))

  # Step 0: Load
  dt <- read_gwas(gwas_file)
  n_total <- nrow(dt)
  qc_log[["0_total"]] <- n_total
  cat(sprintf("    [0] Total variants loaded: %d\n", n_total))

  # Step 1: P-value filter
  dt <- dt[P_VALUE < p_threshold]
  qc_log[["1_pvalue"]] <- nrow(dt)
  cat(sprintf("    [1] After P < %.0e: %d (removed %d)\n",
              p_threshold, nrow(dt), n_total - nrow(dt)))

  # Step 2: MAF filter (skip if MAF column is absent or all NA)
  n_before <- nrow(dt)
  if ("MAF" %in% names(dt) && any(!is.na(dt$MAF))) {
    dt <- dt[is.na(MAF) | MAF >= maf_threshold]
  }
  qc_log[["2_maf"]] <- nrow(dt)
  cat(sprintf("    [2] After MAF >= %.4f: %d (removed %d)\n",
              maf_threshold, nrow(dt), n_before - nrow(dt)))

  # Early return if no variants survive filtering

  if (nrow(dt) == 0L) {
    cat("    No variants remain after filtering — skipping this GWAS.\n")
    qc_report <- data.frame(
      gwas_file = gwas_label,
      step = names(qc_log),
      variants_remaining = unlist(qc_log),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    qc_report$variants_removed <- c(0, -diff(qc_report$variants_remaining))
    return(list(data = dt, qc_report = qc_report))
  }

  # Parse alleles from VAR_ID for subsequent steps
  parts <- tstrsplit(dt$VAR_ID, "_", fixed = TRUE)
  dt[, c("CHR", "POS", "A1", "A2") := .(as.integer(parts[[1]]),
                                          as.integer(parts[[2]]),
                                          parts[[3]], parts[[4]])]

  # Step 3: SNPs only (single-nucleotide alleles)
  n_before <- nrow(dt)
  dt <- dt[nchar(A1) == 1L & nchar(A2) == 1L]
  qc_log[["3_snps_only"]] <- nrow(dt)
  cat(sprintf("    [3] After SNPs only (no indels): %d (removed %d)\n",
              nrow(dt), n_before - nrow(dt)))

  # Step 4: Strand-ambiguous removal
  n_before <- nrow(dt)
  allele_pairs <- paste0(dt$A1, dt$A2)
  ambiguous <- allele_pairs %in% c("AT", "TA", "CG", "GC")
  dt <- dt[!ambiguous]
  qc_log[["4_strand_ambiguous"]] <- nrow(dt)
  cat(sprintf("    [4] After removing strand-ambiguous: %d (removed %d)\n",
              nrow(dt), n_before - nrow(dt)))

  # Step 5: Multiallelic removal (all variants at positions with >1 entry)
  n_before <- nrow(dt)
  dt[, chr_pos := paste(CHR, POS, sep = "_")]
  dup_pos <- dt$chr_pos[duplicated(dt$chr_pos)]
  dt <- dt[!chr_pos %in% dup_pos]
  qc_log[["5_multiallelic"]] <- nrow(dt)
  cat(sprintf("    [5] After removing multiallelic: %d (removed %d)\n",
              nrow(dt), n_before - nrow(dt)))

  # Step 5.5: HapMap3 restriction (optional — reduces ancestry-specific SNP discovery bias)
  n_before <- nrow(dt)
  if (!is.null(hapmap3_file) && file.exists(hapmap3_file)) {
    hm3 <- fread(hapmap3_file, select = "VAR_ID")
    dt <- dt[VAR_ID %in% hm3$VAR_ID]
    qc_log[["5.5_hapmap3"]] <- nrow(dt)
    cat(sprintf("    [5.5] After HapMap3 restriction: %d (removed %d)\n",
                nrow(dt), n_before - nrow(dt)))
  } else {
    qc_log[["5.5_hapmap3"]] <- nrow(dt)
    if (!is.null(hapmap3_file)) {
      cat(sprintf("    [5.5] HapMap3: SKIPPED (file not found: %s)\n", hapmap3_file))
    } else {
      cat("    [5.5] HapMap3: SKIPPED (not configured)\n")
    }
  }

  # Step 6: MHC exclusion (chr6:25,000,000-35,000,000 in GRCh37)
  n_before <- nrow(dt)
  in_mhc <- dt$CHR == 6L & dt$POS >= 25000000L & dt$POS <= 35000000L
  dt <- dt[!in_mhc]
  qc_log[["6_mhc_exclusion"]] <- nrow(dt)
  cat(sprintf("    [6] After MHC exclusion (chr6:25-35Mb): %d (removed %d)\n",
              nrow(dt), n_before - nrow(dt)))

  # Step 7: LD clumping
  n_before <- nrow(dt)
  if (!is.null(ref_panel_prefix)) {
    clumped_varids <- ld_clump(dt, ref_panel_prefix,
                               clump_r2 = clump_r2, clump_kb = clump_kb,
                               plink_bin = plink_bin)
    dt <- dt[VAR_ID %in% clumped_varids]
    qc_log[["7_ld_clump"]] <- nrow(dt)
    cat(sprintf("    [7] After LD clumping (r2<%g, %dkb): %d (removed %d)\n",
                clump_r2, clump_kb, nrow(dt), n_before - nrow(dt)))
  } else {
    qc_log[["7_ld_clump"]] <- nrow(dt)
    cat("    [7] LD clumping: SKIPPED (no reference panel provided)\n")
  }

  # Clean up temporary columns
  dt[, c("CHR", "POS", "A1", "A2", "chr_pos") := NULL]

  cat(sprintf("    Final: %d variants (%.2f%% of total)\n",
              nrow(dt), 100 * nrow(dt) / n_total))

  # Build QC report data.frame
  qc_report <- data.frame(
    gwas_file = gwas_label,
    step = names(qc_log),
    variants_remaining = unlist(qc_log),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  qc_report$variants_removed <- c(0, -diff(qc_report$variants_remaining))

  list(data = dt, qc_report = qc_report)
}


#' QC and union variants from multiple reference GWAS files
#'
#' Calls qc_variants() on each file, then takes the union of unique VAR_IDs.
#' Performs a post-union multiallelic check to remove any positions that
#' appear more than once across the combined set.
#'
#' @param gwas_files Character vector of file paths to reference GWAS files
#' @param p_threshold P-value significance threshold
#' @param maf_threshold Minimum MAF
#' @param ref_panel_prefix Path prefix for per-chromosome plink reference files
#' @param clump_r2 r-squared threshold for LD clumping
#' @param clump_kb Clumping window in kb
#' @param plink_bin Path to plink binary
#' @return List with $data (union data.table), $qc_report (combined QC report data.frame)
qc_variants_multi <- function(gwas_files, p_threshold = 5e-8, maf_threshold = 0.001,
                              ref_panel_prefix = NULL, clump_r2 = 0.05, clump_kb = 500,
                              plink_bin = "plink", hapmap3_file = NULL,
                              union_clump = FALSE, union_hapmap3_file = NULL) {
  all_variants <- data.table()
  all_reports <- data.frame()

  for (i in seq_along(gwas_files)) {
    cat(sprintf("\n--- Reference GWAS %d/%d ---\n", i, length(gwas_files)))
    result <- qc_variants(gwas_files[i], p_threshold = p_threshold,
                          maf_threshold = maf_threshold,
                          ref_panel_prefix = ref_panel_prefix,
                          clump_r2 = clump_r2, clump_kb = clump_kb,
                          plink_bin = plink_bin,
                          hapmap3_file = hapmap3_file)
    all_variants <- rbind(all_variants, result$data, fill = TRUE)
    all_reports <- rbind(all_reports, result$qc_report)
  }

  # Strongest (minimum) P-value per VAR_ID across all reference files, computed
  # from the combined set BEFORE dedup so it reflects every trait the variant
  # was significant in. Used below to rank index SNPs for the union LD clump.
  n_combined <- nrow(all_variants)
  minp <- if (n_combined > 0) {
    all_variants[, .(minP = min(P_VALUE, na.rm = TRUE)), by = VAR_ID]
  } else {
    data.table(VAR_ID = character(), minP = numeric())
  }

  if (n_combined == 0) {
    return(list(data = all_variants, qc_report = all_reports))
  }

  if (union_clump) {
    # ---- Lean union (union-clump ancestries, e.g. META) -------------------
    # Stack the per-file clumped sets, restrict to HapMap3, then LD clump
    # across the pooled set. No standalone VAR_ID dedup / multiallelic check:
    # the union clump (one index SNP per LD region, one RSID per position)
    # collapses duplicate VAR_IDs and multiallelic positions itself, so they
    # are redundant here.
    cat(sprintf("\nUnion (stacked) across %d reference files: %d variants\n",
                length(gwas_files), n_combined))
    union_steps     <- "union_combined"
    union_remaining <- n_combined
    union_removed   <- 0

    # HapMap3 restriction at the union level (before clumping)
    if (!is.null(union_hapmap3_file) && file.exists(union_hapmap3_file)) {
      n_before <- nrow(all_variants)
      hm3 <- fread(union_hapmap3_file, select = "VAR_ID")
      all_variants <- all_variants[VAR_ID %in% hm3$VAR_ID]
      cat(sprintf("Union HapMap3 restriction: %d -> %d (removed %d)\n",
                  n_before, nrow(all_variants), n_before - nrow(all_variants)))
      union_steps     <- c(union_steps, "union_hapmap3")
      union_remaining <- c(union_remaining, nrow(all_variants))
      union_removed   <- c(union_removed, n_before - nrow(all_variants))
    }

    # LD clump across the union (independence across all pooled traits).
    # Per-file clumping only guarantees independence WITHIN a trait.
    if (!is.null(ref_panel_prefix) && nrow(all_variants) > 0) {
      n_before <- nrow(all_variants)
      cat(sprintf("Union LD clump (r2<%g, %dkb): clumping pooled variants...\n",
                  clump_r2, clump_kb))
      # Re-derive CHR/POS (dropped by per-file QC) and rank index SNPs by the
      # strongest P across files (min P per VAR_ID) on a temp copy.
      clump_dt <- copy(all_variants)
      cp <- tstrsplit(clump_dt$VAR_ID, "_", fixed = TRUE)
      clump_dt[, CHR := as.integer(cp[[1]])]
      clump_dt[, POS := as.integer(cp[[2]])]
      clump_dt[minp, P_VALUE := i.minP, on = "VAR_ID"]

      clumped_varids <- ld_clump(clump_dt, ref_panel_prefix,
                                 clump_r2 = clump_r2, clump_kb = clump_kb,
                                 plink_bin = plink_bin)
      # Keep one row per surviving variant (stacking may leave VAR_ID dupes).
      all_variants <- all_variants[VAR_ID %in% clumped_varids][!duplicated(VAR_ID)]
      cat(sprintf("Union LD clump: %d -> %d (removed %d)\n",
                  n_before, nrow(all_variants), n_before - nrow(all_variants)))
      union_steps     <- c(union_steps, "union_ld_clump")
      union_remaining <- c(union_remaining, nrow(all_variants))
      union_removed   <- c(union_removed, n_before - nrow(all_variants))
    } else {
      cat("\nUnion LD clump: SKIPPED (no reference panel provided)\n")
      all_variants <- all_variants[!duplicated(VAR_ID)]  # still ensure unique variants
    }
    n_final <- nrow(all_variants)

  } else {
    # ---- Standard union (non-union-clump ancestries) ----------------------
    # Remove duplicate VAR_IDs, then drop multiallelic positions.
    all_variants <- all_variants[!duplicated(VAR_ID)]
    n_union <- nrow(all_variants)
    cat(sprintf("\nUnion across %d reference files: %d unique variants (from %d total)\n",
                length(gwas_files), n_union, n_combined))

    parts_union <- tstrsplit(all_variants$VAR_ID, "_", fixed = TRUE)
    all_variants[, chr_pos := paste(parts_union[[1]], parts_union[[2]], sep = "_")]
    dup_pos <- all_variants$chr_pos[duplicated(all_variants$chr_pos)]
    n_before_dedup <- nrow(all_variants)
    all_variants <- all_variants[!chr_pos %in% dup_pos]
    all_variants[, chr_pos := NULL]
    n_final <- nrow(all_variants)
    cat(sprintf("Post-union multiallelic check: %d -> %d (removed %d)\n",
                n_before_dedup, n_final, n_before_dedup - n_final))

    union_steps     <- c("union_combined", "union_dedup_varid", "union_multiallelic_qc")
    union_remaining <- c(n_combined, n_union, n_final)
    union_removed   <- c(0, n_combined - n_union, n_union - n_final)
  }

  union_report <- data.frame(
    gwas_file = "UNION",
    step = union_steps,
    variants_remaining = union_remaining,
    variants_removed = union_removed,
    stringsAsFactors = FALSE
  )
  all_reports <- rbind(all_reports, union_report)

  cat(sprintf("\nFinal variant set: %d variants\n", n_final))

  list(data = all_variants, qc_report = all_reports)
}


#' Build Z-score matrix from filtered variants and trait GWAS files
#'
#' @param filtered_variants data.table with VAR_ID, BETA, SE columns from CAD GWAS
#' @param trait_files Named list of file paths (names = trait names)
#' @param proxy_map data.table with target_VAR_ID, proxy_VAR_ID columns (NULL to skip)
#' @param impute_method "median" to fill remaining missing with per-trait median, "zero" for 0
#' @return List with z_matrix (variants x traits), n_matrix (sample sizes), proxy_usage (per-trait counts)
build_z_matrix <- function(filtered_variants, trait_files,
                           proxy_map = NULL, impute_method = "zero") {
  variant_ids <- filtered_variants$VAR_ID
  n_variants <- length(variant_ids)
  trait_names <- names(trait_files)
  n_traits <- length(trait_names)

  z_matrix <- matrix(0, nrow = n_variants, ncol = n_traits,
                     dimnames = list(variant_ids, trait_names))
  n_matrix <- matrix(NA_real_, nrow = n_variants, ncol = n_traits,
                     dimnames = list(variant_ids, trait_names))

  # Track proxy usage per trait
  proxy_counts <- setNames(integer(n_traits), trait_names)

  for (i in seq_along(trait_names)) {
    trait <- trait_names[i]
    trait_file <- trait_files[[trait]]
    cat(sprintf("Processing trait: %s (%s)\n", trait, basename(trait_file)))

    trait_dt <- tryCatch(read_gwas(trait_file), error = function(e) {
      warning(sprintf("Could not read trait file %s: %s", trait_file, e$message))
      return(NULL)
    })
    if (is.null(trait_dt)) next

    # Direct match: variants in the universe
    direct_match <- trait_dt[VAR_ID %in% variant_ids]
    matched <- direct_match$VAR_ID

    # Calculate z-scores for direct matches
    if (length(matched) > 0) {
      z_scores <- direct_match$BETA / direct_match$SE
      sample_sizes <- direct_match$N

      # Align effect alleles with reference GWAS
      if ("Effect_Allele" %in% names(direct_match) && "Effect_Allele" %in% names(filtered_variants)) {
        ref_alleles <- filtered_variants[match(matched, VAR_ID), Effect_Allele]
        trait_alleles <- direct_match$Effect_Allele
        flip <- ref_alleles != trait_alleles
        z_scores[flip] <- -z_scores[flip]
      }

      idx <- match(matched, variant_ids)
      z_matrix[idx, trait] <- z_scores
      n_matrix[idx, trait] <- sample_sizes
    }

    # Proxy lookup: for universe variants NOT matched directly, try proxy
    if (!is.null(proxy_map) && nrow(proxy_map) > 0) {
      unmatched_ids <- setdiff(variant_ids, matched)
      proxy_needed <- proxy_map[target_VAR_ID %in% unmatched_ids]

      if (nrow(proxy_needed) > 0) {
        # Look up proxy VAR_IDs in the trait data
        proxy_hits <- trait_dt[VAR_ID %in% proxy_needed$proxy_VAR_ID]
        if (nrow(proxy_hits) > 0) {
          # Merge to get target_VAR_ID for each proxy hit
          proxy_merged <- merge(proxy_needed[, .(target_VAR_ID, proxy_VAR_ID)],
                                proxy_hits, by.x = "proxy_VAR_ID", by.y = "VAR_ID")

          p_z <- proxy_merged$BETA / proxy_merged$SE
          p_n <- proxy_merged$N

          # Align proxy effect alleles with reference (use proxy's own alleles)
          # Proxy Z-scores are used as-is since they represent the proxy's own effect
          p_idx <- match(proxy_merged$target_VAR_ID, variant_ids)
          z_matrix[p_idx, trait] <- p_z
          n_matrix[p_idx, trait] <- p_n
          proxy_counts[trait] <- nrow(proxy_merged)
        }
      }
    }

    n_filled <- sum(z_matrix[, trait] != 0)
    cat(sprintf("  Matched %d/%d variants (%.1f%%), %d via proxy\n",
                n_filled, n_variants, 100 * n_filled / n_variants, proxy_counts[trait]))
  }

  # Normalize z-scores by sample size
  cat("Normalizing z-scores by sample size...\n")
  sqrt_n <- sqrt(n_matrix)
  mean_sqrt_n <- mean(sqrt_n, na.rm = TRUE)
  z_normalized <- z_matrix
  for (j in seq_len(ncol(z_matrix))) {
    valid <- !is.na(n_matrix[, j]) & n_matrix[, j] > 0
    if (any(valid)) {
      z_normalized[valid, j] <- z_matrix[valid, j] / sqrt(n_matrix[valid, j]) * mean_sqrt_n
    }
  }

  # Replace NAs with 0 first
  z_normalized[is.na(z_normalized)] <- 0

  # Impute remaining zeros with trait median (if configured)
  if (impute_method == "median") {
    cat("Imputing remaining missing values with per-trait median...\n")
    for (j in seq_len(ncol(z_normalized))) {
      nonzero <- z_normalized[, j] != 0
      if (any(nonzero) && !all(nonzero)) {
        trait_median <- median(z_normalized[nonzero, j])
        z_normalized[!nonzero, j] <- trait_median
        cat(sprintf("  %s: imputed %d values with median = %.4f\n",
                    trait_names[j], sum(!nonzero), trait_median))
      }
    }
  }

  list(z_matrix = z_normalized, n_matrix = n_matrix, proxy_usage = proxy_counts)
}


#' Align variant rows of the Z-score matrix to the disease-risk-raising allele
#'
#' The Z-matrix rows are oriented to filtered$Effect_Allele (an arbitrary
#' harmonized allele). Because expand_to_nonneg() routes z>0 -> trait_pos and
#' z<0 -> trait_neg, variants that are biologically identical but carry opposite
#' allele coding get split into mirror clusters. This flips each variant row by a
#' single +/-1 sign so its effect allele is the CAD/T2D-risk-raising allele,
#' collapsing those mirrors. A per-variant scalar preserves the cross-trait
#' pattern and commutes with the (linear) sample-size normalization.
#'
#' @param z_matrix variants x traits matrix, rownames = VAR_ID (oriented to
#'   filtered$Effect_Allele)
#' @param filtered QC'd universe data.table (VAR_ID, Effect_Allele, ...)
#' @param ref_gwas_named_paths NAMED character vector of disease GWAS paths;
#'   names prefixed CAD_ / T2D_ so the two diseases are separable
#' @param conflict_rule "strongest" (orient to disease with smallest p per
#'   variant), "priority_T2D", "priority_CAD", or "first_occurrence"
#'   (sign(filtered$BETA); no file re-reading)
#' @return list(z_matrix = aligned, report = data.table of counts)
align_z_to_disease_risk <- function(z_matrix, filtered, ref_gwas_named_paths,
                                    conflict_rule = "strongest") {
  var_ids <- rownames(z_matrix)
  fmatch  <- match(var_ids, filtered$VAR_ID)
  ref_ea  <- filtered$Effect_Allele[fmatch]   # orientation the matrix already uses

  n_cad_aligned <- 0L; n_t2d_aligned <- 0L; n_discordant <- 0L; n_no_record <- 0L

  # Cheapest rule: filtered already carries the first-occurrence disease BETA
  # oriented to filtered$Effect_Allele (config lists CAD first -> CAD-prioritized).
  if (conflict_rule == "first_occurrence") {
    beta_ref <- filtered$BETA[fmatch]
    sgn <- sign(beta_ref)
    n_no_record <- sum(is.na(beta_ref))
  } else {
    keys <- names(ref_gwas_named_paths)
    if (is.null(keys)) keys <- rep("", length(ref_gwas_named_paths))
    cad_paths <- ref_gwas_named_paths[startsWith(keys, "CAD_") | keys == "CAD"]
    t2d_paths <- ref_gwas_named_paths[startsWith(keys, "T2D_") | keys == "T2D"]

    # Graceful degradation: if we cannot separate the two diseases, fall back.
    if (length(cad_paths) == 0 || length(t2d_paths) == 0) {
      cat("  WARNING: could not split CAD/T2D ref GWAS by key prefix; ",
          "falling back to first_occurrence.\n", sep = "")
      return(align_z_to_disease_risk(z_matrix, filtered, ref_gwas_named_paths,
                                     conflict_rule = "first_occurrence"))
    }

    # For one disease group: per universe variant, oriented BETA + P of the
    # most-significant (min-P) record across that group's files.
    best_for_group <- function(paths) {
      best_p <- rep(NA_real_, length(var_ids))
      best_b <- rep(NA_real_, length(var_ids))
      for (p in paths) {
        d <- tryCatch(fread(p, select = c("VAR_ID", "Effect_Allele", "BETA", "P_VALUE")),
                      error = function(e) NULL)
        if (is.null(d) || nrow(d) == 0) next
        d <- d[VAR_ID %in% var_ids]
        if (nrow(d) == 0) next
        m <- match(d$VAR_ID, var_ids)
        flip <- d$Effect_Allele != ref_ea[m]            # orient to the matrix allele
        bo <- d$BETA; bo[flip] <- -bo[flip]
        better <- is.na(best_p[m]) | (d$P_VALUE < best_p[m])
        upd <- m[better]
        best_p[upd] <- d$P_VALUE[better]
        best_b[upd] <- bo[better]
      }
      list(p = best_p, b = best_b)
    }

    cad <- best_for_group(cad_paths)
    t2d <- best_for_group(t2d_paths)

    has_cad <- !is.na(cad$p); has_t2d <- !is.na(t2d$p)
    # discordance is diagnostic only (sign of best CAD vs best T2D)
    both <- has_cad & has_t2d
    n_discordant <- sum(both & (sign(cad$b) != sign(t2d$b)), na.rm = TRUE)

    # choose which disease's oriented BETA defines the sign, per rule
    use_cad <- rep(FALSE, length(var_ids))
    if (conflict_rule == "priority_T2D") {
      use_cad <- has_cad & !has_t2d
    } else if (conflict_rule == "priority_CAD") {
      use_cad <- has_cad
    } else {  # "strongest": smaller p wins; ties -> CAD
      use_cad <- has_cad & (!has_t2d | (cad$p <= t2d$p))
    }
    chosen_b <- ifelse(use_cad, cad$b, t2d$b)
    # variants present in only one group still use that group regardless of rule
    chosen_b[!has_cad & has_t2d] <- t2d$b[!has_cad & has_t2d]
    chosen_b[has_cad & !has_t2d] <- cad$b[has_cad & !has_t2d]

    sgn <- sign(chosen_b)
    n_cad_aligned <- sum(use_cad & has_cad, na.rm = TRUE)
    n_t2d_aligned <- sum(!use_cad & has_t2d, na.rm = TRUE)
    n_no_record   <- sum(!has_cad & !has_t2d)
  }

  # sign 0 / NA (no usable disease record) -> keep variant as-is (+1)
  sgn[is.na(sgn) | sgn == 0] <- 1
  z_aligned <- sweep(z_matrix, 1, sgn, "*")

  n_flipped <- sum(sgn < 0)
  cat(sprintf("  Allele alignment (%s): %d/%d variants flipped; CAD-aligned=%d, T2D-aligned=%d, discordant=%d, no-record=%d\n",
              conflict_rule, n_flipped, length(var_ids),
              n_cad_aligned, n_t2d_aligned, n_discordant, n_no_record))

  report <- data.table(
    conflict_rule = conflict_rule,
    n_total = length(var_ids),
    n_flipped = n_flipped,
    n_cad_aligned = n_cad_aligned,
    n_t2d_aligned = n_t2d_aligned,
    n_discordant = n_discordant,
    n_no_record = n_no_record
  )
  list(z_matrix = z_aligned, report = report)
}


# =============================================================================
# Trait missingness, proxy variants, and trait filtering (Smith et al. 2024)
# =============================================================================

#' Calculate per-variant trait missingness
#'
#' For each variant in the universe, count how many trait GWAS files contain it.
#' Missingness = n_missing / n_total_traits.
#'
#' @param variant_ids Character vector of VAR_IDs in the variant universe
#' @param trait_files Named list of file paths to trait GWAS files
#' @return data.table with columns: VAR_ID, n_present, n_missing, missingness
calculate_trait_missingness <- function(variant_ids, trait_files) {
  n_traits <- length(trait_files)
  trait_names <- names(trait_files)

  # Build presence matrix: for each trait, which variants are present
  presence <- matrix(FALSE, nrow = length(variant_ids), ncol = n_traits,
                     dimnames = list(variant_ids, trait_names))

  for (i in seq_along(trait_files)) {
    trait_file <- trait_files[[i]]
    trait_varids <- tryCatch({
      fread(trait_file, select = "VAR_ID", data.table = TRUE)$VAR_ID
    }, error = function(e) {
      warning(sprintf("Could not read VAR_IDs from %s: %s", trait_file, e$message))
      character(0)
    })
    presence[variant_ids %in% trait_varids, i] <- TRUE
    cat(sprintf("  Trait missingness scan [%d/%d]: %s — %d/%d variants present\n",
                i, n_traits, trait_names[i],
                sum(presence[, i]), length(variant_ids)))
  }

  n_present <- rowSums(presence)
  data.table(
    VAR_ID = variant_ids,
    n_present = n_present,
    n_missing = n_traits - n_present,
    missingness = (n_traits - n_present) / n_traits
  )
}


#' Find LD proxy candidates for high-missingness variants
#'
#' Uses plink --ld-snp-list --r2 to find variants in LD, batched by chromosome.
#'
#' @param target_var_ids Character vector of VAR_IDs needing proxies
#' @param ref_panel_prefix Path prefix for per-chromosome plink files
#' @param r2_threshold Minimum r2 for proxy candidates (default 0.8)
#' @param ld_window_kb LD window in kb (default 1000)
#' @param plink_bin Path to plink binary
#' @return data.table: target_VAR_ID, proxy_VAR_ID, R2, proxy_A1, proxy_A2
find_proxy_variants <- function(target_var_ids, ref_panel_prefix,
                                r2_threshold = 0.8, ld_window_kb = 1000,
                                plink_bin = "plink") {
  plink_bin <- ensure_plink(plink_bin)
  if (is.null(plink_bin)) {
    cat("    WARNING: plink not found. Proxy search skipped.\n")
    return(data.table(target_VAR_ID = character(), proxy_VAR_ID = character(),
                      R2 = numeric(), proxy_A1 = character(), proxy_A2 = character()))
  }

  # Read reference BIM and build mappings
  ref_bim <- read_ref_bim(ref_panel_prefix)

  # Parse target VAR_IDs to CHR:POS for matching
  parts <- tstrsplit(target_var_ids, "_", fixed = TRUE)
  target_map <- data.table(
    VAR_ID = target_var_ids,
    chr_pos = paste(parts[[1]], parts[[2]], sep = ":"),
    CHR = as.integer(parts[[1]])
  )

  # Match targets to reference RSIDs by position
  target_matched <- merge(target_map, ref_bim[, .(chr_pos, RSID)],
                          by = "chr_pos", all.x = FALSE)
  target_matched <- target_matched[!duplicated(VAR_ID)]

  if (nrow(target_matched) == 0) {
    cat("    WARNING: No target variants matched reference panel for proxy search.\n")
    return(data.table(target_VAR_ID = character(), proxy_VAR_ID = character(),
                      R2 = numeric(), proxy_A1 = character(), proxy_A2 = character()))
  }

  cat(sprintf("    Proxy search: %d/%d targets matched to reference\n",
              nrow(target_matched), length(target_var_ids)))

  tmp_dir <- tempdir()
  all_proxies <- list()

  # Process per chromosome
  for (chr in sort(unique(target_matched$CHR))) {
    chr_targets <- target_matched[CHR == chr]
    if (nrow(chr_targets) == 0) next

    bfile <- paste0(ref_panel_prefix, ".", chr)
    if (!file.exists(paste0(bfile, ".bed"))) next

    # Write target RSIDs to file
    snp_file <- file.path(tmp_dir, sprintf("proxy_targets_chr%d.txt", chr))
    writeLines(chr_targets$RSID, snp_file)

    out_prefix <- file.path(tmp_dir, sprintf("proxy_chr%d", chr))
    cmd <- sprintf(
      "%s --bfile %s --ld-snp-list %s --r2 --ld-window-r2 %g --ld-window-kb %d --ld-window 99999 --out %s --noweb 2>&1",
      plink_bin, bfile, snp_file, r2_threshold, ld_window_kb, out_prefix
    )
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    ld_file <- paste0(out_prefix, ".ld")
    if (file.exists(ld_file)) {
      ld_result <- tryCatch(fread(ld_file), error = function(e) NULL)
      if (!is.null(ld_result) && nrow(ld_result) > 0) {
        # Columns: CHR_A, BP_A, SNP_A, CHR_B, BP_B, SNP_B, R2
        # SNP_A = target RSID, SNP_B = proxy RSID
        ld_result <- ld_result[SNP_A != SNP_B]  # remove self-matches
        if (nrow(ld_result) > 0) {
          all_proxies[[length(all_proxies) + 1]] <- ld_result
        }
      }
    }

    # Clean up temp files
    unlink(snp_file)
    unlink(list.files(tmp_dir, pattern = sprintf("proxy_chr%d", chr), full.names = TRUE))
  }

  if (length(all_proxies) == 0) {
    cat("    No proxy candidates found.\n")
    return(data.table(target_VAR_ID = character(), proxy_VAR_ID = character(),
                      R2 = numeric(), proxy_A1 = character(), proxy_A2 = character()))
  }

  ld_all <- rbindlist(all_proxies)

  # Map target RSIDs back to VAR_IDs
  target_rsid_map <- target_matched[, .(RSID, target_VAR_ID = VAR_ID)]
  ld_all <- merge(ld_all, target_rsid_map, by.x = "SNP_A", by.y = "RSID", all.x = TRUE)

  # Map proxy RSIDs to VAR_ID format using BIM alleles
  proxy_bim <- ref_bim[, .(RSID, proxy_CHR = CHR, proxy_POS = POS, proxy_A1 = A1, proxy_A2 = A2)]
  ld_all <- merge(ld_all, proxy_bim, by.x = "SNP_B", by.y = "RSID", all.x = TRUE)
  ld_all[, proxy_VAR_ID := paste(proxy_CHR, proxy_POS, proxy_A1, proxy_A2, sep = "_")]

  result <- ld_all[, .(target_VAR_ID, proxy_VAR_ID, R2, proxy_A1, proxy_A2)]
  result <- result[!duplicated(paste(target_VAR_ID, proxy_VAR_ID))]

  cat(sprintf("    Found %d proxy candidates for %d targets\n",
              nrow(result), uniqueN(result$target_VAR_ID)))
  result
}


#' Select best proxy for each target variant
#'
#' Filters proxy candidates by missingness, strand ambiguity, and multi-allelicity,
#' then selects the best proxy per target by: lowest missingness, then highest R2.
#'
#' @param proxy_dt data.table from find_proxy_variants()
#' @param variant_ids_universe Character vector of original universe VAR_IDs (to exclude)
#' @param trait_files Named list of trait file paths
#' @param missingness_threshold Max missingness fraction for proxies (default 0.20)
#' @return data.table: target_VAR_ID, proxy_VAR_ID, R2, proxy_missingness (one row per target)
select_best_proxy <- function(proxy_dt, variant_ids_universe, trait_files,
                              missingness_threshold = 0.20) {
  if (nrow(proxy_dt) == 0) return(proxy_dt[0])

  # Exclude proxies already in the variant universe
  proxy_dt <- proxy_dt[!proxy_VAR_ID %in% variant_ids_universe]
  if (nrow(proxy_dt) == 0) {
    cat("    All proxy candidates already in universe, none needed.\n")
    return(data.table(target_VAR_ID = character(), proxy_VAR_ID = character(),
                      R2 = numeric(), proxy_missingness = numeric()))
  }

  # Calculate missingness for unique proxy VAR_IDs
  unique_proxy_ids <- unique(proxy_dt$proxy_VAR_ID)
  cat(sprintf("    Calculating missingness for %d unique proxy candidates...\n",
              length(unique_proxy_ids)))
  proxy_miss <- calculate_trait_missingness(unique_proxy_ids, trait_files)

  # Annotate proxies with missingness
  proxy_dt <- merge(proxy_dt, proxy_miss[, .(VAR_ID, missingness)],
                    by.x = "proxy_VAR_ID", by.y = "VAR_ID", all.x = TRUE)
  setnames(proxy_dt, "missingness", "proxy_missingness")

  # Annotate: is strand-ambiguous?
  proxy_alleles <- paste0(proxy_dt$proxy_A1, proxy_dt$proxy_A2)
  proxy_dt[, is_ambiguous := proxy_alleles %in% c("AT", "TA", "CG", "GC")]

  # Annotate: is multiallelic? (duplicate CHR_POS among proxy candidates)
  proxy_parts <- tstrsplit(proxy_dt$proxy_VAR_ID, "_", fixed = TRUE)
  proxy_dt[, proxy_chr_pos := paste(proxy_parts[[1]], proxy_parts[[2]], sep = "_")]
  dup_pos <- proxy_dt$proxy_chr_pos[duplicated(proxy_dt$proxy_chr_pos)]
  proxy_dt[, is_multiallelic := proxy_chr_pos %in% dup_pos]

  # Filter: remove problematic proxies
  n_before <- nrow(proxy_dt)
  proxy_dt <- proxy_dt[proxy_missingness <= missingness_threshold &
                        is_ambiguous == FALSE &
                        is_multiallelic == FALSE]
  cat(sprintf("    Proxy candidates after filtering: %d/%d\n", nrow(proxy_dt), n_before))

  if (nrow(proxy_dt) == 0) {
    return(data.table(target_VAR_ID = character(), proxy_VAR_ID = character(),
                      R2 = numeric(), proxy_missingness = numeric()))
  }

  # Select best proxy per target: lowest missingness, then highest R2
  setorder(proxy_dt, target_VAR_ID, proxy_missingness, -R2)
  best <- proxy_dt[, .SD[1], by = target_VAR_ID]

  cat(sprintf("    Selected proxies for %d targets\n", nrow(best)))
  best[, .(target_VAR_ID, proxy_VAR_ID, R2, proxy_missingness)]
}


#' Filter traits from Z-matrix by significance and correlation
#'
#' Step 1: Remove traits with zero variants below Bonferroni threshold.
#' Step 2: Iterative greedy removal of traits with Pearson r > threshold.
#'
#' @param z_matrix Z-score matrix (variants x traits)
#' @param trait_files Named list of trait file paths (for reading p-values)
#' @param variant_ids Character vector of VAR_IDs (rows of z_matrix)
#' @param correlation_threshold Max Pearson r between traits (default 0.85)
#' @return List: z_matrix (filtered), removed_traits (character), correlation_matrix
filter_correlated_traits <- function(z_matrix, trait_files, variant_ids,
                                     correlation_threshold = 0.85) {
  trait_names <- colnames(z_matrix)
  n_variants <- nrow(z_matrix)
  removed <- character()

  # Step 1: Bonferroni significance check
  p_bonf <- 0.05 / n_variants
  cat(sprintf("  Bonferroni threshold: p < %.2e (0.05 / %d variants)\n", p_bonf, n_variants))

  nonsig_traits <- character()
  for (trait in trait_names) {
    if (!trait %in% names(trait_files)) next
    trait_dt <- tryCatch(
      fread(trait_files[[trait]], select = c("VAR_ID", "P_VALUE"), data.table = TRUE),
      error = function(e) NULL
    )
    if (is.null(trait_dt)) next

    # Check how many universe variants have p < bonferroni in this trait
    trait_dt <- trait_dt[VAR_ID %in% variant_ids]
    n_sig <- sum(trait_dt$P_VALUE < p_bonf, na.rm = TRUE)
    if (n_sig == 0) {
      nonsig_traits <- c(nonsig_traits, trait)
      cat(sprintf("    Removing trait %s: 0 variants with p < %.2e\n", trait, p_bonf))
    }
  }

  if (length(nonsig_traits) > 0) {
    keep_cols <- !trait_names %in% nonsig_traits
    # Guard: don't remove if it would leave < 3 traits
    if (sum(keep_cols) >= 3) {
      z_matrix <- z_matrix[, keep_cols, drop = FALSE]
      removed <- c(removed, nonsig_traits)
      trait_names <- colnames(z_matrix)
    } else {
      cat("    WARNING: Removing non-significant traits would leave < 3 traits. Skipping.\n")
    }
  }

  # Step 2: Pairwise correlation filter
  if (ncol(z_matrix) < 2) {
    return(list(z_matrix = z_matrix, removed_traits = removed,
                correlation_matrix = matrix(1, 1, 1, dimnames = list(trait_names, trait_names))))
  }

  cor_mat <- cor(z_matrix, use = "pairwise.complete.obs")
  cat(sprintf("  Checking trait correlations (threshold: |r| > %.2f)\n", correlation_threshold))

  # Iterative greedy removal
  repeat {
    # Find all above-threshold pairs (upper triangle only)
    above <- which(abs(cor_mat) > correlation_threshold & upper.tri(cor_mat), arr.ind = TRUE)
    if (nrow(above) == 0) break

    # Guard: don't remove if it would leave < 3 traits
    if (ncol(z_matrix) <= 3) {
      cat("    WARNING: Only 3 traits remain. Stopping correlation filter.\n")
      break
    }

    # Count how many above-threshold pairs each trait is involved in
    involved <- c(above[, 1], above[, 2])
    trait_counts <- sort(table(involved), decreasing = TRUE)

    # The trait with the most violations gets considered for removal
    worst_idx <- as.integer(names(trait_counts)[1])

    # For each pair involving worst_idx, the tiebreaker is non-zero count
    # Remove the trait in the pair with fewer non-zero entries
    partner_indices <- unique(c(above[above[, 1] == worst_idx, 2],
                                 above[above[, 2] == worst_idx, 1]))

    worst_nonzero <- sum(z_matrix[, worst_idx] != 0)
    remove_worst <- TRUE
    for (pi in partner_indices) {
      partner_nonzero <- sum(z_matrix[, pi] != 0)
      if (partner_nonzero < worst_nonzero) {
        remove_worst <- FALSE
        break
      }
    }

    if (remove_worst) {
      remove_idx <- worst_idx
    } else {
      # Find the partner with fewest non-zero entries
      nonzero_counts <- sapply(partner_indices, function(pi) sum(z_matrix[, pi] != 0))
      remove_idx <- partner_indices[which.min(nonzero_counts)]
    }

    remove_name <- colnames(z_matrix)[remove_idx]
    cat(sprintf("    Removing trait %s (|r| > %.2f with %d other traits)\n",
                remove_name, correlation_threshold, length(partner_indices)))
    removed <- c(removed, remove_name)
    z_matrix <- z_matrix[, -remove_idx, drop = FALSE]
    cor_mat <- cor_mat[-remove_idx, -remove_idx, drop = FALSE]
  }

  # Recompute full correlation matrix on final trait set for diagnostics
  final_cor <- cor(z_matrix, use = "pairwise.complete.obs")

  list(z_matrix = z_matrix, removed_traits = removed, correlation_matrix = final_cor)
}


#' Expand Z-score matrix to non-negative format for bNMF
#'
#' For each trait column, creates two columns: trait_pos = max(z, 0) and trait_neg = max(-z, 0).
#'
#' @param z_matrix Matrix of z-scores (variants x traits)
#' @return Non-negative matrix (variants x 2*traits)
expand_to_nonneg <- function(z_matrix) {
  trait_names <- colnames(z_matrix)
  n_variants <- nrow(z_matrix)
  n_traits <- ncol(z_matrix)

  nonneg <- matrix(0, nrow = n_variants, ncol = 2 * n_traits)
  col_names <- character(2 * n_traits)

  for (i in seq_len(n_traits)) {
    pos_idx <- 2 * i - 1
    neg_idx <- 2 * i
    nonneg[, pos_idx] <- pmax(z_matrix[, i], 0)
    nonneg[, neg_idx] <- pmax(-z_matrix[, i], 0)
    col_names[pos_idx] <- paste0(trait_names[i], "_pos")
    col_names[neg_idx] <- paste0(trait_names[i], "_neg")
  }

  colnames(nonneg) <- col_names
  rownames(nonneg) <- rownames(z_matrix)
  nonneg
}
