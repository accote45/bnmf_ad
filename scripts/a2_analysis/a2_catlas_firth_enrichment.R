#!/usr/bin/env Rscript
#
# A2 (SUPPLEMENTARY): CATLAS cell-type enrichment via matched-null Firth regression
# ============================================================================
# This is an ADDITIVE supplementary analysis. It does NOT replace or modify the
# existing permutation-based scripts (a2_catlas_{pancreas,esophagus}_heatmaps.R,
# a2_catlas_wss_summary.R), which remain the primary CATLAS analysis.
#
# It ports the published SNP-level methodology:
#   - Unit of analysis is the SNP (not the nearest gene).
#   - Per cluster, "lead SNPs" (Y=1) are the bNMF variants assigned to that
#     cluster (max-weight hard assignment).
#   - "Null SNPs" (Y=0) are 1000G variants within +/-50 kb of a cluster lead
#     that are NOT in LD (r2 < 0.05) with ANY lead SNP in the cluster.
#   - Each SNP is annotated for direct overlap with cell-type ATAC peaks (X_ij)
#     and with genic elements (CDS / 5'UTR / 3'UTR) from GENCODE v19.
#   - Enrichment per cell type is tested with Firth bias-reduced logistic
#     regression: logit P(Y=1) = a0 + a_exon*CDS + a_5utr*5UTR + a_3utr*3UTR + theta_i*X_i
#     The theta_i p-value (logistf penalized LR test) is the enrichment test.
#   - Bonferroni correction over the cell types tested, per cluster.
#
# Differences from the published method that are intentional here (per scope):
#   - CATLAS panel limited to the cell types already downloaded (pancreas +
#     esophagus), NOT the full ~328-cell-type atlas.
#   - Genic annotation from GENCODE v19 (hg19) instead of Ensembl 104 (GRCh38),
#     because our coordinates are hg19. "Protein-coding exon" => CDS.
#   - LD / null matching from the 1000G EUR panel (single panel).
#
# Usage:
#   Rscript scripts/a2_analysis/a2_catlas_firth_enrichment.R
# (Run under the rocky9 R 4.2.0 recipe; see project memory.)

suppressPackageStartupMessages({
  library(data.table)
  library(logistf)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"

# -- Config -------------------------------------------------------------------

out_dir <- file.path(base_dir, "results/a2_analysis/catlas/catlas_firth")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gtf_path   <- file.path(base_dir, "reference/gencode.v19.annotation.gtf.gz")
w_path     <- file.path(base_dir, "results/a1_analysis/META/W_matrix_META.tsv")
panel_dir  <- file.path(base_dir, "reference/1kg_eur")
panel_pfx  <- "1000G.EUR.QC"             # per-chromosome plink bfile prefix
plink_bin  <- "/hpc/packages/minerva-centos7/plink/1.90b6.21/plink"

work_dir   <- file.path(out_dir, "ld_work")
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

# --- Analysis parameters (reported to stdout at start) -----------------------
null_window  <- 50000   # candidate null must be within +/-50 kb of a cluster lead
ld_window_kb <- 250     # plink LD scan window (>= null_window, to catch LD with nearby leads)
r2_thresh    <- 0.05    # null SNP must have r2 < this to ALL cluster leads
min_leads    <- 10      # exclude clusters with fewer matched lead SNPs (empirical)
bonf_alpha   <- 0.05    # family-wise alpha for Bonferroni (per cluster, over cell types)
plink_mem    <- 4000    # MB
plink_thr    <- 2

cluster_labels <- c(
  K1  = "Lpa",
  K2  = "Adiponectin",
  K3  = "Platelet",
  K4  = "SHBG",
  K5  = "Blood Pressure-Stature",
  K6  = "Metabolic",
  K7  = "Triglycerides-HDL",
  K8  = "ALP-LDL",
  K9  = "Obesity",
  K10 = "Glycemic"
)

# CATLAS cell-type panel: reuse the already-downloaded hg19 BED files only.
tissue_configs <- list(
  Pancreas = list(
    ccre_dir = file.path(base_dir, "data/catlas/pancreas/hg19"),
    labels = c(
      Pancreatic_Acinar_Cell         = "Acinar Cell",
      Pancreatic_Alpha_Cell_1        = "Alpha Cell 1",
      Pancreatic_Alpha_Cell_2        = "Alpha Cell 2",
      Pancreatic_Beta_Cell_1         = "Beta Cell 1",
      Pancreatic_Beta_Cell_2         = "Beta Cell 2",
      `Pancreatic_Delta,Gamma_cell`  = "Delta/Gamma Cell",
      Ductal_Cell_Pancreatic         = "Ductal Cell",
      Fetal_Pancreatic_Acinar_Cell_1 = "Fetal Acinar Cell 1",
      Fetal_Pancreatic_Acinar_Cell_2 = "Fetal Acinar Cell 2",
      Fetal_Pancreatic_Ductal_Cell   = "Fetal Ductal Cell",
      Fetal_Pancreatic_Islet_Cell    = "Fetal Islet Cell"
    )
  ),
  Esophagus = list(
    ccre_dir = file.path(base_dir, "data/catlas/esophagus/hg19"),
    labels = c(
      Esophageal_Epithelial_Cell            = "Epithelial Cell",
      Smooth_Muscle_Esophageal_Mucosal      = "Smooth Muscle (Mucosal)",
      Smooth_Muscle_Esophageal_Muscularis_1 = "Smooth Muscle (Muscularis 1)",
      Smooth_Muscle_Esophageal_Muscularis_2 = "Smooth Muscle (Muscularis 2)",
      Smooth_Muscle_Esophageal_Muscularis_3 = "Smooth Muscle (Muscularis 3)"
    )
  )
)

message("=== A2 SUPPLEMENTARY: CATLAS matched-null Firth enrichment ===\n")
message("Parameters:")
message(sprintf("  null window         : +/-%d bp", null_window))
message(sprintf("  LD scan window      : %d kb", ld_window_kb))
message(sprintf("  r2 threshold (null) : < %.2f", r2_thresh))
message(sprintf("  min lead SNPs/cluster: %d", min_leads))
message(sprintf("  LD/null panel       : %s", panel_dir))
message(sprintf("  Bonferroni alpha    : %.2f (per cluster, over cell types)", bonf_alpha))
message(sprintf("  genic annotation    : GENCODE v19 (CDS / 5'UTR / 3'UTR)"))
message("")

# -- Step 1: Lead SNPs (max-weight hard assignment) ---------------------------

message("Step 1: Assigning lead SNPs to clusters (max-weight)...")

w_mat <- fread(w_path)
k_cols <- grep("^K\\d+$", names(w_mat), value = TRUE)

leads <- w_mat[, .(VAR_ID)]
leads[, cluster := k_cols[max.col(as.matrix(w_mat[, k_cols, with = FALSE]))]]
leads[, chr := sub("_.*$", "", VAR_ID)]
leads[, pos := as.integer(sub("^[^_]+_([0-9]+)_.*$", "\\1", VAR_ID))]
leads <- leads[chr %in% as.character(1:22)]

message(sprintf("  %d lead SNPs across %d clusters", nrow(leads), uniqueN(leads$cluster)))

# -- Step 2: Match leads to 1000G bim by chr:pos (recover plink rsIDs) ---------

message("Step 2: Matching leads to 1000G panel by chr:pos...")

leads[, bim_rsid := NA_character_]
for (cc in sort(unique(leads$chr))) {
  bim_path <- file.path(panel_dir, sprintf("%s.%s.bim", panel_pfx, cc))
  if (!file.exists(bim_path)) next
  bim <- fread(bim_path, header = FALSE,
               col.names = c("chr", "rsid", "cm", "bp", "a1", "a2"))
  bim <- bim[!duplicated(bp)]              # first variant per position
  setkey(bim, bp)
  idx <- leads$chr == cc
  leads[idx, bim_rsid := bim[.(leads$pos[idx]), rsid, on = "bp"]]
}
n_matched <- leads[!is.na(bim_rsid), .N]
message(sprintf("  %d/%d leads matched to a panel variant by position",
                n_matched, nrow(leads)))

all_lead_pos  <- unique(leads[, paste(chr, pos, sep = ":")])
all_lead_rsid <- unique(leads[!is.na(bim_rsid), bim_rsid])

# -- Step 3: plink LD scan per chromosome -------------------------------------

message("Step 3: Running plink --r2 per chromosome...")

ld_list <- list()
for (cc in sort(unique(leads$chr))) {
  lead_rsids <- leads[chr == cc & !is.na(bim_rsid), unique(bim_rsid)]
  if (length(lead_rsids) == 0) next

  out_pfx  <- file.path(work_dir, sprintf("ld_chr%s", cc))
  ld_path  <- paste0(out_pfx, ".ld")

  # Reuse a cached .ld if present (plink scan is the slow step).
  if (!file.exists(ld_path)) {
    snp_file <- file.path(work_dir, sprintf("leads_chr%s.txt", cc))
    writeLines(lead_rsids, snp_file)

    args <- c(
      "--bfile", file.path(panel_dir, sprintf("%s.%s", panel_pfx, cc)),
      "--r2",
      "--ld-snp-list", snp_file,
      "--ld-window-kb", ld_window_kb,
      "--ld-window", 99999,
      "--ld-window-r2", 0,
      "--memory", plink_mem,
      "--threads", plink_thr,
      "--out", out_pfx
    )
    # Clear LD_LIBRARY_PATH for the child so the centos7 plink uses system libs
    # (the R process runs under the rocky9 gcc-14 LD_LIBRARY_PATH).
    system2(plink_bin, args, env = "LD_LIBRARY_PATH=",
            stdout = FALSE, stderr = FALSE)
  }

  if (file.exists(ld_path)) {
    ld_list[[cc]] <- fread(ld_path)
  } else {
    message(sprintf("  chr%s: no .ld output", cc))
  }
}

ld_dt <- rbindlist(ld_list, use.names = TRUE)
message(sprintf("  %d LD pairs across %d chromosomes",
                nrow(ld_dt), length(ld_list)))

# Map lead (SNP_A) -> cluster
lead_rsid_cluster <- leads[!is.na(bim_rsid), .(bim_rsid, cluster)]
lead_rsid_cluster <- lead_rsid_cluster[!duplicated(bim_rsid)]
ld_dt[lead_rsid_cluster, cluster := i.cluster, on = c(SNP_A = "bim_rsid")]

# -- Step 4: Build null SNP sets per cluster ----------------------------------

message("Step 4: Constructing matched null SNP sets per cluster...")

null_list <- list()
cluster_n <- list()
for (cl in names(cluster_labels)) {
  lead_cl <- leads[cluster == cl]
  n_lead  <- nrow(lead_cl)

  pairs <- ld_dt[cluster == cl]
  if (nrow(pairs) == 0) {
    cluster_n[[cl]] <- list(n_lead = n_lead, n_null = 0L)
    next
  }
  pairs[, dist := abs(BP_B - BP_A)]

  # Per candidate (SNP_B): max r2 to any cluster lead, min distance to a lead
  cand <- pairs[, .(maxr2 = max(R2), mindist = min(dist)),
                by = .(chr = CHR_B, pos = BP_B, rsid = SNP_B)]

  # Null: within 50 kb of a lead, r2 < threshold to ALL leads, not itself a lead
  cand[, pos_key := paste(chr, pos, sep = ":")]
  nulls <- cand[mindist <= null_window &
                  maxr2 < r2_thresh &
                  !(rsid %in% all_lead_rsid) &
                  !(pos_key %in% all_lead_pos)]
  nulls <- nulls[!duplicated(pos_key)]

  cluster_n[[cl]] <- list(n_lead = n_lead, n_null = nrow(nulls))

  if (nrow(nulls) > 0) {
    null_list[[cl]] <- nulls[, .(chr = as.character(chr), pos, cluster = cl)]
  }
  message(sprintf("  %s (%s): %d leads, %d matched null SNPs",
                  cl, cluster_labels[cl], n_lead, nrow(nulls)))
}

# Determine which clusters are testable (enough leads AND a null set)
testable <- names(cluster_labels)[
  vapply(names(cluster_labels), function(cl) {
    cn <- cluster_n[[cl]]
    cn$n_lead >= min_leads && cn$n_null >= min_leads
  }, logical(1))
]
# Display columns in the desired META visualization order (Glycemic -> Lpa);
# any testable clusters not in the list fall back to natural K order at the end.
desired_k_order <- c("K10", "K9", "K4", "K2", "K7", "K8", "K6", "K3", "K5", "K1")
testable <- c(desired_k_order[desired_k_order %in% testable],
              sort(setdiff(testable, desired_k_order)))
excluded <- setdiff(names(cluster_labels), testable)
message(sprintf("\n  Testable clusters: %s", paste(testable, collapse = ", ")))
if (length(excluded) > 0) {
  message(sprintf("  Excluded (n_lead<%d or n_null<%d): %s",
                  min_leads, min_leads,
                  paste(sprintf("%s(%s)", excluded, cluster_labels[excluded]),
                        collapse = ", ")))
}

# -- Step 5: Master SNP table (unique positions) ------------------------------

message("\nStep 5: Building master SNP annotation table...")

master <- unique(rbindlist(list(
  leads[, .(chr, pos)],
  if (length(null_list)) rbindlist(null_list)[, .(chr, pos)] else NULL
)))
master <- master[chr %in% as.character(1:22)]
master[, pos_key := paste(chr, pos, sep = ":")]
master[, `:=`(start = pos, end = pos)]
setkey(master, chr, start, end)
message(sprintf("  %d unique SNP positions to annotate", nrow(master)))

# -- Step 6: Genic annotation (GENCODE v19: CDS / 5'UTR / 3'UTR) --------------

message("Step 6: Annotating genic context (GENCODE v19)...")

gtf <- fread(cmd = sprintf("zcat %s | grep -v '^#'", gtf_path),
             sep = "\t", header = FALSE,
             col.names = c("chr", "source", "feature", "start", "end",
                           "score", "strand", "frame", "attr"))
gtf <- gtf[feature %in% c("CDS", "UTR") &
             grepl('gene_type "protein_coding"', attr, fixed = TRUE)]
gtf[, chr := sub("^chr", "", chr)]
gtf <- gtf[chr %in% as.character(1:22)]
gtf[, transcript_id := sub('.*transcript_id "([^"]+)".*', "\\1", attr)]

# CDS bounds per transcript (to orient UTRs as 5' vs 3')
cds <- gtf[feature == "CDS"]
cds_bounds <- cds[, .(cds_min = min(start), cds_max = max(end)), by = transcript_id]

utr <- gtf[feature == "UTR"]
utr[cds_bounds, `:=`(cds_min = i.cds_min, cds_max = i.cds_max),
    on = "transcript_id"]
# Orient: on + strand, UTR upstream of CDS (end <= cds_min) is 5'; downstream is 3'.
# On - strand the assignment is reversed.
utr[, side := fcase(
  is.na(cds_min), "unknown",
  strand == "+" & end   <= cds_min, "utr5",
  strand == "+" & start >= cds_max, "utr3",
  strand == "-" & start >= cds_max, "utr5",
  strand == "-" & end   <= cds_min, "utr3",
  default = "unknown"
)]

annotate_overlap <- function(master_dt, iv_dt) {
  # Returns 0/1 vector over master rows for overlap with interval set iv_dt.
  if (nrow(iv_dt) == 0) return(rep(0L, nrow(master_dt)))
  iv <- iv_dt[, .(chr, start, end)]
  setkey(iv, chr, start, end)
  ov <- foverlaps(master_dt[, .(chr, start, end)], iv,
                  by.x = c("chr", "start", "end"), nomatch = NULL, which = TRUE)
  hit <- rep(0L, nrow(master_dt))
  if (nrow(ov) > 0) hit[unique(ov$xid)] <- 1L
  hit
}

master[, exon := annotate_overlap(master, cds[, .(chr, start, end)])]
master[, utr5 := annotate_overlap(master, utr[side == "utr5", .(chr, start, end)])]
master[, utr3 := annotate_overlap(master, utr[side == "utr3", .(chr, start, end)])]
message(sprintf("  SNPs in CDS: %d | 5'UTR: %d | 3'UTR: %d",
                sum(master$exon), sum(master$utr5), sum(master$utr3)))

# -- Step 7: Peak annotation (X_ij per cell type) -----------------------------

message("Step 7: Annotating ATAC peak overlap per cell type...")

cell_types  <- character(0)
ct_tissue   <- character(0)
ct_label    <- character(0)
for (tissue_name in names(tissue_configs)) {
  tc <- tissue_configs[[tissue_name]]
  bed_files <- list.files(tc$ccre_dir, pattern = "\\.bed$", full.names = TRUE)
  bed_files <- bed_files[!grepl("_unmapped\\.bed$", bed_files)]
  for (bf in bed_files) {
    ct <- tools::file_path_sans_ext(basename(bf))
    bed <- fread(bf, header = FALSE, select = 1:3,
                 col.names = c("chr", "start", "end"))
    bed[, chr := sub("^chr", "", chr)]
    bed <- bed[chr %in% as.character(1:22)]
    bed[, start := start + 1L]   # BED 0-based half-open -> 1-based inclusive
    col <- paste0("X__", ct)
    master[, (col) := annotate_overlap(master, bed)]
    cell_types <- c(cell_types, ct)
    ct_tissue  <- c(ct_tissue, tissue_name)
    nice <- tc$labels[ct]
    ct_label   <- c(ct_label, if (is.na(nice)) gsub("_", " ", ct) else unname(nice))
    message(sprintf("  %-40s %d / %d SNPs overlap a peak",
                    ct, sum(master[[col]]), nrow(master)))
  }
}
ct_meta <- data.table(cell_type = cell_types, tissue = ct_tissue, label = ct_label)

# -- Step 8: Firth regression per cluster x cell type -------------------------

message("\nStep 8: Fitting Firth logistic regressions...")

fit_one <- function(dat, xcol) {
  # Build formula using only genic covariates that vary (avoid collinearity).
  covars <- c("exon", "utr5", "utr3")
  covars <- covars[vapply(covars, function(c) length(unique(dat[[c]])) > 1, logical(1))]
  rhs <- paste(c(covars, "Xv"), collapse = " + ")
  d <- data.frame(Y = dat$Y, Xv = dat[[xcol]], dat[, ..covars])
  fit <- tryCatch(
    logistf(as.formula(paste("Y ~", rhs)), data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(beta = NA_real_, p = NA_real_))
  list(beta = unname(coef(fit)["Xv"]), p = unname(fit$prob["Xv"]))
}

results <- list()
for (cl in testable) {
  lead_cl <- leads[cluster == cl, .(pos_key = paste(chr, pos, sep = ":"), Y = 1L)]
  null_cl <- rbindlist(null_list[cl])[, .(pos_key = paste(chr, pos, sep = ":"), Y = 0L)]
  dat_cl  <- unique(rbindlist(list(lead_cl, null_cl)), by = "pos_key")
  dat_cl[master, `:=`(exon = i.exon, utr5 = i.utr5, utr3 = i.utr3), on = "pos_key"]

  # attach all peak columns
  xcols <- paste0("X__", cell_types)
  dat_cl[master, (xcols) := mget(paste0("i.", xcols)), on = "pos_key"]

  for (i in seq_along(cell_types)) {
    xcol <- xcols[i]
    if (length(unique(dat_cl[[xcol]])) < 2) next   # no variation -> untestable
    fr <- fit_one(dat_cl, xcol)
    results[[length(results) + 1]] <- data.table(
      cluster   = cl,
      cell_type = cell_types[i],
      tissue    = ct_meta$tissue[i],
      label     = ct_meta$label[i],
      n_lead    = sum(dat_cl$Y == 1),
      n_null    = sum(dat_cl$Y == 0),
      n_lead_peak = sum(dat_cl$Y == 1 & dat_cl[[xcol]] == 1),
      beta      = fr$beta,
      p_value   = fr$p
    )
  }
  message(sprintf("  %s (%s): fitted %d cell types",
                  cl, cluster_labels[cl],
                  sum(vapply(results, function(r) r$cluster == cl, logical(1)))))
}

res_dt <- rbindlist(results)

# -- Step 9: Bonferroni correction (per cluster) + save -----------------------

message("\nStep 9: Bonferroni correction and saving results...")

res_dt[, p_bonferroni := pmin(p_value * .N, 1), by = cluster]
res_dt[, log10_p := -log10(pmax(p_value, 1e-300))]
res_dt[, sig := p_bonferroni < bonf_alpha]
res_dt[, cluster_label := cluster_labels[cluster]]
setorder(res_dt, cluster, p_value)

for (cl in testable) {
  ns <- res_dt[cluster == cl & sig == TRUE, .N]
  message(sprintf("  %s (%s): %d/%d cell types significant (Bonferroni < %.2f)",
                  cl, cluster_labels[cl],
                  ns, res_dt[cluster == cl, .N], bonf_alpha))
}

out_csv <- file.path(out_dir, "catlas_firth_enrichment_results.csv")
fwrite(res_dt, out_csv)
message(sprintf("  Saved %s (%d rows)", basename(out_csv), nrow(res_dt)))

# -- Step 10: Heatmap (cell types x clusters, -log10 p) -----------------------

message("\nStep 10: Generating heatmap...")

if (nrow(res_dt) > 0) {
  row_labels <- sprintf("%s — %s", res_dt$tissue, res_dt$label)
  res_dt[, row_label := sprintf("%s — %s", tissue, label)]

  all_rows <- unique(res_dt[, .(row_label, tissue,
                                ord = match(cell_type, ct_meta$cell_type))])
  setorder(all_rows, tissue, ord)
  row_order_lbl <- all_rows$row_label
  col_order_lbl <- cluster_labels[testable]

  mat <- matrix(NA_real_, nrow = length(row_order_lbl), ncol = length(testable),
                dimnames = list(row_order_lbl, col_order_lbl))
  bmat <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(res_dt))) {
    rl <- res_dt$row_label[i]
    cl <- cluster_labels[res_dt$cluster[i]]
    mat[rl, cl]  <- res_dt$log10_p[i]
    bmat[rl, cl] <- res_dt$sig[i]
  }

  max_val <- max(c(mat[is.finite(mat)], 1), na.rm = TRUE)
  col_fun <- colorRamp2(
    seq(0, ceiling(max_val), length.out = 10),
    colorRampPalette(c("#F2F0F7", "#CBC9E2", "#9E9AC8", "#756BB1", "#54278F"))(10)
  )

  tissue_colors <- c("Pancreas" = "#66C2A5", "Esophagus" = "#FC8D62")
  row_tissue <- all_rows$tissue
  row_ha <- rowAnnotation(
    Tissue = row_tissue,
    col = list(Tissue = tissue_colors),
    annotation_legend_param = list(Tissue = list(title = "Tissue"))
  )

  hm <- Heatmap(
    mat,
    col = col_fun,
    na_col = "white",
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = row_order_lbl, column_order = col_order_lbl,
    row_names_side = "left",
    row_names_max_width = unit(12, "cm"),
    row_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 11),
    rect_gp = gpar(col = "white", lwd = 1),
    left_annotation = row_ha,
    name = "firth",
    heatmap_legend_param = list(title = "-log10(p)", direction = "horizontal"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (isTRUE(bmat[i, j]))
        grid.text("*", x, y, gp = gpar(fontsize = 16, col = "black"))
    }
  )

  n_rows  <- nrow(mat)
  panel_h <- max(5, n_rows * 0.3 + 2)
  panel_w <- max(7, length(testable) * 0.9 + 4)

  hm_path <- file.path(out_dir, "catlas_firth_enrichment_heatmap.png")
  png(hm_path, width = panel_w, height = panel_h, units = "in", res = 600)
  draw(hm,
       column_title = sprintf(
         "CATLAS Firth enrichment (matched null, Bonferroni < %.2f marked *)",
         bonf_alpha),
       column_title_gp = gpar(fontsize = 13, fontface = "bold"),
       heatmap_legend_side = "bottom")
  dev.off()
  message(sprintf("  Saved %s", basename(hm_path)))
}

# -- Step 11: Per-tissue, per-cluster Firth heatmaps --------------------------
# Mirrors the per-cluster figure style of the primary CATLAS scripts, but with
# a single enrichment column (Firth theta -log10 p) instead of Abs./Rel.
# Significant cell types (Bonferroni < alpha) are marked with *.

message("\nStep 11: Per-tissue per-cluster Firth heatmaps...")

suppressPackageStartupMessages({
  library(plotgardener)
  library(png)
})

# Shared color scale across all per-cluster panels (global max -log10 p)
glob_max <- max(c(res_dt$log10_p[is.finite(res_dt$log10_p)], 1), na.rm = TRUE)
col_fun_pc <- colorRamp2(
  seq(0, ceiling(glob_max), length.out = 10),
  colorRampPalette(c("#F2F0F7", "#CBC9E2", "#9E9AC8", "#756BB1", "#54278F"))(10)
)

build_firth_col <- function(sub, ct_order_labels, cl_title, col_fun, show_legend) {
  mat  <- matrix(0, nrow = length(ct_order_labels), ncol = 1,
                 dimnames = list(ct_order_labels, ""))
  bvec <- setNames(rep(FALSE, length(ct_order_labels)), ct_order_labels)
  for (i in seq_len(nrow(sub))) {
    mat[sub$label[i], 1] <- sub$log10_p[i]
    bvec[sub$label[i]]   <- isTRUE(sub$sig[i])
  }
  Heatmap(
    mat, col = col_fun, width = unit(10, "mm"),
    cluster_rows = FALSE, row_order = ct_order_labels,
    column_title = cl_title, column_title_rot = 0,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 2),
    show_heatmap_legend = show_legend, column_labels = "",
    row_names_side = "left", row_names_max_width = unit(10, "cm"),
    row_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "-log10(p)", direction = "horizontal"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (isTRUE(bvec[i]))
        grid.text("*", x, y, gp = gpar(fontsize = 14, col = "black"))
    }
  )
}

# Grid placement of per-cluster PNGs (mirrors the primary scripts' layout)
place_grid <- function(imgs, clusters, grid_ncol, panel_w, panel_h,
                       out_path, labels) {
  label_h <- 0.5; gap_x <- 0.25; gap_y <- 0.15; top_margin <- 0.2
  grid_nrow <- ceiling(length(clusters) / grid_ncol)
  slot_h <- label_h + panel_h
  page_w <- panel_w * grid_ncol + gap_x * (grid_ncol - 1)
  page_h <- slot_h * grid_nrow + gap_y * (grid_nrow - 1) + top_margin

  png(out_path, width = page_w, height = page_h, units = "in", res = 600)
  pageCreate(width = page_w, height = page_h, default.units = "inches",
             showGuides = FALSE)
  for (i in seq_along(clusters)) {
    cl  <- clusters[i]
    col <- (i - 1) %% grid_ncol
    row <- (i - 1) %/% grid_ncol
    x_pos  <- col * (panel_w + gap_x)
    slot_y <- top_margin + row * (slot_h + gap_y)
    plotRaster(image = imgs[[cl]], x = x_pos, y = slot_y + label_h,
               width = panel_w, height = panel_h,
               just = c("left", "top"), default.units = "inches")
    plotText(label = letters[i], fontsize = 16, fontface = "bold",
             x = x_pos + 0.05, y = slot_y + 0.12, just = "left",
             default.units = "inches")
    plotText(label = labels[cl], fontsize = 14, fontface = "bold",
             fontcolor = "#54278F",
             x = x_pos + panel_w / 2, y = slot_y + label_h / 2,
             just = "center", default.units = "inches")
  }
  dev.off()
  message(sprintf("    Saved %s", basename(out_path)))
}

for (tissue_name in names(tissue_configs)) {
  ct_t <- ct_meta[tissue == tissue_name][order(match(cell_type, ct_meta$cell_type))]
  ct_order_labels <- ct_t$label

  tissue_out <- file.path(base_dir, "results/a2_analysis/catlas",
                          paste0("catlas_", tolower(tissue_name)), "firth")
  dir.create(tissue_out, showWarnings = FALSE, recursive = TRUE)

  panel_h <- max(2.2, length(ct_order_labels) * 0.32 + 1.2)
  message(sprintf("  %s -> %s", tissue_name, tissue_out))

  for (cl in testable) {
    sub <- res_dt[cluster == cl & tissue == tissue_name]
    hm  <- build_firth_col(sub, ct_order_labels, cluster_labels[cl],
                           col_fun_pc, show_legend = TRUE)
    png_path <- file.path(tissue_out,
                          sprintf("catlas_%s_firth_%s.png", tolower(tissue_name), cl))
    png(png_path, width = 5, height = panel_h, units = "in", res = 600)
    draw(hm, heatmap_legend_side = "bottom")
    dev.off()
  }
  message(sprintf("    Saved %d per-cluster panels", length(testable)))

  # Combined grid across testable clusters
  imgs <- list()
  for (cl in testable) {
    p <- file.path(tissue_out,
                   sprintf("catlas_%s_firth_%s.png", tolower(tissue_name), cl))
    if (file.exists(p)) imgs[[cl]] <- readPNG(p)
  }
  if (length(imgs) > 0) {
    ref <- imgs[[names(imgs)[1]]]
    pw <- 5
    ph <- pw * (dim(ref)[1] / dim(ref)[2])
    place_grid(imgs, names(imgs), grid_ncol = 3, panel_w = pw, panel_h = ph,
               out_path = file.path(tissue_out,
                          sprintf("catlas_%s_firth_combined.png", tolower(tissue_name))),
               labels = cluster_labels)
  }
}

message("\nDone.")
