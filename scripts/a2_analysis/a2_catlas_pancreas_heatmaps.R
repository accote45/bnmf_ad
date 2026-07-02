#!/usr/bin/env Rscript
#
# A2: CATLAS Pancreatic cCRE Enrichment Heatmaps for META Clusters
#
# Tests whether genes near each bNMF cluster's variants are enriched for
# chromatin accessibility (cCREs from sci-ATAC-seq) in pancreatic cell types.
# Uses gene-level cCRE counts (±50kb window) as proxy for regulatory activity.
#
# Produces per-cluster ComplexHeatmap figures with 2 columns:
#   Abs. = cCRE count enrichment
#   Rel. = cell-type specificity enrichment
#
# Usage:
#   Rscript scripts/a2_analysis/a2_catlas_pancreas_heatmaps.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
source(file.path(base_dir, "scripts/a1_analysis/figure_utils.R"))

# -- Config -------------------------------------------------------------------

out_dir <- file.path(base_dir, "results/a2_analysis/catlas/catlas_pancreas")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ccre_dir    <- file.path(base_dir, "data/catlas/pancreas/hg19")
gtf_path    <- file.path(base_dir, "reference/gencode.v19.annotation.gtf.gz")
w_path      <- file.path(base_dir, "results/a1_analysis/META/W_matrix_META.tsv")

n_perm   <- 10000
window   <- 50000  # ±50kb around gene body for cCRE mapping
set.seed(42)       # reproducible permutation p-values

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

# Display labels for pancreatic cell types (BED filename → nice label)
celltype_labels <- c(
  Pancreatic_Acinar_Cell        = "Acinar Cell",
  Pancreatic_Alpha_Cell_1       = "Alpha Cell 1",
  Pancreatic_Alpha_Cell_2       = "Alpha Cell 2",
  Pancreatic_Beta_Cell_1        = "Beta Cell 1",
  Pancreatic_Beta_Cell_2        = "Beta Cell 2",
  `Pancreatic_Delta,Gamma_cell` = "Delta/Gamma Cell",
  Ductal_Cell_Pancreatic        = "Ductal Cell",
  Fetal_Pancreatic_Acinar_Cell_1 = "Fetal Acinar Cell 1",
  Fetal_Pancreatic_Acinar_Cell_2 = "Fetal Acinar Cell 2",
  Fetal_Pancreatic_Ductal_Cell  = "Fetal Ductal Cell",
  Fetal_Pancreatic_Islet_Cell   = "Fetal Islet Cell"
)

message("=== A2: CATLAS Pancreatic cCRE Enrichment Heatmaps ===\n")

# -- Step 1: Parse GTF for protein-coding genes --------------------------------

message("Step 1: Parsing GTF for protein-coding genes...")

gene_df <- read_tsv(
  gtf_path, comment = "#",
  col_names = c("chr", "source", "feature", "start", "end",
                "score", "strand", "frame", "attributes"),
  col_types = cols(chr = col_character(), start = col_double(),
                   end = col_double(), .default = col_character())
) %>%
  filter(feature == "gene",
         str_detect(attributes, 'gene_type "protein_coding"')) %>%
  mutate(
    gene_name = str_replace(attributes, '.*gene_name "([^"]+)".*', "\\1"),
    gene_id   = str_replace(attributes, '.*gene_id "([^"]+)".*', "\\1"),
    gene_id   = str_remove(gene_id, "\\.\\d+$")
  ) %>%
  mutate(chr = str_remove(chr, "^chr")) %>%
  filter(chr %in% as.character(1:22)) %>%
  select(chr, start, end, gene_name, gene_id)

gene_dt <- as.data.table(gene_df)
message(sprintf("  %d protein-coding genes on autosomes", nrow(gene_dt)))

# -- Step 2: Load W matrix and assign variants to clusters --------------------

message("Step 2: Loading W matrix and assigning variants to clusters...")

w_mat <- read_tsv(w_path, show_col_types = FALSE)
k_cols <- grep("^K\\d+$", names(w_mat), value = TRUE)

variant_clusters <- w_mat %>%
  mutate(
    cluster = k_cols[max.col(across(all_of(k_cols)))],
    chr     = str_extract(VAR_ID, "^[^_]+"),
    pos     = as.numeric(str_extract(VAR_ID, "(?<=_)\\d+"))
  ) %>%
  select(VAR_ID, chr, pos, cluster)

var_dt <- as.data.table(variant_clusters)
message(sprintf("  %d variants assigned to %d clusters",
                nrow(var_dt), n_distinct(var_dt$cluster)))

# -- Step 3: Map variants to nearest protein-coding gene ----------------------

message("Step 3: Mapping variants to nearest protein-coding genes...")

nearest_genes <- var_dt[, {
  g <- gene_dt[chr == .BY$chr]
  if (nrow(g) == 0) {
    list(gene_name = NA_character_, gene_id = NA_character_, distance = NA_real_)
  } else {
    dist_to_gene <- pmax(g$start - pos, 0) + pmax(pos - g$end, 0)
    idx <- which.min(dist_to_gene)
    list(gene_name = g$gene_name[idx],
         gene_id   = g$gene_id[idx],
         distance  = dist_to_gene[idx])
  }
}, by = .(VAR_ID, chr, pos, cluster)]

cluster_genes_symbols <- nearest_genes[!is.na(gene_name),
                                       .(genes = list(unique(gene_name))),
                                       by = cluster]

for (i in seq_len(nrow(cluster_genes_symbols))) {
  cl <- cluster_genes_symbols$cluster[i]
  ng <- length(cluster_genes_symbols$genes[[i]])
  message(sprintf("  %s (%s): %d unique genes",
                  cl, cluster_labels[cl], ng))
}

# -- Step 4: Load cCRE BED files and build gene × cell-type matrix ------------

message("Step 4: Loading cCRE BED files and building gene × cell-type matrix...")

bed_files <- list.files(ccre_dir, pattern = "\\.bed$", full.names = TRUE)
bed_files <- bed_files[!grepl("_unmapped\\.bed$", bed_files)]

if (length(bed_files) == 0) {
  stop("No BED files found in ", ccre_dir,
       ". Run download_catlas_pancreas.sh first.")
}

cell_type_ids <- tools::file_path_sans_ext(basename(bed_files))
message(sprintf("  Found %d cell type BED files", length(cell_type_ids)))

# Gene windows: TSS ± window (use gene body boundaries extended by window)
gene_windows <- copy(gene_dt)
gene_windows[, win_start := pmax(start - window, 0L)]
gene_windows[, win_end   := end + window]
setkey(gene_windows, chr, win_start, win_end)

# Count cCREs per gene per cell type
ccre_count_mat <- matrix(0L, nrow = nrow(gene_dt), ncol = length(cell_type_ids),
                         dimnames = list(gene_dt$gene_name, cell_type_ids))

for (i in seq_along(bed_files)) {
  ct <- cell_type_ids[i]
  bed <- fread(bed_files[i], header = FALSE,
               select = 1:3,
               col.names = c("chr", "start", "end"))
  bed[, chr := sub("^chr", "", chr)]
  bed <- bed[chr %in% as.character(1:22)]
  setkey(bed, chr, start, end)

  # foverlaps: find cCREs overlapping each gene window
  overlaps <- foverlaps(bed, gene_windows, by.x = c("chr", "start", "end"),
                        by.y = c("chr", "win_start", "win_end"),
                        nomatch = NULL)

  if (nrow(overlaps) > 0) {
    gene_counts <- overlaps[, .N, by = gene_name]
    idx <- match(gene_counts$gene_name, gene_dt$gene_name)
    ccre_count_mat[idx[!is.na(idx)], i] <- gene_counts$N[!is.na(idx)]
  }

  message(sprintf("  %s: %d cCREs loaded, %d overlapping gene windows",
                  ct, nrow(bed), nrow(overlaps)))
}

# Save intermediate matrix
ccre_csv <- file.path(out_dir, "catlas_pancreas_ccre_gene_matrix.csv")
ccre_out <- data.table(gene_name = gene_dt$gene_name, as.data.table(ccre_count_mat))
fwrite(ccre_out, ccre_csv)
message(sprintf("  Saved gene × cell-type matrix: %s", basename(ccre_csv)))

# -- Step 5: Build specificity matrix -----------------------------------------

message("Step 5: Computing cell-type specificity matrix...")

row_totals <- rowSums(ccre_count_mat)
spec_mat <- ccre_count_mat / row_totals
spec_mat[!is.finite(spec_mat)] <- 0

message(sprintf("  %d genes with at least 1 cCRE across cell types",
                sum(row_totals > 0)))

# Wrap as data.tables for permutation test function
ccre_expr_dt <- data.table(Gene_name = gene_dt$gene_name, as.data.table(ccre_count_mat))
ccre_spec_dt <- data.table(Gene_name = gene_dt$gene_name, as.data.table(spec_mat))

# -- Step 6: Permutation tests ------------------------------------------------

message("Step 6: Running permutation tests (10,000 permutations)...\n")

run_celltype_permtests <- function(expr_dt, test_genes, cell_types, n_perm,
                                   filter_expressed = FALSE) {
  is_test <- expr_dt$Gene_name %in% test_genes
  n1 <- sum(is_test)
  n0 <- sum(!is_test)

  if (n1 < 2 || n0 < 2) {
    return(data.table(
      cell_type = cell_types, p_value = NA_real_,
      original_diff = NA_real_, n_test = n1, n_control = n0
    ))
  }

  expr_mat <- as.matrix(expr_dt[, ..cell_types])
  n_genes <- nrow(expr_mat)
  col_totals <- colSums(expr_mat)
  disease_sums <- colSums(expr_mat[is_test, , drop = FALSE])

  if (filter_expressed) {
    # Absolute test: restrict each cell type's background to genes with a cCRE
    # present (value > 0), matching the GTEx pipeline (filter_expressed = TRUE).
    # Genes with no cCRE contribute 0 to the sum and are dropped from the
    # per-column count, so each mean is taken over only the genes with regulatory
    # activity in that cell type. This removes the trivial inflation from
    # comparing real genes against a background padded with zero-cCRE genes.
    expr_pos    <- expr_mat > 0
    pos_totals  <- colSums(expr_pos)
    disease_pos <- colSums(expr_pos[is_test, , drop = FALSE])
    t_n <- pmax(disease_pos, 1)
    c_n <- pmax(pos_totals - disease_pos, 1)
    original_diffs <- disease_sums / t_n - (col_totals - disease_sums) / c_n

    exceed_count <- rep(0L, length(cell_types))
    for (p in seq_len(n_perm)) {
      idx <- sample.int(n_genes, n1)
      perm_sums <- colSums(expr_mat[idx, , drop = FALSE])
      perm_pos  <- colSums(expr_pos[idx, , drop = FALSE])
      pt_n <- pmax(perm_pos, 1)
      pc_n <- pmax(pos_totals - perm_pos, 1)
      perm_diffs <- perm_sums / pt_n - (col_totals - perm_sums) / pc_n
      exceed_count <- exceed_count + (perm_diffs >= original_diffs)
    }
  } else {
    # Relative (specificity) test: all genes form the background (unchanged).
    original_diffs <- disease_sums / n1 - (col_totals - disease_sums) / n0
    exceed_count <- rep(0L, length(cell_types))
    for (p in seq_len(n_perm)) {
      idx <- sample.int(n_genes, n1)
      perm_sums <- colSums(expr_mat[idx, , drop = FALSE])
      perm_diffs <- perm_sums / n1 - (col_totals - perm_sums) / n0
      exceed_count <- exceed_count + (perm_diffs >= original_diffs)
    }
  }

  p_values <- exceed_count / n_perm

  data.table(
    cell_type = cell_types,
    p_value = p_values,
    original_diff = as.numeric(original_diffs),
    n_test = n1, n_control = n0
  )
}

# Color scale: purple/violet gradient (distinct from GTEx green and TS orange)
col_fun <- colorRamp2(
  seq(0, 10, length.out = 10),
  colorRampPalette(c("#F2F0F7", "#CBC9E2", "#9E9AC8", "#756BB1", "#54278F"))(10)
)

build_hm_column <- function(pval_vec, cell_type_labels, col_title,
                            show_row, show_legend, col_fun,
                            nice_labels) {
  mat <- matrix(NA_real_, nrow = length(cell_type_labels), ncol = 1,
                dimnames = list(cell_type_labels, ""))
  for (ct in names(pval_vec)) {
    ct_label <- nice_labels[ct]
    if (!is.na(ct_label) && ct_label %in% cell_type_labels) {
      mat[ct_label, 1] <- pval_vec[ct]
    }
  }
  mat[mat > 10] <- 10
  mat[is.na(mat)] <- 0

  Heatmap(
    mat,
    col = col_fun,
    width = unit(10, "mm"),
    show_row_dend = FALSE,
    row_order = cell_type_labels,
    column_title = col_title,
    column_title_rot = 0,
    column_title_gp = gpar(fontsize = 10),
    rect_gp = gpar(col = "white", lwd = 2),
    show_heatmap_legend = show_legend,
    column_labels = "",
    row_names_side = "left",
    row_names_max_width = unit(10, "cm"),
    row_names_gp = gpar(fontsize = if (show_row) 10 else 0,
                         col = "black"),
    show_row_names = show_row,
    heatmap_legend_param = list(
      title = "-logP", direction = "horizontal",
      at = seq(-log10(1), -log10(0.0001), length = 5)
    )
  )
}

all_results <- list()
cell_type_display <- sort(unname(celltype_labels))

for (cl in names(cluster_labels)) {
  cl_genes <- cluster_genes_symbols[cluster == cl]$genes[[1]]
  cl_genes_in_data <- cl_genes[cl_genes %in% ccre_expr_dt$Gene_name]

  message(sprintf("  %s (%s): %d/%d cluster genes in cCRE matrix",
                  cl, cluster_labels[cl],
                  length(cl_genes_in_data), length(cl_genes)))

  if (length(cl_genes_in_data) < 2) {
    message(sprintf("    Skipping %s: too few genes", cl))
    next
  }

  # Absolute: cCRE count enrichment (background = genes with a cCRE present)
  res_abs <- run_celltype_permtests(ccre_expr_dt, cl_genes_in_data,
                                    cell_type_ids, n_perm, filter_expressed = TRUE)
  res_abs[, `:=`(cluster = cl, expr_type = "abs")]

  # Relative: specificity enrichment
  res_rel <- run_celltype_permtests(ccre_spec_dt, cl_genes_in_data,
                                    cell_type_ids, n_perm)
  res_rel[, `:=`(cluster = cl, expr_type = "rel")]

  all_results <- c(all_results, list(res_abs, res_rel))

  make_pval_vec <- function(res) {
    setNames(-log10(pmax(res$p_value, 1e-300)), res$cell_type)
  }

  pv_abs <- make_pval_vec(res_abs)
  pv_rel <- make_pval_vec(res_rel)

  hm1 <- build_hm_column(pv_abs, cell_type_display, "Abs.",
                          show_row = TRUE, show_legend = FALSE,
                          col_fun, celltype_labels)
  hm2 <- build_hm_column(pv_rel, cell_type_display, "Rel.",
                          show_row = FALSE, show_legend = TRUE,
                          col_fun, celltype_labels)

  hm_list <- hm1 + hm2

  n_ct <- length(cell_type_display)
  panel_h <- max(4, n_ct * 0.3 + 1.5)

  png_path <- file.path(out_dir, sprintf("catlas_pancreas_%s.png", cl))
  png(png_path, width = 6, height = panel_h, units = "in", res = 600)
  draw(hm_list,
       ht_gap = unit(4, "mm"),
       column_title = sprintf("Pancreas — %s", cluster_labels[cl]),
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       show_heatmap_legend = FALSE)
  dev.off()
  message(sprintf("    Saved %s", basename(png_path)))
}

# -- Step 7: Save results CSV -------------------------------------------------

message("\nStep 7: Saving results CSV...")
results_dt <- rbindlist(all_results, use.names = TRUE)
results_dt[, log10_p := -log10(pmax(p_value, 1e-300))]

out_csv <- file.path(out_dir, "catlas_pancreas_permtest_results.csv")
fwrite(results_dt, out_csv)
message(sprintf("  Saved %s (%d rows)", basename(out_csv), nrow(results_dt)))

# -- Step 8: Combined grid figure ---------------------------------------------

message("\nStep 8: Assembling combined grid figure...")

suppressPackageStartupMessages({
  library(plotgardener)
  library(png)
})

panel_w  <- 6
label_h  <- 0.6
gap_x    <- 0.25
gap_y    <- 0.15
top_margin <- 0.25

place_panels <- function(cluster_imgs, clusters, grid_ncol, grid_nrow,
                         panel_w, panel_h, out_path, cluster_labels,
                         letter_offset = 0) {
  slot_h <- label_h + panel_h
  page_w <- panel_w * grid_ncol + gap_x * (grid_ncol - 1)
  page_h <- slot_h * grid_nrow + gap_y * (grid_nrow - 1) + top_margin

  png(out_path, width = page_w, height = page_h, units = "in", res = 600)
  pageCreate(width = page_w, height = page_h, default.units = "inches",
             showGuides = FALSE)

  for (i in seq_along(clusters)) {
    cl <- clusters[i]
    col <- (i - 1) %% grid_ncol
    row <- (i - 1) %/% grid_ncol
    x_pos    <- col * (panel_w + gap_x)
    slot_y   <- top_margin + row * (slot_h + gap_y)
    raster_y <- slot_y + label_h

    plotRaster(
      image = cluster_imgs[[cl]],
      x = x_pos, y = raster_y,
      width = panel_w, height = panel_h,
      just = c("left", "top"),
      default.units = "inches"
    )

    plotText(
      label = letters[i + letter_offset], fontsize = 18, fontface = "bold",
      fontfamily = "arial",
      x = x_pos + 0.05, y = slot_y + 0.15,
      just = "left", default.units = "inches"
    )

    plotText(
      label = cluster_labels[cl], fontsize = 16, fontface = "bold",
      fontcolor = "#54278F",
      x = x_pos + panel_w / 2, y = slot_y + label_h / 2,
      just = "center", default.units = "inches"
    )
  }

  dev.off()
  message(sprintf("  Saved %s", basename(out_path)))
}

# Load all individual PNGs
imgs <- list()
for (cl in names(cluster_labels)) {
  png_path <- file.path(out_dir, sprintf("catlas_pancreas_%s.png", cl))
  if (!file.exists(png_path)) {
    message(sprintf("  Warning: %s not found, skipping", basename(png_path)))
    next
  }
  imgs[[cl]] <- readPNG(png_path)
}

if (length(imgs) > 0) {
  ref_img <- imgs[[names(imgs)[1]]]
  panel_h_px <- dim(ref_img)[1]
  panel_w_px <- dim(ref_img)[2]
  panel_h <- panel_h_px / (panel_w_px / panel_w)

  # Part 1: K1–K6 (3×2)
  part1_cls <- intersect(names(cluster_labels)[1:6], names(imgs))
  if (length(part1_cls) > 0) {
    place_panels(
      cluster_imgs = imgs,
      clusters = part1_cls,
      grid_ncol = 3, grid_nrow = 2,
      panel_w = panel_w, panel_h = panel_h,
      out_path = file.path(out_dir, "catlas_pancreas_combined_part1.png"),
      cluster_labels = cluster_labels,
      letter_offset = 0
    )
  }

  # Part 2: K7–K10 (2×2)
  part2_cls <- intersect(names(cluster_labels)[7:10], names(imgs))
  if (length(part2_cls) > 0) {
    place_panels(
      cluster_imgs = imgs,
      clusters = part2_cls,
      grid_ncol = 2, grid_nrow = 2,
      panel_w = panel_w, panel_h = panel_h,
      out_path = file.path(out_dir, "catlas_pancreas_combined_part2.png"),
      cluster_labels = cluster_labels,
      letter_offset = 6
    )
  }
}

message("\nDone.")
