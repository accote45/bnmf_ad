#!/usr/bin/env Rscript
#
# A2: GTEx Tissue Expression Heatmaps for META Clusters
#
# For each cluster, maps variants to nearest protein-coding genes, then tests
# whether those genes are more highly expressed than background genes in each
# GTEx tissue (one-sided t-test, alternative = "greater").
# Produces per-cluster ComplexHeatmap figures with 2 columns:
#   Absolute Expression, Relative Expression (Specificity)
# Arranged in a 4x3 global layout via plotgardener.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(png)
  library(plotgardener)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
source(file.path(base_dir, "scripts/a1_analysis/figure_utils.R"))

# -- CLI args -----------------------------------------------------------------
# Optional single-gene-set mode: when --gene_list_file is supplied, the script
# skips the W-matrix cluster assignment / window gene-mapping (Steps 2-3) and the
# multi-panel grid (Step 7), and instead runs the SAME Abs/Rel t-tests + heatmap
# on one externally-supplied gene set (e.g. a published trait's nearest genes).
# Default (no flags) behaviour is unchanged.
args <- commandArgs(trailingOnly = TRUE)
gene_list_file <- NULL; single_label <- NULL; out_dir_arg <- NULL
ai <- 1
while (ai <= length(args)) {
  if (args[ai] == "--gene_list_file") { gene_list_file <- args[ai + 1]; ai <- ai + 2 }
  else if (args[ai] == "--label")     { single_label   <- args[ai + 1]; ai <- ai + 2 }
  else if (args[ai] == "--out_dir")   { out_dir_arg    <- args[ai + 1]; ai <- ai + 2 }
  else ai <- ai + 1
}
single_mode <- !is.null(gene_list_file)
if (single_mode && is.null(single_label)) single_label <- "GENESET"

# -- Config -------------------------------------------------------------------

out_dir  <- if (!is.null(out_dir_arg)) out_dir_arg else file.path(base_dir, "results/a2_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gtf_path       <- file.path(base_dir, "reference/gencode.v19.annotation.gtf.gz")
w_path         <- file.path(base_dir, "results/a1_analysis/META/W_matrix_META.tsv")
gtex_tar       <- "/sc/arion/projects/paul_oreilly/data/Functional_Genomics/gtex/gtex_rnaseq/median_expression/qced_data/gtex.tar.gz"
gtex_spec_path <- file.path(base_dir, "data/gtex/GTEx_specificity.csv")

# Ordered in the desired META display order (Glycemic -> Lpa). The combined-panel
# figures split this vector as [1:6] (part1) and [7:10] (part2), so the order
# here directly sets the panel layout order.
cluster_labels <- c(
  K10 = "Glycemic",
  K9  = "Obesity",
  K4  = "SHBG",
  K2  = "Adiponectin",
  K7  = "Triglycerides-HDL",
  K8  = "ALP-LDL",
  K6  = "Metabolic",
  K3  = "Platelet",
  K5  = "Blood Pressure-Stature",
  K1  = "Lpa"
)

message("=== A2: GTEx Tissue Expression Heatmaps ===\n")

# -- Step 1: Parse GTF for protein-coding genes with Ensembl IDs --------------

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

message(sprintf("  %d protein-coding genes on autosomes", nrow(gene_df)))

# -- Steps 2-3: Build cluster_genes -------------------------------------------
# Default mode: assign W-matrix variants to clusters, then map to genes in a
# +/-50kb window. Single-gene-set mode (--gene_list_file): read one external
# ENSG list as a single "cluster" and skip the W-matrix / window mapping.
if (!single_mode) {

# -- Step 2: Assign variants to max-weight cluster ----------------------------

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

message(sprintf("  %d variants assigned to %d clusters",
                nrow(variant_clusters), n_distinct(variant_clusters$cluster)))

# -- Step 3: Map variants to protein-coding genes within ±500 kb -------------

gene_window <- 50000  # bp
message(sprintf("Step 3: Mapping variants to genes within +/-%d kb...", gene_window / 1000))

var_dt  <- as.data.table(variant_clusters)
gene_dt <- as.data.table(gene_df)

window_genes <- var_dt[, {
  g <- gene_dt[chr == .BY$chr]
  if (nrow(g) == 0) {
    list(gene_name = NA_character_, gene_id = NA_character_, distance = NA_real_)
  } else {
    dist_to_gene <- pmax(g$start - pos, 0) + pmax(pos - g$end, 0)
    in_window <- dist_to_gene <= gene_window
    if (any(in_window)) {
      list(gene_name = g$gene_name[in_window],
           gene_id   = g$gene_id[in_window],
           distance  = dist_to_gene[in_window])
    } else {
      list(gene_name = NA_character_, gene_id = NA_character_, distance = NA_real_)
    }
  }
}, by = .(VAR_ID, chr, pos, cluster)]

cluster_genes <- window_genes[!is.na(gene_id),
                              .(genes = list(unique(gene_id))),
                              by = cluster]

for (i in seq_len(nrow(cluster_genes))) {
  cl <- cluster_genes$cluster[i]
  ng <- length(cluster_genes$genes[[i]])
  message(sprintf("  %s (%s): %d unique genes",
                  cl, cluster_labels[cl], ng))
}

# Gene mapping summary: n_genes, mean distance, sd distance per cluster
gene_mapping_summary <- window_genes[!is.na(gene_id),
  .(n_genes      = uniqueN(gene_id),
    mean_dist_bp = mean(distance),
    sd_dist_bp   = sd(distance)),
  by = cluster
][order(cluster)]
gene_mapping_summary[, cluster_label := cluster_labels[cluster]]
setcolorder(gene_mapping_summary, c("cluster", "cluster_label"))

gene_map_csv <- file.path(out_dir, "gene_mapping_summary.csv")
fwrite(gene_mapping_summary, gene_map_csv)
message(sprintf("  Saved %s", basename(gene_map_csv)))

} else {

# -- Steps 2-3 (single-gene-set mode): one external gene set ------------------
message(sprintf("Single-gene-set mode: reading gene list for '%s'", single_label))
gl_raw <- fread(gene_list_file, header = TRUE)
# Accept a column named ensemblID/gene_id/ENSG/Name, else fall back to column 1
ensg_col <- intersect(c("ensemblID", "gene_id", "ENSG", "Name"), names(gl_raw))
ensg <- if (length(ensg_col) > 0) gl_raw[[ensg_col[1]]] else gl_raw[[1]]
ensg <- unique(str_remove(trimws(as.character(ensg)), "\\.\\d+$"))  # strip ENSG version
ensg <- ensg[!is.na(ensg) & ensg != ""]
cluster_genes <- data.table(cluster = single_label, genes = list(ensg))
# Reassign cluster_labels to a one-element named vector so the Step 6 loop
# `for (cl in names(cluster_labels))` and label lookups run once for this set.
cluster_labels <- setNames(single_label, single_label)
message(sprintf("  %s: %d unique genes from %s",
                single_label, length(ensg), basename(gene_list_file)))
}

# -- Step 4: Load and transform GTEx expression -------------------------------

message("Step 4: Loading GTEx expression data...")

gtex_tmp <- tempdir()
untar(gtex_tar, exdir = gtex_tmp)
gtex_expr <- fread(file.path(gtex_tmp, "gtex/GTEx_expression.csv"))

tissues <- setdiff(names(gtex_expr), "Name")
message(sprintf("  %d genes x %d tissues", nrow(gtex_expr), length(tissues)))

# log2(TPM + 1) transform
gtex_expr[, (tissues) := lapply(.SD, function(x) log2(x + 1)), .SDcols = tissues]

# Keep only protein-coding genes
gtex_expr <- gtex_expr[Name %in% gene_df$gene_id]
message(sprintf("  %d protein-coding genes with GTEx data", nrow(gtex_expr)))

# Load relative expression (Skene procedure specificity)
message("  Loading GTEx specificity data...")
gtex_spec <- fread(gtex_spec_path)
gtex_spec <- gtex_spec[Name %in% gene_df$gene_id]
message(sprintf("  %d protein-coding genes with specificity values", nrow(gtex_spec)))

# In single-gene-set mode, report which input genes are absent from the GTEx
# protein-coding matrices (e.g. non-coding entries) so they aren't mistaken for
# a bug during visual comparison. run_tissue_ttests silently uses only matches.
if (single_mode) {
  ts_genes <- cluster_genes$genes[[1]]
  missing_dt <- data.table(
    gene_id        = ts_genes,
    in_gtex_expr   = ts_genes %in% gtex_expr$Name,
    in_gtex_spec   = ts_genes %in% gtex_spec$Name
  )[in_gtex_expr == FALSE | in_gtex_spec == FALSE]
  miss_csv <- file.path(out_dir, sprintf("%s_genes_missing_from_gtex.csv",
                                         tolower(single_label)))
  fwrite(missing_dt, miss_csv)
  message(sprintf("  %d/%d input genes matched GTEx expr; %d missing (see %s)",
                  sum(ts_genes %in% gtex_expr$Name), length(ts_genes),
                  nrow(missing_dt), basename(miss_csv)))
}

# -- Step 5: T-tests for all gene set x expression type combinations ----------

message("Step 5: Running one-sided t-tests...")

run_tissue_ttests <- function(expr_dt, test_genes, tissues, filter_expressed = FALSE) {
  results <- vector("list", length(tissues))
  for (ti in seq_along(tissues)) {
    tis <- tissues[ti]
    vals <- expr_dt[[tis]]
    ids  <- expr_dt$Name

    if (filter_expressed) {
      keep <- vals > 0
      vals <- vals[keep]
      ids  <- ids[keep]
    }

    is_test  <- ids %in% test_genes
    test_val <- vals[is_test]
    ctrl_val <- vals[!is_test]
    n_t <- length(test_val)
    n_c <- length(ctrl_val)

    if (n_t >= 2 && n_c >= 2) {
      tt <- t.test(test_val, ctrl_val, alternative = "greater")
      results[[ti]] <- data.table(
        tissue = tis, t_statistic = as.numeric(tt$statistic),
        p_value = tt$p.value, n_test = n_t, n_control = n_c,
        mean_test = mean(test_val), mean_control = mean(ctrl_val)
      )
    } else {
      results[[ti]] <- data.table(
        tissue = tis, t_statistic = NA_real_,
        p_value = NA_real_, n_test = n_t, n_control = n_c,
        mean_test = NA_real_, mean_control = NA_real_
      )
    }
  }
  rbindlist(results)
}

all_results <- list()

for (i in seq_len(nrow(cluster_genes))) {
  cl       <- cluster_genes$cluster[i]
  cl_genes <- cluster_genes$genes[[i]]

  # Absolute expression: nearest genes
  res <- run_tissue_ttests(gtex_expr, cl_genes, tissues, filter_expressed = TRUE)
  res[, `:=`(cluster = cl, gene_set = "nearest", expr_type = "abs")]
  all_results[[length(all_results) + 1]] <- res

  # Relative expression: nearest genes
  res_r <- run_tissue_ttests(gtex_spec, cl_genes, tissues, filter_expressed = FALSE)
  res_r[, `:=`(cluster = cl, gene_set = "nearest", expr_type = "rel")]
  all_results[[length(all_results) + 1]] <- res_r

  message(sprintf("  %s done", cl))
}

results_dt <- rbindlist(all_results, use.names = TRUE)
results_dt[, log10_p := -log10(pmax(p_value, 1e-300))]

out_csv <- file.path(out_dir, "gtex_ttest_results_nearest.csv")
fwrite(results_dt, out_csv)
message(sprintf("  Saved %s (%d rows)", basename(out_csv), nrow(results_dt)))

# -- Step 6: ComplexHeatmap generation -----------------------------------------

message("Step 6: Generating ComplexHeatmap figures...")

# Format tissue names for display
format_tissue <- function(x) str_replace_all(x, "([a-z])([A-Z])", "\\1 \\2")

# Alphabetical tissue order (A at top)
tissue_labels <- sort(unique(format_tissue(tissues)))

# Shared color scale: -log10(p), range set to observed max (ceiling)
max_log10p <- ceiling(max(results_dt$log10_p, na.rm = TRUE))
max_log10p <- max(max_log10p, 1)  # ensure minimum range of 1
message(sprintf("  Color scale range: 0 to %d", max_log10p))

col_fun <- colorRamp2(
  seq(0, max_log10p, length.out = 10),
  colorRampPalette(c("#F0FAF0", "#B2E2B2", "#66C266", "#2E8B57", "#1B4D2E"))(10)
)

#' Build a single-column Heatmap for one gene_set x expr_type combination
build_hm_column <- function(results_dt, cl, gene_set_val, expr_type_val,
                            col_title, show_row, show_legend,
                            tissue_labels, col_fun) {
  sub <- results_dt[cluster == cl & gene_set == gene_set_val & expr_type == expr_type_val]
  sub[, tissue_label := format_tissue(tissue)]

  mat <- matrix(NA_real_, nrow = length(tissue_labels), ncol = 1,
                dimnames = list(tissue_labels, ""))
  for (r in seq_len(nrow(sub))) {
    tl <- sub$tissue_label[r]
    if (tl %in% tissue_labels) mat[tl, 1] <- sub$log10_p[r]
  }

  mat[mat > max_log10p] <- max_log10p
  mat[is.na(mat)] <- 0

  Heatmap(
    mat,
    col = col_fun,
    width = unit(10, "mm"),
    show_row_dend = FALSE,
    row_order = tissue_labels,
    column_title = col_title,
    column_title_rot = 0,
    column_title_gp = gpar(fontsize = 11),
    rect_gp = gpar(col = "white", lwd = 2),
    show_heatmap_legend = show_legend,
    column_labels = "",
    show_row_names = show_row,
    row_names_side = "left",
    row_names_max_width = unit(10, "cm"),
    row_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "-logP", direction = "horizontal"),
    name = paste(cl, gene_set_val, expr_type_val, sep = "_")
  )
}

# -- Generate per-cluster PNGs ------------------------------------------------

cluster_png_paths <- list()

for (cl in names(cluster_labels)) {
  hm1 <- build_hm_column(results_dt, cl, "nearest", "abs", "Abs.",
                          show_row = TRUE, show_legend = FALSE,
                          tissue_labels, col_fun)
  hm2 <- build_hm_column(results_dt, cl, "nearest", "rel", "Rel.",
                          show_row = FALSE, show_legend = TRUE,
                          tissue_labels, col_fun)

  hm_list <- hm1 + hm2

  png_path <- file.path(out_dir, sprintf("gtex_heatmap_CH_%s.png", cl))
  png(png_path, width = 6, height = 12, units = "in", res = 600)
  draw(hm_list,
       ht_gap = unit(4, "mm"),
       show_heatmap_legend = FALSE)

  dev.off()
  cluster_png_paths[[cl]] <- png_path
  message(sprintf("  Saved gtex_heatmap_CH_%s.png", cl))
}

# -- Step 7: Assemble two grid layouts with plotgardener ----------------------
# Skipped in single-gene-set mode (the K1-K6 / K7-K10 grid is cluster-specific).
if (!single_mode) {

message("Step 7: Assembling split grid layouts (K1-K6 and K7-K10)...")

cluster_imgs <- lapply(names(cluster_labels), function(cl) {
  readPNG(cluster_png_paths[[cl]])
})
names(cluster_imgs) <- names(cluster_labels)

panel_w <- 6
panel_h <- 12
label_h <- 0.6
slot_h  <- label_h + panel_h
gap_x   <- 0.25
gap_y   <- 0.15
top_margin <- 0.25

abs_x0 <- 3.80
abs_x1 <- 4.10
rel_x0 <- 4.20
rel_x1 <- 4.50
label_dy <- -0.15

add_exp_labels <- function(x_off, raster_y) {
  seg_y <- raster_y - 0.05
  txt_y <- seg_y + label_dy

  plotSegments(
    x0 = x_off + abs_x0, x1 = x_off + abs_x1,
    y0 = seg_y, y1 = seg_y,
    default.units = "inches", lwd = 1.25, linecolor = "#37a7db"
  )
  plotText(
    label = "Abs. Exp", fontsize = 7, fontcolor = "#37a7db",
    x = x_off + (abs_x0 + abs_x1) / 2, y = txt_y,
    just = "center", default.units = "inches"
  )

  plotSegments(
    x0 = x_off + rel_x0, x1 = x_off + rel_x1,
    y0 = seg_y, y1 = seg_y,
    default.units = "inches", lwd = 1.25, linecolor = "#37a7db"
  )
  plotText(
    label = "Rel. Exp", fontsize = 7, fontcolor = "#37a7db",
    x = x_off + (rel_x0 + rel_x1) / 2, y = txt_y,
    just = "center", default.units = "inches"
  )
}

place_panels <- function(clusters, grid_ncol, grid_nrow, out_path, letter_offset = 0) {
  page_w <- panel_w * grid_ncol + gap_x * (grid_ncol - 1)
  page_h <- slot_h * grid_nrow + gap_y * (grid_nrow - 1) + top_margin

  png(out_path, width = page_w, height = page_h, units = "in", res = 600)
  pageCreate(width = page_w, height = page_h, default.units = "inches",
             showGuides = FALSE)

  for (i in seq_along(clusters)) {
    cl   <- clusters[i]
    col  <- (i - 1) %% grid_ncol
    row  <- (i - 1) %/% grid_ncol
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
      fontcolor = "darkblue",
      x = x_pos + panel_w / 2, y = slot_y + label_h / 2,
      just = "center", default.units = "inches"
    )
  }

  dev.off()
  message(sprintf("  Saved %s", basename(out_path)))
}

# Part 1: K1-K6 in 2x3 grid
place_panels(
  clusters = names(cluster_labels)[1:6],
  grid_ncol = 3, grid_nrow = 2,
  out_path = file.path(out_dir, "gtex_heatmap_combined_CH_part1.png"),
  letter_offset = 0
)

# Part 2: K7-K10 in 2x2 grid
place_panels(
  clusters = names(cluster_labels)[7:10],
  grid_ncol = 2, grid_nrow = 2,
  out_path = file.path(out_dir, "gtex_heatmap_combined_CH_part2.png"),
  letter_offset = 6
)

}  # end if (!single_mode): Step 7 grid assembly

message("\nDone.")
