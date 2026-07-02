#!/usr/bin/env Rscript
#
# A2: Tabula Sapiens Cell-Type Expression Heatmaps for META Clusters
#
# For each cluster, maps variants to nearest protein-coding genes, then tests
# whether those genes are more highly expressed than background genes in each
# cell type per tissue using a permutation approach (10,000 permutations).
# Produces per-tissue × per-cluster ComplexHeatmap figures with 2 columns:
#   Absolute Expression, Relative Expression (Specificity)
# With a cell count (NCells) barplot annotation on the left.

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

out_dir <- file.path(base_dir, "results/a2_analysis/tabula_sapiens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ts_data_dir <- file.path(base_dir,
  "data/tabula_sapiens/tabula_sapiens_cell_types/tabula_sapiens_processed_data")
gtf_path    <- file.path(base_dir, "reference/gencode.v19.annotation.gtf.gz")
w_path      <- file.path(base_dir, "results/a1_analysis/META/W_matrix_META.tsv")
all_meta_path <- file.path(ts_data_dir, "all_metadata.csv")

n_perm <- 10000  # number of permutations
set.seed(42)     # reproducible permutation p-values

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

# Discover available tissues from expression CSV files
tissue_files <- list.files(ts_data_dir, pattern = "_expression\\.csv$", full.names = FALSE)
tissues <- str_remove(tissue_files, "_expression\\.csv$")
message(sprintf("Found %d tissues: %s", length(tissues), paste(tissues, collapse = ", ")))

stopifnot(
  "all_metadata.csv not found" = file.exists(all_meta_path)
)

# Load all_metadata for nice cell type display labels
all_meta_dt <- fread(all_meta_path)

# -- Duplicate cell type merge rules -----------------------------------------
# Each entry: list of expression column names to merge → kept column name
# Expression/specificity values are averaged; NCells are summed.
merge_rules <- list(
  Blood = list(
    list(cols = c("NK_cell", "nk_cell"), keep = "NK_cell", label = "NK Cells")
  ),
  Large_Intestine = list(
    list(cols = c("cd8pos__T_cell", "cd8pos_t_cell"),
         keep = "cd8pos__T_cell", label = "CD8-positive alpha-beta T Cells"),
    list(cols = c("intestinal_transient_amplifying_cell", "transit_amplifying_cell_of_lg_int"),
         keep = "intestinal_transient_amplifying_cell", label = "Intestinal Transient Amplifying Cells")
  ),
  Small_Intestine = list(
    list(cols = c("Cd8pos__T_cell", "cd8pos_t_cell"),
         keep = "Cd8pos__T_cell", label = "CD8 alpha-beta T Cells"),
    list(cols = c("intestinal_transient_amplifying_cell", "transit_amplifying_cell_of_small_int"),
         keep = "intestinal_transient_amplifying_cell", label = "Intestinal Transient Amplifying Cells"),
    list(cols = c("paneth_cell", "paneth_cell_of_epithelium_of_small_int"),
         keep = "paneth_cell", label = "Paneth Cells")
  ),
  Lymph_Node = list(
    list(cols = c("CD141positive_myeloid_dendritic_cell", "cd141positive_myeloid_dendritic_cell"),
         keep = "CD141positive_myeloid_dendritic_cell", label = "CD141 Myeloid Dendritic Cells"),
    list(cols = c("CD1cpositive_myeloid_dendritic_cell", "cd1cpositive_myeloid_dendritic_cell"),
         keep = "CD1cpositive_myeloid_dendritic_cell", label = "CD1c Myeloid Dendritic Cells"),
    list(cols = c("regulatory_T_cell", "regulatory_t_cell"),
         keep = "regulatory_T_cell", label = "Regulatory T Cells")
  ),
  Thymus = list(
    list(cols = c("Macrophages", "macrophage"),
         keep = "Macrophages", label = "Macrophages"),
    list(cols = c("Monocytes", "monocyte"),
         keep = "Monocytes", label = "Monocytes"),
    list(cols = c("Natural_killer_cells", "nk_cell"),
         keep = "Natural_killer_cells", label = "NK Cells"),
    list(cols = c("Plasma_cells", "plasma_cell"),
         keep = "Plasma_cells", label = "Plasma Cells")
  )
)

#' Apply merge rules to expression/specificity data.tables and metadata
#' Averages expression values for merged columns, sums NCells in metadata.
apply_merges <- function(expr_dt, spec_dt, meta_dt, cell_types, rules) {
  drop_cols <- character(0)
  for (rule in rules) {
    present <- rule$cols[rule$cols %in% cell_types]
    if (length(present) < 2) next
    keep <- rule$keep
    others <- setdiff(present, keep)

    # Average expression across duplicate columns
    for (dt in list(expr_dt, spec_dt)) {
      avg <- rowMeans(as.matrix(dt[, ..present]), na.rm = TRUE)
      dt[, (keep) := avg]
    }

    # Sum NCells for the merged rows in metadata (positional)
    keep_idx <- match(keep, cell_types)
    other_idx <- match(others, cell_types)
    valid_idx <- c(keep_idx, other_idx)
    valid_idx <- valid_idx[!is.na(valid_idx) & valid_idx <= nrow(meta_dt)]
    if (length(valid_idx) > 0 && !is.na(keep_idx) && keep_idx <= nrow(meta_dt)) {
      meta_dt[keep_idx, NCells := sum(meta_dt$NCells[valid_idx])]
    }

    drop_cols <- c(drop_cols, others)
    message(sprintf("    Merged: %s → %s", paste(present, collapse = " + "), rule$label))
  }

  # Drop merged-away columns
  if (length(drop_cols) > 0) {
    drop_present <- drop_cols[drop_cols %in% names(expr_dt)]
    if (length(drop_present) > 0) {
      expr_dt[, (drop_present) := NULL]
      spec_dt[, (drop_present) := NULL]
    }
    # Drop corresponding metadata rows
    drop_idx <- match(drop_cols, cell_types)
    drop_idx <- drop_idx[!is.na(drop_idx) & drop_idx <= nrow(meta_dt)]
    if (length(drop_idx) > 0) {
      meta_dt <- meta_dt[-drop_idx]
    }
    cell_types <- setdiff(cell_types, drop_cols)
  }

  list(expr_dt = expr_dt, spec_dt = spec_dt, meta_dt = meta_dt, cell_types = cell_types)
}

# Label overrides for cases where all_metadata has incorrect labels
# Key = expression column name, value = corrected display label
label_overrides <- c(
  memory_b_cell = "Memory B Cells"  # Thymus: all_metadata mislabels as "B Cells"
)

message("=== A2: Tabula Sapiens Cell-Type Expression Heatmaps ===\n")

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

message(sprintf("  %d protein-coding genes on autosomes", nrow(gene_df)))

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

# -- Step 3: Map variants to nearest protein-coding gene ----------------------

message("Step 3: Mapping variants to nearest protein-coding genes...")

var_dt  <- as.data.table(variant_clusters)
gene_dt <- as.data.table(gene_df)

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

# Store both gene symbols and Ensembl IDs per cluster
cluster_genes_symbols <- nearest_genes[!is.na(gene_name),
                                       .(genes = list(unique(gene_name))),
                                       by = cluster]

for (i in seq_len(nrow(cluster_genes_symbols))) {
  cl <- cluster_genes_symbols$cluster[i]
  ng <- length(cluster_genes_symbols$genes[[i]])
  message(sprintf("  %s (%s): %d unique genes",
                  cl, cluster_labels[cl], ng))
}

# -- Step 4: Permutation test function (vectorized) ---------------------------

message("Step 4: Setting up permutation test...")

#' Vectorized permutation tests across all cell types simultaneously
#'
#' For each permutation, shuffles the disease/control labels and computes
#' the difference in means for ALL cell types at once using matrix operations.
#' This is orders of magnitude faster than looping per cell type.
#'
#' @param expr_dt data.table with Gene_name + cell type columns
#' @param test_genes Character vector of gene symbols to test
#' @param cell_types Character vector of cell type column names
#' @param n_perm Number of permutations
#' @return data.table with cell_type, p_value, original_diff, n_test, n_control
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

  # Convert to matrix for fast column operations
  expr_mat <- as.matrix(expr_dt[, ..cell_types])
  n_genes <- nrow(expr_mat)
  col_totals <- colSums(expr_mat)
  disease_sums <- colSums(expr_mat[is_test, , drop = FALSE])

  if (filter_expressed) {
    # Absolute test: restrict each cell type's background to EXPRESSED genes
    # (value > 0), matching the GTEx pipeline (filter_expressed = TRUE). Silent
    # genes contribute 0 to the sum and are dropped from the per-column count, so
    # every mean is taken over only the genes expressed in that cell type. This
    # removes the trivial inflation from comparing real (expressed) genes against
    # a background padded with non-expressed genes.
    expr_pos    <- expr_mat > 0
    pos_totals  <- colSums(expr_pos)                            # expressed genes / cell type
    disease_pos <- colSums(expr_pos[is_test, , drop = FALSE])   # expressed test genes / cell type
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

# -- Step 5: Process each tissue ----------------------------------------------

message("Step 5: Running permutation tests across all tissues...\n")

all_results <- list()

# Shared color scale: -log10(p) from 0 to 10
col_fun <- colorRamp2(
  seq(0, 10, length.out = 10),
  colorRampPalette(c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"))(10)
)

# Format cell type names for display (underscores to spaces)
format_celltype <- function(x) str_replace_all(x, "_", " ")

#' Build a single-column Heatmap for one gene_set x expr_type combination
#' @param pval_vec Named numeric vector (-log10 p), names are expression column names
#' @param cell_type_labels Character vector of nice display labels (alphabetical)
#' @param nice_labels Named character vector: expr col name → nice label
build_hm_column <- function(pval_vec, cell_type_labels, col_title,
                            show_row, show_legend, col_fun,
                            nice_labels, row_ha = NULL) {
  # Build matrix: rows = nice labels, 1 column
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

  hm_args <- list(
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

  # Add left annotation only to first column
  if (!is.null(row_ha)) {
    hm_args$left_annotation <- row_ha
  }

  do.call(Heatmap, hm_args)
}

for (tissue in tissues) {
  message(sprintf("--- Processing tissue: %s ---", tissue))

  # Load expression, specificity, and metadata
  expr_path <- file.path(ts_data_dir, sprintf("%s_expression.csv", tissue))
  spec_path <- file.path(ts_data_dir, sprintf("%s_specificity.csv", tissue))
  meta_path <- file.path(ts_data_dir, sprintf("%s_metadata.csv", tissue))

  expr_dt <- fread(expr_path)
  spec_dt <- fread(spec_path)
  meta_dt <- fread(meta_path)

  cell_types_raw <- setdiff(names(expr_dt), "Gene_name")

  # Filter out placeholder columns (e.g., "notexisting" in Prostate)
  placeholder <- grepl("^not.?exist", cell_types_raw, ignore.case = TRUE)
  if (any(placeholder)) {
    message(sprintf("  Dropping placeholder column(s): %s",
                    paste(cell_types_raw[placeholder], collapse = ", ")))
    expr_dt[, (cell_types_raw[placeholder]) := NULL]
    spec_dt[, (cell_types_raw[placeholder]) := NULL]
  }
  cell_types <- cell_types_raw[!placeholder]

  # Filter out placeholder/extra rows from per-tissue metadata to align
  # with expression columns (e.g., Prostate has "Mast Cell", "not-existing")
  meta_dt <- meta_dt[!grepl("not.?exist", Cell_Type, ignore.case = TRUE)]
  if (nrow(meta_dt) > length(cell_types)) {
    meta_dt <- meta_dt[seq_along(cell_types)]
  }
  message(sprintf("  %d genes x %d cell types", nrow(expr_dt), length(cell_types)))

  # Build nice display labels from all_metadata.csv BEFORE merging
  # (positional match: all_metadata rows align with expression columns)
  tissue_meta <- all_meta_dt[Tissue == tissue]
  if (nrow(tissue_meta) == length(cell_types)) {
    nice_labels <- setNames(tissue_meta$Cell_Type, cell_types)
  } else {
    message(sprintf("  Warning: all_metadata count (%d) != expression cols (%d), using formatted names",
                    nrow(tissue_meta), length(cell_types)))
    nice_labels <- setNames(format_celltype(cell_types), cell_types)
  }

  # Apply label overrides for mislabeled entries in all_metadata
  for (col_name in names(label_overrides)) {
    if (col_name %in% names(nice_labels)) {
      nice_labels[col_name] <- label_overrides[col_name]
    }
  }

  # Apply merge rules for duplicate cell types (before transforms)
  if (tissue %in% names(merge_rules)) {
    message(sprintf("  Applying merge rules for %s...", tissue))
    merged <- apply_merges(expr_dt, spec_dt, meta_dt, cell_types, merge_rules[[tissue]])
    expr_dt    <- merged$expr_dt
    spec_dt    <- merged$spec_dt
    meta_dt    <- merged$meta_dt
    cell_types <- merged$cell_types

    # Update nice_labels: use merge rule labels for kept columns, drop merged-away
    for (rule in merge_rules[[tissue]]) {
      if (rule$keep %in% cell_types) {
        nice_labels[rule$keep] <- rule$label
      }
    }
    nice_labels <- nice_labels[cell_types]
    message(sprintf("  After merging: %d cell types", length(cell_types)))
  }

  # log2(CPM + 1) transform for absolute expression
  expr_dt[, (cell_types) := lapply(.SD, function(x) log2(x + 1)), .SDcols = cell_types]

  # Alphabetically sorted display labels
  cell_type_labels <- sort(unname(nice_labels))

  # Run permutation tests for each cluster
  for (cl in names(cluster_labels)) {
    cl_genes <- cluster_genes_symbols[cluster == cl]$genes[[1]]
    cl_genes_in_tissue <- cl_genes[cl_genes %in% expr_dt$Gene_name]

    message(sprintf("  %s (%s): %d/%d cluster genes in %s",
                    cl, cluster_labels[cl],
                    length(cl_genes_in_tissue), length(cl_genes), tissue))

    # Skip if fewer than 2 cluster genes present
    if (length(cl_genes_in_tissue) < 2) {
      message(sprintf("    Skipping %s_%s: too few genes", tissue, cl))
      next
    }

    # Permutation tests: nearest genes x abs/rel
    # Nearest - Absolute
    res_na <- run_celltype_permtests(expr_dt, cl_genes_in_tissue, cell_types, n_perm,
                                     filter_expressed = TRUE)
    res_na[, `:=`(tissue = tissue, cluster = cl, gene_set = "nearest", expr_type = "abs")]

    # Nearest - Relative (specificity)
    res_nr <- run_celltype_permtests(spec_dt, cl_genes_in_tissue, cell_types, n_perm)
    res_nr[, `:=`(tissue = tissue, cluster = cl, gene_set = "nearest", expr_type = "rel")]

    all_results <- c(all_results, list(res_na, res_nr))

    # Compute -log10(p) vectors for heatmap
    make_pval_vec <- function(res) {
      pv <- setNames(-log10(pmax(res$p_value, 1e-300)), res$cell_type)
      pv
    }

    pv_na <- make_pval_vec(res_na)
    pv_nr <- make_pval_vec(res_nr)

    # Build 2-column heatmap
    hm1 <- build_hm_column(pv_na, cell_type_labels, "Abs.",
                            show_row = TRUE, show_legend = FALSE,
                            col_fun, nice_labels)
    hm2 <- build_hm_column(pv_nr, cell_type_labels, "Rel.",
                            show_row = FALSE, show_legend = TRUE,
                            col_fun, nice_labels)

    hm_list <- hm1 + hm2

    # Determine panel height based on number of cell types
    n_ct <- length(cell_type_labels)
    panel_h <- max(4, n_ct * 0.3 + 1.5)

    png_path <- file.path(out_dir,
                          sprintf("tabula_sapiens_%s_%s.png", tissue, cl))
    png(png_path, width = 6, height = panel_h, units = "in", res = 600)
    draw(hm_list,
         ht_gap = unit(4, "mm"),
         column_title = sprintf("%s — %s", tissue, cluster_labels[cl]),
         column_title_gp = gpar(fontsize = 14, fontface = "bold"),
         show_heatmap_legend = FALSE)
    dev.off()
    message(sprintf("    Saved %s", basename(png_path)))
  }
}

# -- Step 6: Save results CSV -------------------------------------------------

message("\nStep 6: Saving results CSV...")
results_dt <- rbindlist(all_results, use.names = TRUE)
results_dt[, log10_p := -log10(pmax(p_value, 1e-300))]

out_csv <- file.path(out_dir, "tabula_sapiens_permtest_results.csv")
fwrite(results_dt, out_csv)
message(sprintf("  Saved %s (%d rows)", basename(out_csv), nrow(results_dt)))

# -- Step 7: Weighted Specificity Enrichment (WSS) summary heatmap ------------
#
# For each cluster k and cell type c, compute:
#   WSS(k, c) = sum( W[i,k] * specificity[gene_i, c] )
# where i indexes all variants mapped to a gene present in the tissue.
# Permutation null: shuffle gene labels in the specificity matrix (10,000x).
# FDR: per-cluster BH correction at q < 0.05.

message("\nStep 7: Weighted Specificity Enrichment (WSS) summary heatmap...")

# Build variant-level table with W matrix weights and nearest gene
w_long <- w_mat %>%
  select(VAR_ID, all_of(k_cols)) %>%
  as.data.table()

# Merge with nearest gene mapping (one gene per variant)
var_gene <- nearest_genes[!is.na(gene_name), .(VAR_ID, gene_name)]
var_gene <- var_gene[!duplicated(VAR_ID)]
w_gene <- merge(w_long, var_gene, by = "VAR_ID")

message(sprintf("  %d variants with gene mappings and W weights", nrow(w_gene)))

# Nice label builder (reused for display)
format_celltype <- function(x) str_replace_all(x, "_", " ")

build_nice_label_map <- function(tissue_name, all_meta_dt, merge_rules,
                                 label_overrides, ts_data_dir) {
  expr_path <- file.path(ts_data_dir, sprintf("%s_expression.csv", tissue_name))
  header <- names(fread(expr_path, nrows = 0))
  cell_types_raw <- setdiff(header, "Gene_name")

  placeholder <- grepl("^not.?exist", cell_types_raw, ignore.case = TRUE)
  cell_types <- cell_types_raw[!placeholder]

  tissue_meta <- all_meta_dt[Tissue == tissue_name]
  if (nrow(tissue_meta) == length(cell_types)) {
    nice <- setNames(tissue_meta$Cell_Type, cell_types)
  } else {
    nice <- setNames(format_celltype(cell_types), cell_types)
  }

  for (col_name in names(label_overrides)) {
    if (col_name %in% names(nice)) {
      nice[col_name] <- label_overrides[col_name]
    }
  }

  if (tissue_name %in% names(merge_rules)) {
    for (rule in merge_rules[[tissue_name]]) {
      if (rule$keep %in% cell_types) {
        nice[rule$keep] <- rule$label
      }
      others <- setdiff(rule$cols, rule$keep)
      cell_types <- setdiff(cell_types, others)
    }
    nice <- nice[cell_types]
  }

  nice
}

# Process each tissue: compute WSS and permutation p-values
wss_results <- list()

for (tissue in tissues) {
  message(sprintf("  --- WSS for tissue: %s ---", tissue))

  spec_path <- file.path(ts_data_dir, sprintf("%s_specificity.csv", tissue))
  spec_dt <- fread(spec_path)

  cell_types_raw <- setdiff(names(spec_dt), "Gene_name")
  placeholder <- grepl("^not.?exist", cell_types_raw, ignore.case = TRUE)
  if (any(placeholder)) spec_dt[, (cell_types_raw[placeholder]) := NULL]
  cell_types <- cell_types_raw[!placeholder]

  # Apply merge rules to specificity data
  if (tissue %in% names(merge_rules)) {
    drop_cols <- character(0)
    for (rule in merge_rules[[tissue]]) {
      present <- rule$cols[rule$cols %in% cell_types]
      if (length(present) < 2) next
      avg <- rowMeans(as.matrix(spec_dt[, ..present]), na.rm = TRUE)
      spec_dt[, (rule$keep) := avg]
      drop_cols <- c(drop_cols, setdiff(present, rule$keep))
    }
    if (length(drop_cols) > 0) {
      drop_present <- drop_cols[drop_cols %in% names(spec_dt)]
      if (length(drop_present) > 0) spec_dt[, (drop_present) := NULL]
      cell_types <- setdiff(cell_types, drop_cols)
    }
  }

  # Subset to variants whose nearest gene is in this tissue's specificity data
  var_in_tissue <- w_gene[gene_name %in% spec_dt$Gene_name]
  if (nrow(var_in_tissue) < 2) {
    message(sprintf("    Skipping %s: too few variants with genes", tissue))
    next
  }

  # Build specificity matrix for mapped genes (row order matches var_in_tissue)
  gene_idx <- match(var_in_tissue$gene_name, spec_dt$Gene_name)
  spec_mat <- as.matrix(spec_dt[gene_idx, ..cell_types])
  n_variants <- nrow(spec_mat)
  n_celltypes <- ncol(spec_mat)

  # W matrix weights for these variants (n_variants x n_clusters)
  w_weights <- as.matrix(var_in_tissue[, ..k_cols])

  # Observed WSS: t(W) %*% spec_mat → (n_clusters x n_celltypes)
  obs_wss <- crossprod(w_weights, spec_mat)

  # All gene names in specificity data (pool for permutation)
  all_genes <- spec_dt$Gene_name
  n_all <- length(all_genes)
  all_spec_mat <- as.matrix(spec_dt[, ..cell_types])

  # Permutation test: shuffle gene assignments, recompute WSS
  exceed_count <- matrix(0L, nrow = length(k_cols), ncol = n_celltypes)
  for (p in seq_len(n_perm)) {
    perm_idx <- sample.int(n_all, n_variants)
    perm_spec <- all_spec_mat[perm_idx, , drop = FALSE]
    perm_wss <- crossprod(w_weights, perm_spec)
    exceed_count <- exceed_count + (perm_wss >= obs_wss)
  }
  p_mat <- exceed_count / n_perm

  # Store results
  for (ki in seq_along(k_cols)) {
    cl <- k_cols[ki]
    for (ci in seq_along(cell_types)) {
      wss_results[[length(wss_results) + 1]] <- data.table(
        tissue = tissue,
        cell_type = cell_types[ci],
        cluster = cl,
        wss_observed = obs_wss[ki, ci],
        p_value = p_mat[ki, ci]
      )
    }
  }

  message(sprintf("    %d variants x %d cell types x %d clusters",
                  n_variants, n_celltypes, length(k_cols)))
}

wss_dt <- rbindlist(wss_results)
wss_dt[, log10_p := -log10(pmax(p_value, 1 / (n_perm + 1)))]

# Per-cluster BH FDR correction
wss_dt[, fdr := p.adjust(p_value, method = "BH"), by = cluster]

# Report test counts and FDR results
for (cl in k_cols) {
  n_cl <- sum(wss_dt$cluster == cl)
  n_sig <- sum(wss_dt$cluster == cl & wss_dt$fdr < 0.05)
  message(sprintf("  %s (%s): %d tests, %d significant at FDR < 0.01",
                  cl, cluster_labels[cl], n_cl, n_sig))
}

# Save full WSS results
wss_csv <- file.path(out_dir, "tabula_sapiens_wss_results.csv")
fwrite(wss_dt, wss_csv)
message(sprintf("  Saved %s (%d rows)", basename(wss_csv), nrow(wss_dt)))

sig_wss <- wss_dt[fdr < 0.01]
message(sprintf("  Total FDR < 0.01: %d", nrow(sig_wss)))

if (nrow(sig_wss) == 0) {
  message("  No significant results — skipping summary heatmap.")
} else {
  # Build nice display labels
  nice_lookup <- list()
  for (tis in tissues) {
    lmap <- build_nice_label_map(tis, all_meta_dt, merge_rules,
                                 label_overrides, ts_data_dir)
    for (ct in names(lmap)) {
      nice_lookup[[paste(tis, ct, sep = "|")]] <- lmap[ct]
    }
  }

  format_tissue <- function(x) str_replace_all(x, "_", " ")

  sig_wss[, display_label := {
    key <- paste(tissue, cell_type, sep = "|")
    nice <- nice_lookup[key]
    ifelse(
      vapply(nice, is.null, logical(1)),
      paste0(format_tissue(tissue), " — ", format_celltype(cell_type)),
      paste0(format_tissue(tissue), " — ", unlist(nice))
    )
  }]

  # Build heatmap matrix
  sig_clusters <- sort(unique(sig_wss$cluster))
  sig_labels <- sort(unique(sig_wss$display_label))

  mat <- matrix(NA_real_, nrow = length(sig_labels), ncol = length(sig_clusters),
                dimnames = list(sig_labels, cluster_labels[sig_clusters]))

  for (i in seq_len(nrow(sig_wss))) {
    rl <- sig_wss$display_label[i]
    cl <- cluster_labels[sig_wss$cluster[i]]
    mat[rl, cl] <- sig_wss$log10_p[i]
  }

  # Tissue annotation color strip
  row_tissues <- str_extract(sig_labels, "^[^—]+") |> str_trim()
  unique_tissues <- sort(unique(row_tissues))
  tissue_colors <- setNames(
    hcl.colors(length(unique_tissues), palette = "Set 2"),
    unique_tissues
  )

  row_ha <- rowAnnotation(
    Tissue = row_tissues,
    col = list(Tissue = tissue_colors),
    show_legend = TRUE,
    annotation_legend_param = list(Tissue = list(title = "Tissue"))
  )

  # Color scale
  max_val <- ceiling(max(mat, na.rm = TRUE))
  max_val <- max(max_val, 1)
  summary_col_fun <- colorRamp2(
    seq(0, max_val, length.out = 10),
    colorRampPalette(c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"))(10)
  )

  mat[is.na(mat)] <- 0

  n_rows <- length(sig_labels)
  panel_h <- max(6, n_rows * 0.25 + 2)

  hm <- Heatmap(
    mat,
    col = summary_col_fun,
    na_col = "white",
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_order = sig_labels,
    column_order = cluster_labels[sig_clusters],
    row_names_side = "left",
    row_names_max_width = unit(12, "cm"),
    row_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    rect_gp = gpar(col = "white", lwd = 1),
    heatmap_legend_param = list(title = "-log10(p)", direction = "horizontal"),
    left_annotation = row_ha,
    name = "summary"
  )

  summary_path <- file.path(out_dir, "tabula_sapiens_summary_heatmap.png")
  png(summary_path, width = 12, height = panel_h, units = "in", res = 600)
  draw(hm,
       column_title = "Tabula Sapiens: FDR-Significant Cell Types (WSS, per-cluster FDR < 0.01)",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       heatmap_legend_side = "bottom")
  dev.off()
  message(sprintf("  Saved %s", basename(summary_path)))
}

# -- Step 8: Combined grid figures per tissue ----------------------------------

message("\nStep 8: Assembling combined grid figures for focus tissues...")

suppressPackageStartupMessages({
  library(plotgardener)
  library(png)
})

combined_tissues <- c("Blood", "Fat", "Heart", "Small_Intestine", "Liver", "Vasculature")
panel_w <- 6
label_h <- 0.6
gap_x   <- 0.25
gap_y   <- 0.15
top_margin <- 0.25

ts_place_panels <- function(cluster_imgs, clusters, grid_ncol, grid_nrow,
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
      fontcolor = "#D94801",
      x = x_pos + panel_w / 2, y = slot_y + label_h / 2,
      just = "center", default.units = "inches"
    )
  }

  dev.off()
  message(sprintf("  Saved %s", basename(out_path)))
}

for (tis in combined_tissues) {
  message(sprintf("--- Combined figures for %s ---", tis))

  imgs <- list()
  for (cl in names(cluster_labels)) {
    png_path <- file.path(out_dir, sprintf("tabula_sapiens_%s_%s.png", tis, cl))
    if (!file.exists(png_path)) {
      message(sprintf("  Warning: %s not found, skipping", basename(png_path)))
      next
    }
    imgs[[cl]] <- readPNG(png_path)
  }

  ref_img <- imgs[[names(imgs)[1]]]
  panel_h_px <- dim(ref_img)[1]
  panel_w_px <- dim(ref_img)[2]
  panel_h <- panel_h_px / (panel_w_px / panel_w)

  ts_place_panels(
    cluster_imgs = imgs,
    clusters = names(cluster_labels)[1:6],
    grid_ncol = 3, grid_nrow = 2,
    panel_w = panel_w, panel_h = panel_h,
    out_path = file.path(out_dir, sprintf("tabula_sapiens_combined_%s_part1.png", tis)),
    cluster_labels = cluster_labels,
    letter_offset = 0
  )

  ts_place_panels(
    cluster_imgs = imgs,
    clusters = names(cluster_labels)[7:10],
    grid_ncol = 2, grid_nrow = 2,
    panel_w = panel_w, panel_h = panel_h,
    out_path = file.path(out_dir, sprintf("tabula_sapiens_combined_%s_part2.png", tis)),
    cluster_labels = cluster_labels,
    letter_offset = 6
  )
}

message("\nDone.")
