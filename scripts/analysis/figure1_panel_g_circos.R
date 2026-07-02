# figure1_panel_g_circos.R
# Panel G: Circos plot of gene-level contributions to each bNMF cluster.
# Maps SNPs from W_matrix_{ANC}.tsv to genes via GTF, aggregates mean normalized
# weight per gene, and displays one concentric track per ancestry x cluster.

# Assumes figure_utils.R has been sourced (ggplot2, tidyverse already loaded).


#' Parse a GTF file to extract protein-coding gene coordinates
#'
#' @param gtf_path Path to an Ensembl GTF file (GRCh37)
#' @return tibble with columns: chr, start, end, gene_name
parse_gtf_genes <- function(gtf_path) {
  gtf <- read_tsv(
    gtf_path, comment = "#",
    col_names = c("chr", "source", "feature", "start", "end",
                  "score", "strand", "frame", "attributes"),
    col_types = cols(
      chr = col_character(), start = col_double(), end = col_double(),
      .default = col_character()
    )
  )

  gtf %>%
    filter(feature == "gene", str_detect(attributes, 'gene_biotype "protein_coding"')) %>%
    mutate(gene_name = str_replace(attributes, '.*gene_name "([^"]+)".*', "\\1")) %>%
    filter(chr %in% c(as.character(1:22), "X", "Y")) %>%
    select(chr, start, end, gene_name)
}


#' Map SNPs to genes using gene-body overlap
#'
#' @param snp_df tibble with columns: chr, pos, VAR_ID, plus cluster weight columns
#' @param gene_df tibble from parse_gtf_genes()
#' @return tibble with gene_name joined to SNPs (SNPs not in any gene are dropped)
map_snps_to_genes <- function(snp_df, gene_df) {
  snp_df %>%
    inner_join(gene_df, by = "chr", relationship = "many-to-many") %>%
    filter(pos >= start, pos <= end)
}


#' Create Panel G: Grid of individual Circos plots per ancestry x cluster
#'
#' Each cluster gets its own minimalistic circos plot with chromosome ideogram
#' and green bars showing average gene weight (zero-weight genes excluded).
#'
#' @param results_dir Path to results directory
#' @param ancestries Character vector of ancestry codes
#' @param gtf_path Path to GRCh37 GTF file
#' @return patchwork object combining all individual circos plots
plot_panel_g <- function(results_dir = "results",
                         ancestries = c("EUR", "AFR", "META"),
                         gtf_path = "/sc/arion/projects/paul_oreilly/lab/kestir01/Homo_sapiens.GRCh37.75.gtf") {

  # --- Step 1: Parse GTF ---
  cat("  Parsing GTF for protein-coding genes...\n")
  gene_df <- parse_gtf_genes(gtf_path)
  cat(sprintf("  Found %d protein-coding genes\n", nrow(gene_df)))

  # --- Step 2 & 3: Map SNPs to genes, aggregate per gene x cluster ---
  track_data_list <- list()

  for (anc in ancestries) {
    w_path <- file.path(results_dir, anc, sprintf("W_matrix_%s.tsv", anc))
    if (!file.exists(w_path)) {
      warning(sprintf("W_matrix_%s.tsv not found for %s", anc, anc))
      next
    }

    w_df <- read_tsv(w_path, show_col_types = FALSE) %>%
      mutate(
        chr = str_replace(VAR_ID, "^([^_]+)_.*", "\\1"),
        pos = as.numeric(str_replace(VAR_ID, "^[^_]+_([^_]+)_.*", "\\1"))
      )
    cluster_cols <- setdiff(colnames(w_df), c("VAR_ID", "chr", "pos"))

    # Map to genes
    mapped <- map_snps_to_genes(w_df, gene_df)
    cat(sprintf("  %s: %d SNPs mapped to genes (of %d total)\n",
                anc, n_distinct(mapped$VAR_ID), nrow(w_df)))

    # Aggregate mean weight per gene per cluster
    for (kcol in cluster_cols) {
      gene_weights <- mapped %>%
        group_by(gene_name) %>%
        summarise(
          mean_weight = mean(.data[[kcol]], na.rm = TRUE),
          chr = first(chr), start = first(start), end = first(end),
          .groups = "drop"
        ) %>%
        filter(mean_weight > 0) %>%
        slice_max(mean_weight, n = 10) %>%
        mutate(norm_weight = if_else(max(mean_weight) > 0,
                                     mean_weight / max(mean_weight), 0))

      track_label <- paste(anc, kcol)
      track_data_list[[track_label]] <- list(
        data  = gene_weights,
        label = track_label
      )
    }
  }

  if (length(track_data_list) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data for Circos",
                 size = 6) +
        theme_big_text()
    )
  }

  # --- Step 4: Build one horizontal bar chart per ancestry x cluster ---
  bar_color <- "#50B878"

  bar_plots <- map(track_data_list, function(track_info) {
    td <- track_info$data
    label <- track_info$label

    # Order by weight ascending so highest is at top after coord_flip
    td <- td %>%
      arrange(norm_weight) %>%
      mutate(gene_name = fct_inorder(gene_name))

    ggplot(td, aes(x = gene_name, y = norm_weight)) +
      geom_col(fill = bar_color, width = 0.7) +
      coord_flip() +
      theme_big_text() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ggtitle(label) +
      labs(x = NULL, y = "Normalized Weight")
  })

  # --- Step 5: Combine into grid (as single element for one panel tag) ---
  ncols <- 2
  patchwork::wrap_elements(patchwork::wrap_plots(bar_plots, ncol = ncols))
}
