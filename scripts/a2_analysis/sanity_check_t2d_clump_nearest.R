#!/usr/bin/env Rscript
#
# Sanity check (secondary): reproduce the paper's "nearest-to-hit" T2D gene list
# from Suzuki 2024 META summary statistics, using our own pipeline, and report
# concordance against the paper-provided list.
#
# Paper definition: clump with PLINK 1.9 + 1000G EUR LD panel, keep P<=5e-8,
# remove variants within 250kb correlated r2>=0.5 with the index (or P>=0.01);
# for each clump take the single nearest protein-coding gene to the index SNP.
#
# Reuse: ld_clump()/ensure_plink()/read_ref_bim() from prep_bnmf.R (a pure
# function library). ld_clump hardcodes --clump-p1 1 --clump-p2 1; pre-filtering
# to P<=5e-8 makes this identical to the paper's p1=5e-8/p2=0.01 (every survivor
# is < 0.01). Gene universe is parsed identically to a2_gtex_heatmaps.R Step 1
# (gencode.v19, protein-coding autosomes, ENSG version stripped) so the recovered
# genes are directly comparable to the heatmap input.
#
# Usage:
#   Rscript scripts/a2_analysis/sanity_check_t2d_clump_nearest.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"

# Reuse clumping helpers (pure function library; no top-level execution)
source(file.path(base_dir, "scripts/a1_analysis/prep_bnmf.R"))

# -- Parameters (reported up front) -------------------------------------------
sumstats_path <- file.path(base_dir,
  "sumstats/harmonized/Suzuki_Nature_2024.t2d.META.GRCh37.processed.txt.gz")
ref_panel     <- file.path(base_dir, "reference/1kg_eur/1000G.EUR.QC")  # per-chr bfiles
gtf_path      <- file.path(base_dir, "reference/gencode.v19.annotation.gtf.gz")
paper_path    <- file.path(base_dir,
  "results/a2_analysis/sanity_check/t2d_paper_nearest_genes.csv")
out_dir       <- file.path(base_dir, "results/a2_analysis/sanity_check")

p_threshold <- 5e-8   # index variant significance (paper)
clump_r2    <- 0.5    # paper r2
clump_kb    <- 250    # paper window (kb)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Sanity check: reproduce T2D nearest-to-hit genes (Suzuki META) ===\n")
cat(sprintf("  sumstats : %s\n", sumstats_path))
cat(sprintf("  LD panel : %s.<chr> (1000G EUR)\n", ref_panel))
cat(sprintf("  GTF      : %s (gencode.v19, protein-coding)\n", basename(gtf_path)))
cat(sprintf("  P<=%.0e, clump_r2=%g, clump_kb=%d\n", p_threshold, clump_r2, clump_kb))

# -- Step 1: Parse GTF for protein-coding genes (mirrors a2 Step 1) ------------
message("\nStep 1: Parsing GTF for protein-coding genes...")
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

# -- Step 2: Load sumstats, derive CHR/POS, filter to genome-wide sig ----------
message("Step 2: Loading Suzuki META sumstats and filtering P<=5e-8...")
ss <- fread(sumstats_path)
# VAR_ID is chr_pos_ref_alt; RSID column is unusable (mostly "."). Derive CHR/POS.
ss[, c("CHR", "POS") := tstrsplit(VAR_ID, "_", keep = 1:2)]
ss[, CHR := as.integer(CHR)]
ss[, POS := as.integer(POS)]
dt <- ss[!is.na(P_VALUE) & P_VALUE <= p_threshold & CHR %in% 1:22,
         .(VAR_ID, CHR, POS, P_VALUE)]
message(sprintf("  %d variants with P<=%.0e", nrow(dt), p_threshold))

# -- Step 3: LD clump (reuse ld_clump with paper r2/kb) -----------------------
message("Step 3: LD clumping (PLINK 1.9, 1000G EUR)...")
idx_varids <- ld_clump(dt, ref_panel_prefix = ref_panel,
                       clump_r2 = clump_r2, clump_kb = clump_kb)
idx <- dt[VAR_ID %in% idx_varids]
message(sprintf("  %d independent index variants", nrow(idx)))

# -- Step 4: Map each index variant to the single nearest protein-coding gene --
message("Step 4: Mapping index variants to single nearest protein-coding gene...")
nearest_for <- function(chr_i, pos_i) {
  g <- gene_dt[chr == as.character(chr_i)]
  if (nrow(g) == 0) return(list(gene_id = NA_character_, gene_name = NA_character_,
                                distance = NA_real_))
  d <- pmax(g$start - pos_i, 0) + pmax(pos_i - g$end, 0)  # same metric as a2
  j <- which.min(d)
  list(gene_id = g$gene_id[j], gene_name = g$gene_name[j], distance = d[j])
}
recovered <- idx[, {
  nn <- nearest_for(CHR, POS)
  list(gene_id = nn$gene_id, gene_name = nn$gene_name, distance = nn$distance)
}, by = .(VAR_ID, CHR, POS, P_VALUE)][order(P_VALUE)]

rec_csv <- file.path(out_dir, "t2d_clump_recovered_genes.csv")
fwrite(recovered, rec_csv)
message(sprintf("  Saved %s (%d index variants, %d unique genes)",
                basename(rec_csv), nrow(recovered), uniqueN(recovered$gene_id)))

# -- Step 5: Concordance vs paper-provided nearest-to-hit list -----------------
message("Step 5: Concordance vs paper list...")
paper <- fread(paper_path)
paper_ensg  <- unique(str_remove(paper$ensemblID, "\\.\\d+$"))
paper_names <- unique(paper$Gene_Name)
rec_ensg    <- unique(recovered$gene_id[!is.na(recovered$gene_id)])
rec_names   <- unique(recovered$gene_name[!is.na(recovered$gene_name)])

concordance <- data.table(
  Gene_Name   = union(paper_names, rec_names)
)[, in_paper := Gene_Name %in% paper_names
][, in_ours  := Gene_Name %in% rec_names
][, status   := fifelse(in_paper & in_ours, "recovered",
                 fifelse(in_paper & !in_ours, "missed (paper only)",
                         "extra (ours only)"))][order(-in_paper, Gene_Name)]
conc_csv <- file.path(out_dir, "t2d_clump_concordance_report.csv")
fwrite(concordance, conc_csv)

n_recov  <- sum(concordance$status == "recovered")
n_missed <- sum(concordance$status == "missed (paper only)")
n_extra  <- sum(concordance$status == "extra (ours only)")

cat("\n--- Concordance summary (by gene symbol) ---\n")
cat(sprintf("  paper unique genes : %d\n", length(paper_names)))
cat(sprintf("  our unique genes   : %d\n", length(rec_names)))
cat(sprintf("  recovered          : %d\n", n_recov))
cat(sprintf("  missed (paper only): %d  [%s]\n", n_missed,
            paste(concordance[status == "missed (paper only)", Gene_Name], collapse = ", ")))
cat(sprintf("  extra (ours only)  : %d  [%s]\n", n_extra,
            paste(concordance[status == "extra (ours only)", Gene_Name], collapse = ", ")))
cat(sprintf("  Saved %s\n", basename(conc_csv)))

cat("\nDone.\n")
