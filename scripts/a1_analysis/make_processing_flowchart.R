#!/usr/bin/env Rscript
# make_processing_flowchart.R
# Regenerate the META bNMF GWAS-processing / variant-universe construction
# flowchart directly from the live config + QC result files, so it always
# reflects the current pipeline -- including the union LD-clump step.
#
# Emits a graphviz .dot source and renders it to PNG via the `dot` binary.
# All counts are read from result files, so the figure auto-corrects after a
# rerun. Steps whose result rows do not yet exist (e.g. union_ld_clump before
# the next rerun) are shown as "pending rerun".
#
# Usage:
#   Rscript scripts/a1_analysis/make_processing_flowchart.R
#   Rscript scripts/a1_analysis/make_processing_flowchart.R --config config/a1_config.yaml

suppressPackageStartupMessages({ library(yaml); library(data.table) })

# --- Args ---
args <- commandArgs(trailingOnly = TRUE)
config_path <- "config/a1_config.yaml"
if ("--config" %in% args) config_path <- args[which(args == "--config") + 1]
ancestry <- "META"

cfg <- read_yaml(config_path)
results_dir <- cfg$results_dir
anc_dir <- file.path(results_dir, ancestry)
out_dot <- file.path(results_dir, sprintf("%s_gwas_processing_flowchart.dot", ancestry))
out_png <- file.path(results_dir, sprintf("%s_gwas_processing_flowchart.png", ancestry))

# --- Pipeline parameters from config ---
ref_gwas <- names(unlist(cfg$ref_gwas[[ancestry]]))
n_ref <- length(ref_gwas)
n_cad <- sum(grepl("^CAD", ref_gwas, ignore.case = TRUE))
n_t2d <- sum(grepl("^T2D", ref_gwas, ignore.case = TRUE))
anc_tags <- paste(sort(unique(sub(".*_", "", ref_gwas))), collapse = " / ")  # study-group ancestries
p_thr    <- cfg$bnmf$p_threshold[[ancestry]]
maf      <- cfg$bnmf$maf_threshold
r2       <- cfg$ld_clump$r2
kb       <- cfg$ld_clump$kb
panel    <- basename(cfg$ld_clump$ref_panel[[ancestry]])
do_union <- ancestry %in% unlist(cfg$ld_clump$union_clump_ancestries)
miss     <- cfg$trait_missingness$missingness_threshold
proxy_r2 <- cfg$trait_missingness$proxy_r2_threshold
nreps    <- cfg$bnmf$nreps
K        <- cfg$bnmf$K
phi      <- cfg$bnmf$phi
hm3_on   <- !is.null(cfg$hapmap3_snp_file) && nzchar(cfg$hapmap3_snp_file)

# --- Counts from result files (graceful fallback to "pending rerun") ---
fmt <- function(x) if (length(x) == 0 || is.na(x)) "pending rerun" else format(as.integer(x), big.mark = ",")

qc_path <- file.path(anc_dir, sprintf("qc_report_%s.tsv", ancestry))
qc <- if (file.exists(qc_path)) fread(qc_path) else data.table(gwas_file = character(), step = character(), variants_remaining = integer())
ucount <- function(s) qc[gwas_file == "UNION" & step == s, variants_remaining]

n_combined <- fmt(ucount("union_combined"))
n_uhapmap3 <- fmt(ucount("union_hapmap3"))      # union-level HapMap3 (lean flow)
n_uclump   <- fmt(ucount("union_ld_clump"))     # absent until the union-clump rerun
n_dedup    <- fmt(ucount("union_dedup_varid"))  # standard flow only
n_multi    <- fmt(ucount("union_multiallelic_qc"))

fv_path <- file.path(anc_dir, sprintf("filtered_variants_%s.tsv", ancestry))
n_final <- if (file.exists(fv_path)) fmt(nrow(fread(fv_path, select = 1L))) else "pending rerun"
pm_path <- file.path(anc_dir, sprintf("prepared_matrix_%s.tsv", ancestry))
n_traitcol <- if (file.exists(pm_path)) fmt(ncol(fread(pm_path, nrows = 0)) - 1L) else "pending rerun"

# --- Build graphviz .dot ---
# \l = left-justified line break within a node label.
qc_steps <- paste(
  sprintf("1. P-value filter: P < %s", p_thr),
  "2. MAF filter: MAF >= 0.01",
  "3. SNPs only: drop indels",
  "4. Strand-ambiguous removal: drop A/T & C/G",
  "5. Multiallelic removal",
  "6. MHC exclusion: chr6 25-35 Mb",
  sprintf("7. LD clumping: r2 < %g, %d kb (%s)", r2, kb, panel),
  sep = "\\l")

graph_header <- sprintf('digraph META_gwas_processing {
  graph [rankdir=TB, fontname="Helvetica", labelloc=t, fontsize=15,
         label="META bNMF pipeline - GWAS processing & variant-universe construction\\nconfig/a1_config.yaml | scripts/a1_analysis/{01_run_bnmf.R, prep_bnmf.R} | build GRCh37\\nCounts reflect current result files; union-LD-clump and all downstream counts update after the next META rerun"];
  node [shape=box, style="filled,rounded", fontname="Helvetica", fontsize=11, margin="0.18,0.10"];
  edge [color="#666666"];

  input    [fillcolor="#cfe2f3", label="INPUT: %d reference GWAS = %d CAD + %d T2D\\nancestries: %s"];
  qc       [fillcolor="#fff2cc", label="PER-FILE QC (applied independently to each of the %d GWAS)\\l%s\\l"];',
  n_ref, n_cad, n_t2d, anc_tags, n_ref, qc_steps)

if (do_union) {
  # Lean union flow: stack -> HapMap3 -> LD clump across union.
  uclump_label <- sprintf(
    "NEW: LD clump across UNION\\l(r2 < %g, %d kb, %s)\\lranks index SNPs by min P across files\\l-> %s variants\\l",
    r2, kb, panel, n_uclump)
  body <- sprintf('
  uni      [fillcolor="#d9ead3", label="UNION of the %d per-file post-clump sets (stacked)\\l-> %s variants\\l"];
  hm3      [fillcolor="#d9ead3", label="HapMap3 restriction\\l-> %s variants\\l"];
  uclump   [fillcolor="#f9cb9c", color="#b45f06", penwidth=2, label="%s"];
  postuni  [fillcolor="#ead1dc", label="POST-UNION FILTERS\\ltrait-missingness filter (>%g%% missing,\\lproxy r2 >= %g substitution)\\l"];
  final    [fillcolor="#d9d2e9", label="FINAL VARIANT UNIVERSE: %s variants"];
  zmat     [fillcolor="#cfe2f3", label="Z-SCORE MATRIX: %s variants x %s trait columns"];
  bnmf     [fillcolor="#d5a6bd", label="bNMF FACTORIZATION\\nnreps = %d, K = %d, phi = %s"];

  input -> qc -> uni -> hm3 -> uclump -> postuni -> final -> zmat -> bnmf;
}',
    n_ref, n_combined, n_uhapmap3, uclump_label,
    miss * 100, proxy_r2, n_final, n_final, n_traitcol, nreps, K, phi)
} else {
  # Standard flow: VAR_ID dedup + multiallelic check; HapMap3 in post-union.
  hm3_line <- if (hm3_on) "HapMap3 restriction + " else ""
  body <- sprintf('
  uni      [fillcolor="#d9ead3", label="UNION of the %d per-file post-clump sets (%s variants stacked)\\ldedup VAR_ID -> %s unique\\lmultiallelic check (drop CHR:POS appearing >1x) -> %s\\l"];
  postuni  [fillcolor="#ead1dc", label="POST-UNION FILTERS\\l%strait-missingness filter (>%g%% missing,\\lproxy r2 >= %g substitution)\\l"];
  final    [fillcolor="#d9d2e9", label="FINAL VARIANT UNIVERSE: %s variants"];
  zmat     [fillcolor="#cfe2f3", label="Z-SCORE MATRIX: %s variants x %s trait columns"];
  bnmf     [fillcolor="#d5a6bd", label="bNMF FACTORIZATION\\nnreps = %d, K = %d, phi = %s"];

  input -> qc -> uni -> postuni -> final -> zmat -> bnmf;
}',
    n_ref, n_combined, n_dedup, n_multi,
    hm3_line, miss * 100, proxy_r2, n_final, n_final, n_traitcol, nreps, K, phi)
}

dot <- paste0(graph_header, body)

writeLines(dot, out_dot)
cat(sprintf("Wrote dot source: %s\n", out_dot))

# --- Render PNG via graphviz ---
dot_bin <- Sys.which("dot")
if (!nzchar(dot_bin)) stop("graphviz 'dot' not found on PATH; cannot render PNG.")
status <- system(sprintf("%s -Tpng %s -o %s", dot_bin, shQuote(out_dot), shQuote(out_png)))
if (status != 0) stop("dot failed to render the flowchart PNG.")
cat(sprintf("Rendered flowchart: %s\n", out_png))
