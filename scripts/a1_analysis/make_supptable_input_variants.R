#!/usr/bin/env Rscript
# make_supptable_input_variants.R
# Create supplementary table of 1,354 bNMF input variants with per-disease
# p-values and risk allele frequencies from CAD and T2D GWAS.

library(data.table)
library(openxlsx)
library(yaml)

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

cfg <- read_yaml("config/a1_config.yaml")
gtf_path <- cfg$gtf_path

# --- 1. Load the 1,354 variant IDs from the W matrix ---
w <- fread("results/a1_analysis/META/W_matrix_META.tsv", select = "VAR_ID")
target_ids <- w$VAR_ID
cat(sprintf("Target variants from W matrix: %d\n", length(target_ids)))

# --- 2. Get rsID and effect allele from filtered_variants ---
fv <- fread("results/a1_analysis/META/filtered_variants_META.tsv",
            select = c("VAR_ID", "RSID", "Effect_Allele"))
fv <- fv[VAR_ID %in% target_ids]
# Missing rsIDs are encoded as "." in the harmonized sumstats; leave the cell
# blank rather than printing a period as an empty-value placeholder.
fv[RSID == ".", RSID := NA_character_]

# --- 3. Parse VAR_ID into chr, pos, A1, A2 ---
fv[, c("chr", "pos", "A1", "A2") := {
  parts <- tstrsplit(VAR_ID, "_", fixed = TRUE)
  list(as.integer(parts[[1]]), as.integer(parts[[2]]), parts[[3]], parts[[4]])
}]

# --- 4. Look up min p-value and corresponding EAF from all CAD and T2D GWAS ---
ref_gwas <- cfg$ref_gwas$META
cad_files <- ref_gwas[grep("^CAD_", names(ref_gwas))]
t2d_files <- ref_gwas[grep("^T2D_", names(ref_gwas))]

lookup_min_pval <- function(gwas_files, target_ids, label) {
  results <- data.table()
  for (i in seq_along(gwas_files)) {
    fname <- gwas_files[[i]]
    cat(sprintf("  [%s %d/%d] Reading %s ...\n", label, i, length(gwas_files), basename(fname)))
    dt <- fread(fname, select = c("VAR_ID", "P_VALUE", "EAF"))
    dt <- dt[VAR_ID %in% target_ids]
    if (nrow(dt) > 0) {
      results <- rbind(results, dt)
    }
  }
  if (nrow(results) == 0) return(data.table(VAR_ID = character(), P_VALUE = numeric(), EAF = numeric()))
  results[, .SD[which.min(P_VALUE)], by = VAR_ID]
}

cat("\nLooking up CAD p-values...\n")
cad_lookup <- lookup_min_pval(cad_files, target_ids, "CAD")
setnames(cad_lookup, c("P_VALUE", "EAF"), c("P_value_CAD", "RAF_CAD"))

cat("\nLooking up T2D p-values...\n")
t2d_lookup <- lookup_min_pval(t2d_files, target_ids, "T2D")
setnames(t2d_lookup, c("P_VALUE", "EAF"), c("P_value_T2D", "RAF_T2D"))

# --- 5. Map to nearest gene (Locus) from GTF ---
cat("\nParsing GTF for gene annotations...\n")
gtf <- fread(cmd = paste("grep -v '^#'", gtf_path), sep = "\t", header = FALSE,
             select = c(1, 3, 4, 5, 9),
             col.names = c("chr", "feature", "start", "end", "attributes"))
gtf <- gtf[feature == "gene" & grepl('gene_biotype "protein_coding"', attributes)]
gtf[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", attributes)]
gtf[, chr := as.integer(sub("^chr", "", chr))]
gtf <- gtf[!is.na(chr), .(chr, start, end, gene_name)]
setkey(gtf, chr, start, end)

nearest_gene <- function(chrom, position, gene_dt) {
  genes_chr <- gene_dt[chr == chrom]
  if (nrow(genes_chr) == 0) return(NA_character_)
  inside <- genes_chr[start <= position & end >= position]
  if (nrow(inside) > 0) return(inside$gene_name[1])
  genes_chr[, dist := pmin(abs(start - position), abs(end - position))]
  genes_chr$gene_name[which.min(genes_chr$dist)]
}

cat("Mapping variants to nearest gene...\n")
fv[, Locus := mapply(nearest_gene, chr, pos, MoreArgs = list(gene_dt = gtf))]

# --- 6. Merge everything ---
result <- merge(fv, cad_lookup[, .(VAR_ID, P_value_CAD, RAF_CAD)], by = "VAR_ID", all.x = TRUE)
result <- merge(result, t2d_lookup[, .(VAR_ID, P_value_T2D, RAF_T2D)], by = "VAR_ID", all.x = TRUE)

setorder(result, chr, pos)
result[, Index := .I]

out <- result[, .(Index, `Variant ID (hg19)` = VAR_ID, `Risk Allele` = A1,
                  rsID = RSID, Locus,
                  `P-value_CAD` = P_value_CAD, `P-value_T2D` = P_value_T2D,
                  `RAF_CAD` = RAF_CAD, `RAF_T2D` = RAF_T2D)]

cat(sprintf("\nFinal table: %d rows, %d columns\n", nrow(out), ncol(out)))
cat(sprintf("  Variants with CAD p-value: %d\n", sum(!is.na(out$`P-value_CAD`))))
cat(sprintf("  Variants with T2D p-value: %d\n", sum(!is.na(out$`P-value_T2D`))))
cat(sprintf("  Variants with Locus: %d\n", sum(!is.na(out$Locus))))

# --- 7. Write Excel ---
out_dir <- "results/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "cad_t2d_input_variants.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Input Variants")

header_text <- "Supplementary Table X: CAD- and T2D-associated input variant set for clustering analysis"
header_style <- createStyle(textDecoration = "bold", fontSize = 12)
writeData(wb, 1, x = header_text, startRow = 1, startCol = 1)
addStyle(wb, 1, style = header_style, rows = 1, cols = 1)
mergeCells(wb, 1, cols = 1:ncol(out), rows = 1)

col_header_style <- createStyle(textDecoration = "bold", border = "Bottom",
                                halign = "center", wrapText = TRUE)
writeData(wb, 1, x = out, startRow = 3, headerStyle = col_header_style)

pval_style <- createStyle(numFmt = "0.00E+00")
addStyle(wb, 1, style = pval_style,
         rows = 4:(nrow(out) + 3), cols = which(names(out) == "P-value_CAD"),
         gridExpand = TRUE)
addStyle(wb, 1, style = pval_style,
         rows = 4:(nrow(out) + 3), cols = which(names(out) == "P-value_T2D"),
         gridExpand = TRUE)

raf_style <- createStyle(numFmt = "0.000")
addStyle(wb, 1, style = raf_style,
         rows = 4:(nrow(out) + 3), cols = which(names(out) == "RAF_CAD"),
         gridExpand = TRUE)
addStyle(wb, 1, style = raf_style,
         rows = 4:(nrow(out) + 3), cols = which(names(out) == "RAF_T2D"),
         gridExpand = TRUE)

# Build the reference-study lists dynamically from the config so the legend
# always reflects the actual CAD/T2D GWAS used (currently a 4-CAD / 5-T2D subset).
study_names <- function(gwas_named) {
  s <- sub(".*?([A-Za-z]+)([0-9]{4}).*", "\\1 \\2", names(gwas_named))
  sort(unique(s))
}
cad_studies <- study_names(cad_files)
t2d_studies <- study_names(t2d_files)

legend_start <- nrow(out) + 5
legend_lines <- c(
  "Legend:",
  "Variant ID (hg19) = chr_position_A1_A2 in GRCh37/hg19 coordinates.",
  "Risk Allele = effect allele (A1) from harmonized GWAS summary statistics.",
  "rsID = dbSNP reference SNP identifier; left blank where no rsID is available.",
  "Locus = nearest protein-coding gene (Ensembl GRCh37.75 GTF annotation).",
  sprintf("P-value_CAD = minimum P-value across the CAD reference GWAS (%s; multiple ancestries).",
          paste(cad_studies, collapse = ", ")),
  sprintf("P-value_T2D = minimum P-value across the T2D reference GWAS (%s; multiple ancestries).",
          paste(t2d_studies, collapse = ", ")),
  "RAF_CAD = risk allele frequency from the CAD GWAS with the minimum P-value.",
  "RAF_T2D = risk allele frequency from the T2D GWAS with the minimum P-value.",
  "Blank P-value/RAF = variant not present in any GWAS for that disease."
)
legend_style <- createStyle(fontSize = 10, wrapText = TRUE)
bold_legend <- createStyle(textDecoration = "bold", fontSize = 10)

for (i in seq_along(legend_lines)) {
  writeData(wb, 1, x = legend_lines[i], startRow = legend_start + i - 1, startCol = 1)
  sty <- if (i == 1) bold_legend else legend_style
  addStyle(wb, 1, style = sty, rows = legend_start + i - 1, cols = 1)
  mergeCells(wb, 1, cols = 1:ncol(out), rows = legend_start + i - 1)
}

setColWidths(wb, 1, cols = 1:ncol(out),
             widths = c(8, 22, 12, 14, 14, 14, 14, 10, 10))

saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("\nSaved: %s\n", out_file))
