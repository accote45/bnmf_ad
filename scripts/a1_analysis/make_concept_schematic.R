#!/usr/bin/env Rscript
# make_concept_schematic.R
# Presentation schematic of the AD bNMF subtyping pipeline (non-technical).
#   AD GWAS loci  +  trait GWAS panel  ->  Z-score matrix  ->  bNMF  ->  4 subtypes
# Pure diagram (no data inputs). Output: figures/concept_schematic.png
#
# Usage: Rscript scripts/a1_analysis/make_concept_schematic.R [out.png]

suppressPackageStartupMessages({ library(ggplot2) })

out <- commandArgs(trailingOnly = TRUE)
out <- if (length(out)) out[1] else "figures/concept_schematic.png"
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito colorblind-safe palette for the 4 subtypes
col_imm <- "#E69F00"  # Immune  (orange)
col_lip <- "#009E73"  # Lipid   (green)
col_met <- "#0072B2"  # Metabolic (blue)
col_neu <- "#CC79A7"  # Neurodegeneration (purple)
ink     <- "#222222"
soft    <- "#5a5a5a"

# --- box + arrow helpers ---------------------------------------------------
box <- function(xc, yc, w, h, fill, border = ink, r = 0.12) {
  # simple rounded-ish rectangle via a filled rect (ggplot has no native rounding)
  annotate("rect", xmin = xc - w/2, xmax = xc + w/2, ymin = yc - h/2, ymax = yc + h/2,
           fill = fill, color = border, linewidth = 0.6)
}
txt <- function(xc, yc, label, size = 4.2, col = ink, face = "plain", lineheight = 0.95) {
  annotate("text", x = xc, y = yc, label = label, size = size, color = col,
           fontface = face, lineheight = lineheight)
}
arrow_seg <- function(x0, y0, x1, y1) {
  annotate("segment", x = x0, y = y0, xend = x1, yend = y1,
           linewidth = 0.9, color = soft,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"))
}

p <- ggplot() +
  # ---- inputs -------------------------------------------------------------
  box(2.4, 9.1, 3.4, 1.5, "#EDEDED") +
  txt(2.4, 9.45, "AD GWAS", size = 5, face = "bold") +
  txt(2.4, 8.9, "105 genome-wide loci\n(APOE region excluded)", size = 3.6, col = soft) +

  box(7.6, 9.1, 3.7, 1.5, "#EDEDED") +
  txt(7.6, 9.45, "21 trait GWAS", size = 5, face = "bold") +
  txt(7.6, 8.9, "lipids, immune, metabolic,\nneuro, psychiatric", size = 3.6, col = soft) +

  arrow_seg(2.4, 8.35, 4.3, 7.35) +
  arrow_seg(7.6, 8.35, 5.7, 7.35) +

  # ---- Z matrix -----------------------------------------------------------
  box(5.0, 6.6, 5.2, 1.4, "#DCE6F1") +
  txt(5.0, 6.9, "Variant  x  Trait  Z-score matrix", size = 4.6, face = "bold") +
  txt(5.0, 6.4, "105 loci   x   21 traits", size = 3.8, col = soft) +
  arrow_seg(5.0, 5.85, 5.0, 5.15) +

  # ---- bNMF ---------------------------------------------------------------
  box(5.0, 4.5, 3.2, 1.1, "#CBB7D8") +
  txt(5.0, 4.72, "bNMF", size = 5, face = "bold") +
  txt(5.0, 4.28, "Bayesian matrix factorization", size = 3.4, col = soft) +
  arrow_seg(5.0, 3.9, 5.0, 3.25) +

  txt(5.0, 3.02, "Four genetic subtypes of Alzheimer's disease", size = 4.6, face = "bold") +

  # ---- 4 subtype outputs --------------------------------------------------
  box(1.9, 1.6, 2.05, 1.5, col_imm) +
  txt(1.9, 1.95, "Immune", size = 4.3, face = "bold", col = "white") +
  txt(1.9, 1.35, "CD33\nSCIMP", size = 3.5, col = "white") +

  box(4.05, 1.6, 2.05, 1.5, col_lip) +
  txt(4.05, 1.95, "Lipid", size = 4.3, face = "bold", col = "white") +
  txt(4.05, 1.35, "ABCA1", size = 3.5, col = "white") +

  box(6.2, 1.6, 2.05, 1.5, col_met) +
  txt(6.2, 1.95, "Metabolic", size = 4.0, face = "bold", col = "white") +
  txt(6.2, 1.35, "SPI1\nKAT8", size = 3.5, col = "white") +

  box(8.5, 1.6, 2.4, 1.5, col_neu) +
  txt(8.5, 1.95, "Neurodegen.", size = 3.8, face = "bold", col = "white") +
  txt(8.5, 1.35, "GRN, MAPT,\nTMEM106B", size = 3.3, col = "white") +

  coord_cartesian(xlim = c(0.3, 9.9), ylim = c(0.6, 10.1), expand = FALSE) +
  theme_void()

ggsave(out, p, width = 10, height = 7, dpi = 300, bg = "white")
cat("Wrote", out, "\n")
