#!/usr/bin/env Rscript
# make_supptable_trait_gwas.R
# Create supplementary table listing all trait GWAS used in the META bNMF
# clustering analysis.

library(data.table)
library(openxlsx)
library(yaml)

base_dir <- "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
setwd(base_dir)

cfg <- read_yaml("config/a1_config.yaml")

# --- 1. Get the traits actually used from the prepared matrix columns ---
cols <- names(fread("results/a1_analysis/META/prepared_matrix_META.tsv", nrows = 0))
cols <- cols[cols != "VAR_ID"]
traits_used <- unique(sub("_(pos|neg)$", "", cols))
cat(sprintf("Traits in prepared matrix: %d\n", length(traits_used)))

# --- 2. Match each trait to its GWAS file from config (META then EUR fallback) ---
meta_gwas <- cfg$trait_gwas$META
eur_gwas  <- cfg$trait_gwas$EUR

file_map <- data.table(
  trait_key = character(),
  file_path = character(),
  source_ancestry = character()
)

for (tk in traits_used) {
  if (tk %in% names(meta_gwas)) {
    file_map <- rbind(file_map, data.table(trait_key = tk, file_path = meta_gwas[[tk]],
                                           source_ancestry = "META"))
  } else if (tk %in% names(eur_gwas)) {
    file_map <- rbind(file_map, data.table(trait_key = tk, file_path = eur_gwas[[tk]],
                                           source_ancestry = "EUR"))
  } else {
    file_map <- rbind(file_map, data.table(trait_key = tk, file_path = NA_character_,
                                           source_ancestry = NA_character_))
    cat(sprintf("  WARNING: no config entry for %s\n", tk))
  }
}

# Extract ancestry from filename (e.g., "...ALP.META.GRCh37..." -> "META")
file_map[, ancestry_from_file := sub(".*\\.([A-Z]+)\\.GRCh37.*", "\\1", file_path)]

# --- 3. Map trait categories ---
trait_categories <- list(
  Lipids = c("HDL", "LDL", "TotalCholesterol", "Triglycerides", "ApoA", "ApoB", "Lpa"),
  Anthropometric = c("BMI", "Height", "WaistCircum", "WaistHipRatio", "BMR"),
  `Body Composition` = c("VAT", "ASAT", "GFAT", "VATASAT", "VATGFAT", "ASATGFAT"),
  `Blood Pressure` = c("SBP", "DBP", "MAP", "HR"),
  Glycemic = c("FastingGlucose", "FastingInsulin", "RandomGlucose", "Hba1c", "HOMAB", "HOMAIR"),
  `Liver / Hepatic` = c("ALT", "AST", "GGT", "ALP", "Bilirubin"),
  Inflammatory = c("CRP"),
  Renal = c("uACR", "Creatinine", "CystatinC", "Urate", "Urea", "Phosphate"),
  Hematologic = c("Hgb", "WBC", "RBC", "MCV", "LymphocyteCount", "MonocyteCount",
                   "NeutrophilCount", "PltCount", "MeanPltVol", "PltDistWidth",
                   "RBCDistWidth", "Eosinophilcount", "MeanReticVol",
                   "ReticulocyteCount", "Basophilcount"),
  `Endocrine / Hormonal` = c("Testosterone", "SHBG", "IGF1"),
  `Mineral / Vitamin` = c("Calcium", "VitD", "Protein"),
  Adipokine = c("Adiponectin"),
  Atherosclerosis = c("atherosclerosis")
)

cat_lookup <- rbindlist(lapply(names(trait_categories), function(cat) {
  data.table(trait_base = trait_categories[[cat]], category = cat)
}))

# Parse trait base name (everything before the last _AuthorYear)
file_map[, trait_base := sub("_[A-Za-z]+\\d{4}$", "", trait_key)]

# Clean trait name for display (expand abbreviations)
trait_display <- c(
  HDL = "HDL cholesterol", LDL = "LDL cholesterol",
  TotalCholesterol = "Total cholesterol", Triglycerides = "Triglycerides",
  ApoA = "Apolipoprotein A", ApoB = "Apolipoprotein B", Lpa = "Lipoprotein(a)",
  BMI = "Body mass index", Height = "Height",
  WaistCircum = "Waist circumference", WaistHipRatio = "Waist-hip ratio",
  BMR = "Basal metabolic rate",
  VAT = "Visceral adipose tissue", ASAT = "Abdominal subcutaneous adipose tissue",
  GFAT = "Gluteofemoral adipose tissue",
  VATASAT = "VAT/ASAT ratio", VATGFAT = "VAT/GFAT ratio", ASATGFAT = "ASAT/GFAT ratio",
  SBP = "Systolic blood pressure", DBP = "Diastolic blood pressure",
  MAP = "Mean arterial pressure", HR = "Heart rate",
  FastingGlucose = "Fasting glucose", FastingInsulin = "Fasting insulin",
  RandomGlucose = "Random glucose", Hba1c = "HbA1c",
  HOMAB = "HOMA-B", HOMAIR = "HOMA-IR",
  ALT = "Alanine aminotransferase", AST = "Aspartate aminotransferase",
  GGT = "Gamma-glutamyl transferase", ALP = "Alkaline phosphatase",
  Bilirubin = "Bilirubin",
  CRP = "C-reactive protein",
  uACR = "Urinary albumin-to-creatinine ratio", Creatinine = "Creatinine",
  CystatinC = "Cystatin C", Urate = "Urate", Urea = "Urea", Phosphate = "Phosphate",
  Hgb = "Hemoglobin", WBC = "White blood cell count", RBC = "Red blood cell count",
  MCV = "Mean corpuscular volume",
  LymphocyteCount = "Lymphocyte count", MonocyteCount = "Monocyte count",
  NeutrophilCount = "Neutrophil count", PltCount = "Platelet count",
  MeanPltVol = "Mean platelet volume", PltDistWidth = "Platelet distribution width",
  RBCDistWidth = "Red cell distribution width",
  Eosinophilcount = "Eosinophil count", MeanReticVol = "Mean reticulocyte volume",
  ReticulocyteCount = "Reticulocyte count", Basophilcount = "Basophil count",
  Testosterone = "Testosterone", SHBG = "Sex hormone-binding globulin",
  IGF1 = "Insulin-like growth factor 1",
  Calcium = "Calcium", VitD = "Vitamin D", Protein = "Total protein",
  Adiponectin = "Adiponectin",
  atherosclerosis = "Atherosclerosis"
)

file_map <- merge(file_map, cat_lookup, by = "trait_base", all.x = TRUE)
file_map[, trait_display := trait_display[trait_base]]

# Order by category then trait
cat_order <- names(trait_categories)
file_map[, category := factor(category, levels = cat_order)]
setorder(file_map, category, trait_display)

file_map[, filename := basename(file_path)]

out <- file_map[, .(`Trait Category` = as.character(category),
                     Trait = trait_display,
                     `Sample Size (N)` = NA_character_,
                     `Ancestry Group` = ancestry_from_file,
                     Reference = NA_character_,
                     `File Name` = filename)]

cat(sprintf("\nFinal table: %d rows\n", nrow(out)))
cat(sprintf("  Ancestry breakdown: %s\n",
            paste(names(table(out$`Ancestry Group`)),
                  table(out$`Ancestry Group`), sep = "=", collapse = ", ")))

# --- 4. Write Excel ---
out_dir <- "results/supplementary_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "trait_gwas_used.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Trait GWAS")

header_text <- "Supplementary Table X: Trait GWAS summary statistics used for clustering analysis"
header_style <- createStyle(textDecoration = "bold", fontSize = 12)
writeData(wb, 1, x = header_text, startRow = 1, startCol = 1)
addStyle(wb, 1, style = header_style, rows = 1, cols = 1)
mergeCells(wb, 1, cols = 1:ncol(out), rows = 1)

col_header_style <- createStyle(textDecoration = "bold", border = "Bottom",
                                halign = "center", wrapText = TRUE)
writeData(wb, 1, x = out, startRow = 3, headerStyle = col_header_style)

setColWidths(wb, 1, cols = 1:ncol(out), widths = c(22, 38, 16, 18, 30, 55))

legend_start <- nrow(out) + 5
legend_lines <- c(
  "Legend:",
  "Trait Category = functional grouping of the measured biomarker or phenotype.",
  "Trait = full name of the measured trait.",
  "Sample Size (N) = total number of individuals in the GWAS.",
  "Ancestry Group = population ancestry of the GWAS cohort (META = multi-ancestry meta-analysis; EUR = European).",
  "Reference = citation for the GWAS publication.",
  "File Name = harmonized summary statistics file name used in the analysis.",
  "",
  "Abbreviations:",
  "HDL = high-density lipoprotein; LDL = low-density lipoprotein; HbA1c = glycated hemoglobin; HOMA-B = homeostatic model assessment of beta-cell function; HOMA-IR = homeostatic model assessment of insulin resistance; ALT = alanine aminotransferase; AST = aspartate aminotransferase; GGT = gamma-glutamyl transferase; ALP = alkaline phosphatase; CRP = C-reactive protein; uACR = urinary albumin-to-creatinine ratio; VAT = visceral adipose tissue; ASAT = abdominal subcutaneous adipose tissue; GFAT = gluteofemoral adipose tissue; SHBG = sex hormone-binding globulin; IGF-1 = insulin-like growth factor 1; MCV = mean corpuscular volume; RBC = red blood cell; WBC = white blood cell."
)
legend_style <- createStyle(fontSize = 10, wrapText = TRUE)
bold_legend <- createStyle(textDecoration = "bold", fontSize = 10)

for (i in seq_along(legend_lines)) {
  writeData(wb, 1, x = legend_lines[i], startRow = legend_start + i - 1, startCol = 1)
  sty <- if (i %in% c(1, 8)) bold_legend else legend_style
  addStyle(wb, 1, style = sty, rows = legend_start + i - 1, cols = 1)
  mergeCells(wb, 1, cols = 1:ncol(out), rows = legend_start + i - 1)
}

saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("\nSaved: %s\n", out_file))
