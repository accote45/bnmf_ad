# generate_synthetic_data.R
# Creates realistic but fake GWAS summary statistics for testing the bNMF pipeline.
# Outputs: one CAD GWAS file + 5 trait files per ancestry, written to tests/test_data/

library(data.table)

generate_synthetic_gwas <- function(output_dir = "tests/test_data",
                                    n_variants = 500,
                                    n_significant = 50,
                                    seed = 42) {
  set.seed(seed)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  alleles <- c("A", "T", "C", "G")
  chromosomes <- 1:22

  # Generate variant IDs in CHR_POS_REF_ALT format
  chrs <- sample(chromosomes, n_variants, replace = TRUE)
  positions <- sort(sample(1e5:1e8, n_variants, replace = FALSE))
  ref_alleles <- sample(alleles, n_variants, replace = TRUE)
  alt_alleles <- sapply(ref_alleles, function(r) sample(setdiff(alleles, r), 1))
  var_ids <- paste(chrs, positions, ref_alleles, alt_alleles, sep = "_")

  # Helper: generate one GWAS sumstats file
  make_gwas <- function(var_ids, n_variants, n_significant, trait_name, sample_size) {
    betas <- rnorm(n_variants, mean = 0, sd = 0.02)
    ses <- runif(n_variants, min = 0.005, max = 0.03)
    ses <- pmax(ses, 0.005)  # ensure positive

    # Make some variants significant by giving them large effects
    sig_idx <- sample(seq_len(n_variants), n_significant)
    betas[sig_idx] <- rnorm(n_significant, mean = 0, sd = 0.15)
    ses[sig_idx] <- runif(n_significant, min = 0.005, max = 0.02)

    p_values <- 2 * pnorm(-abs(betas / ses))

    maf <- runif(n_variants, min = 0.005, max = 0.5)
    # Add a few low-MAF variants that should be filtered
    low_maf_idx <- sample(setdiff(seq_len(n_variants), sig_idx), 5)
    maf[low_maf_idx] <- runif(5, min = 0.0001, max = 0.0009)

    n_ph <- rep(sample_size, n_variants) + sample(-500:500, n_variants, replace = TRUE)

    data.table(
      VAR_ID = var_ids,
      Effect_Allele = alt_alleles,
      P_VALUE = p_values,
      BETA = betas,
      SE = ses,
      N = n_ph,
      MAF = maf
    )
  }

  # Trait configurations: trait_name -> sample_size
  trait_configs <- list(
    T2D    = 200000,
    LDL    = 350000,
    Stroke = 150000,
    BMI    = 500000,
    HF     = 100000
  )
  cad_sample_size <- 250000

  ancestries <- c("EUR", "AFR")

  for (anc in ancestries) {
    # CAD GWAS (variant backbone) - more significant variants
    cad <- make_gwas(var_ids, n_variants, n_significant = n_significant, "CAD", cad_sample_size)
    cad_file <- file.path(output_dir, paste0("synthetic_CAD_", anc, ".txt.gz"))
    fwrite(cad, cad_file, sep = "\t")

    # Trait GWAS files - each has partially overlapping significant variants
    for (trait_name in names(trait_configs)) {
      # Vary the number of significant variants per trait
      n_sig_trait <- max(10, n_significant + sample(-20:20, 1))
      trait_dt <- make_gwas(var_ids, n_variants, n_significant = n_sig_trait,
                            trait_name, trait_configs[[trait_name]])
      trait_file <- file.path(output_dir, paste0("synthetic_", trait_name, "_", anc, ".txt.gz"))
      fwrite(trait_dt, trait_file, sep = "\t")
    }
  }

  cat("Synthetic data generated in:", output_dir, "\n")
  cat("Files created:\n")
  cat(paste(" ", list.files(output_dir, pattern = "synthetic_"), collapse = "\n"), "\n")

  invisible(output_dir)
}

# Run if executed directly
if (!interactive() && identical(sys.frame(1)$ofile, NULL) ||
    (!interactive() && length(sys.frames()) == 0)) {
  generate_synthetic_gwas()
}
