"""
Snakemake pipeline for GWAS harmonization, bNMF analysis, and figure assembly.

Steps:
  0. gwas_extractor — (optional) download raw GWAS files from Google Sheet
  1. check_build     — infer actual genome build of each raw GWAS file (-> JSON)
  2. harmonize       — standardize columns, liftover, allele alignment to T2D, QC report
  3. run_bnmf        — run bNMF factorization per ancestry
  4. assemble_figure — produce multi-panel figure from all ancestries

Prerequisites:
  - Raw GWAS files in sumstats/ (downloaded via gwas_extractor.py or manually)
  - Python environment with gwaslab, pandas, numpy, scipy
  - R 4.2.0 with data.table, tidyverse, patchwork, pheatmap

Usage:
  snakemake -n                                       # dry run
  snakemake --cores 32                               # full run (local)
  snakemake --profile profiles/lsf                   # HPC with LSF submission
  snakemake --profile profiles/lsf -n                # HPC dry run
  snakemake --profile profiles/lsf -R harmonize      # re-run harmonization only
  snakemake gwas_extractor                           # download raw GWAS (optional)
"""

import glob
import os
import re

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SUMSTATS_DIR   = config["sumstats_dir"]
HARMONIZED_DIR = config["harmonized_dir"]
BUILD_CHECK_DIR = config["build_check_dir"]
PREFERRED_BUILD = str(config["preferred_build"])
ANCESTRIES     = config["ancestries"]
N_OVERRIDES    = config.get("n_overrides", {})
T2D_REFERENCE  = config.get("t2d_reference", {})
QC_DIR         = config.get("qc_dir", "sumstats/qc")
SKIP_PREFIXES  = config.get("skip_files", [])

# bNMF parameters
RESULTS_DIR      = config.get("results_dir", "results")
FIGURES_DIR      = config.get("figures_dir", "results/figures")
BNMF_CONFIG      = config.get("bnmf", {})
BNMF_NREPS       = BNMF_CONFIG.get("nreps", 10)
BNMF_K           = BNMF_CONFIG.get("K", 20)
BNMF_PHI         = BNMF_CONFIG.get("phi", 1.0)
BNMF_P_STRICT    = BNMF_CONFIG.get("p_threshold_strict", 5e-8)
BNMF_P_RELAXED   = BNMF_CONFIG.get("p_threshold_relaxed", 1e-5)
BNMF_STRICT_ANCS = BNMF_CONFIG.get("strict_ancestries", ["EUR", "META"])


def get_p_threshold(ancestry):
    """Return p-value threshold based on ancestry."""
    return BNMF_P_STRICT if ancestry in BNMF_STRICT_ANCS else BNMF_P_RELAXED


def get_ancestry(base):
    """Extract ancestry from a file base string like '...CAD.EUR.GRCh37'."""
    for ancestry in ANCESTRIES:
        if f".{ancestry}." in base:
            return ancestry
    return None


# Map numeric build -> GRCh label for output filenames
BUILD_LABEL = {"19": "GRCh37", "38": "GRCh38"}
PREFERRED_LABEL = BUILD_LABEL.get(PREFERRED_BUILD, f"GRCh{PREFERRED_BUILD}")

SCRIPTS_DIR = "scripts/gwas_processing"
DATA_EXTS = (".txt", ".txt.gz", ".tsv", ".tsv.gz")

# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def get_file_base(filename):
    """Extract base name through GRCh##.

    E.g. 'Tcheandjieu_NatureMed_2023.CAD.AFR.GRCh37.txt.gz'
      -> 'Tcheandjieu_NatureMed_2023.CAD.AFR.GRCh37'
    """
    parts = filename.split(".")
    for i, part in enumerate(parts):
        if re.match(r"GRCh\d+", part):
            return ".".join(parts[: i + 1])
    return parts[0]


def rewrite_build_label(base, target_label):
    """Replace GRCh## in a file base with the target label.

    E.g. ('...AFR.GRCh38', 'GRCh37') -> '...AFR.GRCh37'
    """
    return re.sub(r"GRCh\d+", target_label, base)


# Build set of T2D reference basenames to exclude from processing
T2D_REF_BASES = set()
for anc, ref_path in T2D_REFERENCE.items():
    ref_fname = os.path.basename(ref_path)
    T2D_REF_BASES.add(get_file_base(ref_fname))

# Discover raw GWAS files and build lookup tables
# RAW_LOOKUP: output_base -> raw filename (for input lookup functions)
RAW_LOOKUP = {}       # output_base -> raw_filename
RAW_BASE_MAP = {}     # raw_base -> raw_filename
OUTPUT_BASES = []     # list of output bases (with preferred build label)

for ancestry in ANCESTRIES:
    pattern = os.path.join(SUMSTATS_DIR, f"*.{ancestry}.*")
    for filepath in sorted(glob.glob(pattern)):
        fname = os.path.basename(filepath)
        if not any(fname.endswith(ext) for ext in DATA_EXTS):
            continue
        raw_base = get_file_base(fname)
        out_base = rewrite_build_label(raw_base, PREFERRED_LABEL)

        # Skip T2D reference files (static, not reprocessed)
        if raw_base in T2D_REF_BASES or out_base in T2D_REF_BASES:
            continue

        # Skip files listed in config skip_files (e.g. no usable beta/SE)
        if any(raw_base.startswith(prefix) for prefix in SKIP_PREFIXES):
            continue

        RAW_LOOKUP[out_base] = fname
        RAW_BASE_MAP[raw_base] = fname
        OUTPUT_BASES.append(out_base)

# Discover harmonized files per ancestry (for bNMF inputs).
# Uses existing files if present, otherwise falls back to Snakemake output paths
# so the dependency chain (harmonize -> run_bnmf) is properly established.
HARMONIZED_BY_ANCESTRY = {}
for ancestry in ANCESTRIES:
    pattern = os.path.join(HARMONIZED_DIR, f"*.{ancestry}.*")
    files = sorted(glob.glob(pattern))
    existing = [f for f in files if any(f.endswith(ext) for ext in DATA_EXTS)]
    if existing:
        HARMONIZED_BY_ANCESTRY[ancestry] = existing
    else:
        # No harmonized files on disk yet -- depend on Snakemake outputs
        HARMONIZED_BY_ANCESTRY[ancestry] = [
            os.path.join(HARMONIZED_DIR, f"{base}.processed.txt.gz")
            for base in OUTPUT_BASES
            if f".{ancestry}." in base
        ]


def get_n_override(raw_filename):
    """Look up N override for a raw filename from config, or return empty string."""
    for key, value in N_OVERRIDES.items():
        if raw_filename.startswith(key):
            return str(value)
    return ""


def get_ref_gwas_flag(wildcards):
    """Return --ref-gwas flag for ancestry-matched T2D reference, or empty string."""
    ancestry = get_ancestry(wildcards.base)
    if ancestry and ancestry in T2D_REFERENCE:
        return f"--ref-gwas {T2D_REFERENCE[ancestry]}"
    return ""


# Lightweight rules run on the head node; heavy rules are submitted to LSF
localrules: all, gwas_extractor, qc_report, assemble_figure


# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------

rule all:
    input:
        expand(
            os.path.join(HARMONIZED_DIR, "{base}.processed.txt.gz"),
            base=OUTPUT_BASES,
        ),
        os.path.join(QC_DIR, "harmonization_qc_report.csv"),
        os.path.join(FIGURES_DIR, "figure_main.pdf"),
        os.path.join(FIGURES_DIR, "figure_main.png"),


# ---------------------------------------------------------------------------
# Rule: gwas_extractor (optional, not in default target)
# ---------------------------------------------------------------------------

rule gwas_extractor:
    """Download GWAS files from Google Sheet. Run explicitly: snakemake gwas_extractor"""
    output:
        touch("sumstats/.gwas_extractor_done"),
    params:
        scripts=SCRIPTS_DIR,
        sumstats_dir=SUMSTATS_DIR,
    shell:
        """
        python {params.scripts}/gwas_extractor.py {params.sumstats_dir}
        """


# ---------------------------------------------------------------------------
# Rule: check_build
# ---------------------------------------------------------------------------

rule check_build:
    """Infer the actual genome build of a raw GWAS file."""
    input:
        gwas=lambda wildcards: os.path.join(
            SUMSTATS_DIR, RAW_BASE_MAP[wildcards.raw_base]
        ),
    output:
        build_json=os.path.join(BUILD_CHECK_DIR, "{raw_base}.build.json"),
    params:
        scripts=SCRIPTS_DIR,
    threads: 1
    resources:
        mem_mb=32000,
        runtime=30,
    shell:
        """
        OPENBLAS_NUM_THREADS=1 python {params.scripts}/check_build.py \
            --input-file {input.gwas} \
            --output-file {output.build_json}
        """


# ---------------------------------------------------------------------------
# Rule: harmonize
# ---------------------------------------------------------------------------

def harmonize_input(wildcards):
    """Resolve inputs for the harmonize rule from the output base."""
    out_base = wildcards.base
    raw_fname = RAW_LOOKUP[out_base]
    raw_base = get_file_base(raw_fname)
    return {
        "gwas": os.path.join(SUMSTATS_DIR, raw_fname),
        "build_json": os.path.join(BUILD_CHECK_DIR, f"{raw_base}.build.json"),
    }


rule harmonize:
    """Harmonize a single GWAS file: column detection, liftover, OR->BETA, allele alignment."""
    input:
        unpack(harmonize_input),
    output:
        processed=os.path.join(HARMONIZED_DIR, "{base}.processed.txt.gz"),
        qc_json=os.path.join(QC_DIR, "{base}.qc.json"),
    params:
        scripts=SCRIPTS_DIR,
        preferred_build=PREFERRED_BUILD,
        n_override_flag=lambda wildcards: (
            f"--n-override {get_n_override(RAW_LOOKUP[wildcards.base])}"
            if get_n_override(RAW_LOOKUP[wildcards.base])
            else ""
        ),
        ref_gwas_flag=lambda wildcards: get_ref_gwas_flag(wildcards),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        """
        OPENBLAS_NUM_THREADS=1 python {params.scripts}/harmonize_sumstats.py \
            --input-file {input.gwas} \
            --build-info {input.build_json} \
            --output-file {output.processed} \
            --qc-json {output.qc_json} \
            --preferred-build {params.preferred_build} \
            --n-cores {threads} \
            {params.n_override_flag} \
            {params.ref_gwas_flag}
        """


# ---------------------------------------------------------------------------
# Rule: qc_report
# ---------------------------------------------------------------------------

rule qc_report:
    """Aggregate per-file QC JSONs into a single CSV report."""
    input:
        qc_jsons=expand(
            os.path.join(QC_DIR, "{base}.qc.json"),
            base=OUTPUT_BASES,
        ),
    output:
        csv=os.path.join(QC_DIR, "harmonization_qc_report.csv"),
    run:
        import json
        import csv
        rows = []
        for qc_file in input.qc_jsons:
            with open(qc_file) as f:
                rows.append(json.load(f))
        if rows:
            fieldnames = list(rows[0].keys())
            with open(output.csv, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows)


# ---------------------------------------------------------------------------
# Rule: run_bnmf
# ---------------------------------------------------------------------------

rule run_bnmf:
    """Run bNMF analysis for a single ancestry."""
    input:
        harmonized=lambda wildcards: HARMONIZED_BY_ANCESTRY[wildcards.ancestry],
    output:
        h_matrix    = os.path.join(RESULTS_DIR, "{ancestry}", "H_matrix.tsv"),
        w_matrix    = os.path.join(RESULTS_DIR, "{ancestry}", "W_matrix.tsv"),
        filtered    = os.path.join(RESULTS_DIR, "{ancestry}", "filtered_variants_{ancestry}.tsv"),
        prepared    = os.path.join(RESULTS_DIR, "{ancestry}", "prepared_matrix.tsv"),
        heatmap_w   = os.path.join(RESULTS_DIR, "{ancestry}", "heatmap_W.pdf"),
        heatmap_h   = os.path.join(RESULTS_DIR, "{ancestry}", "heatmap_H.pdf"),
        run_summary = os.path.join(RESULTS_DIR, "{ancestry}", "run_summary.csv"),
        summary_txt = os.path.join(RESULTS_DIR, "{ancestry}", "summary.txt"),
    params:
        nreps       = BNMF_NREPS,
        K           = BNMF_K,
        p_threshold = lambda wildcards: get_p_threshold(wildcards.ancestry),
        phi         = BNMF_PHI,
    threads: 1
    resources:
        mem_mb=32000,
        runtime=720,
    shell:
        """
        module load R/4.2.0
        export OMP_NUM_THREADS=1
        Rscript scripts/pipeline/run_toy_pipeline.R \
            --ancestry {wildcards.ancestry} \
            --nreps {params.nreps} \
            --K {params.K} \
            --p-threshold {params.p_threshold} \
            --phi {params.phi}
        """


# ---------------------------------------------------------------------------
# Rule: assemble_figure
# ---------------------------------------------------------------------------

rule assemble_figure:
    """Assemble multi-panel figure from bNMF results across all ancestries."""
    input:
        h_matrices=expand(
            os.path.join(RESULTS_DIR, "{ancestry}", "H_matrix.tsv"),
            ancestry=ANCESTRIES,
        ),
    output:
        pdf = os.path.join(FIGURES_DIR, "figure_main.pdf"),
        png = os.path.join(FIGURES_DIR, "figure_main.png"),
    params:
        results_dir = RESULTS_DIR,
        ancestries  = ",".join(ANCESTRIES),
    shell:
        """
        module load R/4.2.0
        Rscript scripts/analysis/assemble_figure.R \
            --results-dir {params.results_dir} \
            --ancestries {params.ancestries}
        """
