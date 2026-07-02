"""
column_map.py
Shared column-name mappings and header-reading utilities for GWAS processing scripts.

Maps source column names found in various GWAS file formats to gwaslab
Sumstats() parameter names. Used by both check_build.py and harmonize_sumstats.py.
"""

import gzip

import pandas as pd

# Maps source column names → list of gwaslab Sumstats() parameter(s).
# When a column maps to multiple params (e.g. SNP_ID → snpid + rsid),
# all are set unless already claimed by an earlier column.
COLUMN_MAP = {
    # SNP identifiers
    "SNP_ID":     ["snpid", "rsid"],
    "VAR_ID":     ["snpid"],
    "snpid":      ["snpid", "rsid"],   # Gurdasani (lowercase, compound CHR:POS:REF:ALT)
    "rsid":       ["rsid"],
    "rsID":       ["snpid", "rsid"],
    "RSID":       ["rsid"],
    "rs_id":      ["rsid"],
    "SNP":        ["snpid", "rsid"],  # Kanai (BOLT-LMM output)
    "ID":         ["snpid", "rsid"],  # PGC sumstats-VCF (MDD 2025, SCZ wave3): rsID column
    "MarkerName": ["snpid", "rsid"],  # Evangelou, Sinnott-Armstrong
    "snpname":    ["snpid", "rsid"],  # Justice NatGenet 2019 (WHR)
    "variant":    ["snpid"],          # NealeLab UKB (CHR:POS:REF:ALT format)
    # Chromosome
    "chrom":      ["chrom"],
    "Chromsome":  ["chrom"],
    "chromosome": ["chrom"],
    "CHROM":      ["chrom"],
    "#CHROM":     ["chrom"],    # VCF-style header (e.g. Sinnott-Armstrong)
    "#chrom":     ["chrom"],    # FinnGen
    "CHR":        ["chrom"],    # Sakaue (Biobank Japan SAIGE output)
    "chr":        ["chrom"],
    "Chr":        ["chrom"],      # Teumer uACR
    "CHROMOSOME":          ["chrom"],   # Morris NatGenet 2012
    "chromosome(b37)":     ["chrom"],   # DIAMANTE consortium (Mahajan 2022 EUR/EAS/SAS)
    # Position
    "pos":                ["pos"],
    "POS":                ["pos"],  # Sinnott-Armstrong AST
    "Position":           ["pos"],
    "base_pair_location": ["pos"],
    "BP":                 ["pos"],  # Kanai (BOLT-LMM output)
    "POS_b37":            ["pos"],
    "pos_b37":            ["pos"],
    "Pos":                ["pos"],  # Mahajan NatGenet 2022
    "Pos_b37":            ["pos"],  # Teumer uACR
    "POSITION":           ["pos"],  # Morris NatGenet 2012
    "position(b37)":      ["pos"],  # DIAMANTE consortium (Mahajan 2022 EUR/EAS/SAS)
    # Effect allele
    "ea":             ["ea"],
    "EffectAllele":   ["ea"],
    "Effect_Allele":  ["ea"],
    "effect_allele":  ["ea"],
    "ALT":            ["ea"],
    "alt":            ["ea"],       # Justice NatGenet 2019 (WHR)
    "A1":             ["ea"],
    "Allele1":        ["ea"],       # Evangelou, Teumer
    "ALLELE1":        ["ea"],       # Kanai (BOLT-LMM output)
    "EA":             ["ea"],       # Mahajan NatGenet 2022
    "EFFECT_ALLELE":  ["ea"],       # Mahajan NatGenet 2018 exome
    "RISK_ALLELE":    ["ea"],       # Morris NatGenet 2012
    "minor_allele":   ["ea"],       # NealeLab UKB
    # Non-effect allele
    "ref":              ["nea"],
    "NonEffectAllele":  ["nea"],
    "other_allele":     ["nea"],
    "REF":              ["nea"],
    "A2":               ["nea"],
    "Allele2":          ["nea"],    # Evangelou, Teumer
    "NEA":              ["nea"],    # Mahajan NatGenet 2022
    "OTHER_ALLELE":     ["nea"],    # Morris NatGenet 2012
    "ALLELE0":          ["nea"],    # Kanai (BOLT-LMM output)
    # Effect allele frequency
    "af":                       ["eaf"],
    "EAF":                      ["eaf"],
    "effect_allele_frequency":  ["eaf"],
    "POOLED_ALT_AF":            ["eaf"],
    "A1_freq":                  ["eaf"],
    "Freq1":                    ["eaf"],  # Evangelou, Teumer
    "A1FREQ":                   ["eaf"],  # Kanai (BOLT-LMM output)
    "minor_AF":                 ["eaf"],  # NealeLab UKB
    # Effect size
    "beta":        ["beta"],
    "Beta":        ["beta"],
    "BETA":        ["beta"],
    "EFFECT_SIZE": ["beta"],
    "A1_beta":     ["beta"],
    "METAL_Effect": ["beta"],
    "Effect":       ["beta"],      # Evangelou, Teumer, Sinnott-Armstrong
    "beta_fe":      ["beta"],      # Gurdasani (fixed-effects meta-analysis)
    "Fixed-effects_beta": ["beta"],  # DIAMANTE consortium (Mahajan 2022 EUR/EAS/SAS)
    "or":          ["OR"],
    "OR":          ["OR"],           # Morris NatGenet 2012
    "odds_ratio":  ["OR"],           # GWAS Catalog harmonized (Verma, Nikpay, Aragam)
    # Standard error
    "SE":             ["se"],
    "se":             ["se"],
    "sebeta":         ["se"],
    "METAL_StdErr":   ["se"],
    "StdErr":         ["se"],       # Evangelou, Teumer, Sinnott-Armstrong
    "se_fe":          ["se"],       # Gurdasani (fixed-effects meta-analysis)
    "standard_error": ["se"],
    "Fixed-effects_SE": ["se"],    # DIAMANTE consortium (Mahajan 2022 EUR/EAS/SAS)
    # P-value
    "pval":    ["p"],
    "Pval":    ["p"],
    "P_VALUE": ["p"],
    "p_value": ["p"],
    "pvalue":        ["p"],
    "Pvalue":        ["p"],        # Mahajan NatGenet 2022
    "METAL_Pvalue":  ["p"],
    "P":             ["p"],         # Evangelou
    "PVAL":          ["p"],         # PGC sumstats-VCF (MDD 2025, SCZ wave3)
    "P-value":       ["p"],         # Teumer, Sinnott-Armstrong
    "pval_fe":       ["p"],         # Gurdasani (fixed-effects meta-analysis)
    "p.value":       ["p"],         # Sakaue (SAIGE output)
    "P_BOLT_LMM_INF": ["p"],       # Kanai (BOLT-LMM output)
    "Fixed-effects_p-value": ["p"],  # DIAMANTE consortium (Mahajan 2022 EUR/EAS/SAS)
    "P_VALUE":       ["p"],          # Morris NatGenet 2012
    # Sample size
    "num_samples": ["n"],
    "Neff":        ["n"],
    "N_effective": ["n"],           # van Rheenen ALS 2021
    "NEFF":        ["n"],           # PGC MDD 2025 (effective N)
    "N":           ["n"],
    "n":           ["n"],
    "N_total":          ["n"],
    "TotalSampleSize":  ["n"],      # Evangelou
    "n_total_sum":      ["n"],      # Teumer
    "n_complete_samples": ["n"],    # NealeLab UKB
    "sample_size":      ["n"],      # Chen META (FastingGlucose/Insulin)
    "Weight":           ["n"],      # METAL output (Chen Diabetologia)
    "totalsamplesize":  ["n"],      # Graff AJHG 2023 (PAGE consortium)
    # SNP identifier — GWAS Catalog harmonized format (Loh, HuertaChagoya)
    "variant_id":       ["snpid"],
}


def read_header(filepath):
    """Read the first line of a file (handling .gz) and return column names and separator.

    Auto-detects delimiter: tab-separated first, falls back to space-separated
    (e.g. Evangelou, Teumer, Gurdasani files use spaces).

    Returns:
        tuple: (columns list, separator string)
    """
    opener = gzip.open if filepath.endswith(".gz") else open
    with opener(filepath, "rt") as f:
        header = f.readline().strip()
    # Auto-detect delimiter: tab first, fall back to space
    if "\t" in header:
        columns = header.split("\t")
        sep = "\t"
    else:
        columns = header.split()
        sep = r"\s+"
    return columns, sep


def split_compound_id(df):
    """Split compound variant IDs (CHR:POS:REF:ALT or CHR_POS_REF_ALT) into columns.

    Only acts when SNPID or rsID contains compound IDs but CHR/POS are missing.
    Also extracts REF→NEA and ALT→EA if those columns are missing.

    Args:
        df: gwaslab Sumstats .data DataFrame (modified in-place).

    Returns:
        bool: True if splitting was performed, False otherwise.
    """
    # Find the column containing the variant ID (gwaslab may rename to rsID)
    if "SNPID" in df.columns:
        id_col = "SNPID"
    elif "rsID" in df.columns:
        id_col = "rsID"
    else:
        return False
    if "CHR" in df.columns and "POS" in df.columns:
        return False

    # Detect separator from first non-null value
    non_null = df[id_col].dropna()
    if len(non_null) == 0:
        return False
    sample = str(non_null.iloc[0])
    if ":" in sample:
        sep = ":"
    elif "_" in sample:
        sep = "_"
    else:
        return False

    parts = df[id_col].str.split(sep, n=3, expand=True)
    if parts.shape[1] < 2:
        return False

    df["CHR"] = parts[0]
    df["POS"] = pd.to_numeric(parts[1], errors="coerce")

    # Extract REF→NEA and ALT→EA if available and not already present
    if parts.shape[1] >= 4:
        if "NEA" not in df.columns:
            df["NEA"] = parts[2]   # REF
        if "EA" not in df.columns:
            df["EA"] = parts[3]    # ALT

    print(f"  Extracted CHR/POS from compound {id_col} (sep='{sep}')", flush=True)
    return True


def detect_params(columns, build="99"):
    """Map file column names to gwaslab parameters using COLUMN_MAP.

    Args:
        columns: list of column name strings from the file header
        build: genome build to pass to gwaslab (default "99" = dummy for infer_build)

    Returns:
        dict of gwaslab Sumstats() keyword arguments
    """
    params = {}
    for col in columns:
        if col in COLUMN_MAP:
            for gwas_param in COLUMN_MAP[col]:
                if gwas_param not in params:
                    params[gwas_param] = col
    # Drop OR when beta is already mapped — having both causes gwaslab's
    # basic_check to fail when the OR column contains non-numeric values
    # (e.g. literal "NA" strings in Sakaue/GWAS Catalog files)
    if "beta" in params and "OR" in params:
        del params["OR"]

    params["build"] = str(build)
    return params
