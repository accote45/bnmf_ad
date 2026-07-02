#!/usr/bin/env python3
"""
extract_gwas_catalog.py
Extract and standardize GWAS Catalog top-hit TSV files into harmonized format.

These files contain only a handful of SNPs in NHGRI-EBI GWAS Catalog format.
Parses: CHR_ID, CHR_POS, SNPS (rsID), STRONGEST SNP-RISK ALLELE, OR or BETA,
95% CI (TEXT), RISK ALLELE FREQUENCY, P-VALUE, INITIAL SAMPLE SIZE.

Output columns match harmonize_sumstats.py:
    VAR_ID, RSID, Effect_Allele, P_VALUE, BETA, SE, N, MAF, EAF

Usage:
    python extract_gwas_catalog.py
"""

import os
import re

import numpy as np
import pandas as pd

PROJECT_ROOT = "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic"
SUMSTATS = os.path.join(PROJECT_ROOT, "sumstats")
HARMONIZED = os.path.join(SUMSTATS, "harmonized")

# GWAS Catalog top-hit files to process
CATALOG_FILES = [
    "Flannick_Nature_2019.t2d.META.GRCh37.tsv",
    "Mercader_Diabetes_2017.t2d.AMR.GRCh37.tsv",
    "Ng_PlosGenetics_2014.t2d.AFR.GRCh37.tsv",
    "Sigma_Nature_2014.t2d.AMR.GRCh37.tsv",
]


def extract_allele(snp_risk_allele):
    """Extract risk allele from 'rsXXXXX-A' format."""
    if pd.isna(snp_risk_allele):
        return np.nan
    parts = str(snp_risk_allele).rsplit("-", 1)
    if len(parts) == 2 and parts[1].strip():
        allele = parts[1].strip().upper()
        if allele == "?":
            return np.nan
        return allele
    return np.nan


def parse_or_beta(value):
    """Parse 'OR or BETA' column — convert OR to log(OR) for BETA."""
    if pd.isna(value):
        return np.nan
    val = float(value)
    if val <= 0:
        return np.nan
    # Values > 0 and < 5 are likely ORs; take log
    # Values that are already on log scale (negative or very small) would be unusual
    # for GWAS Catalog which reports ORs for binary traits
    return np.log(val)


def parse_se_from_ci(ci_text, beta):
    """Derive SE from '95% CI (TEXT)' field and BETA.

    CI format is typically '[lower-upper]' or 'lower-upper'.
    SE ≈ (log(upper) - log(lower)) / (2 * 1.96) for OR-based CIs.
    """
    if pd.isna(ci_text) or pd.isna(beta):
        return np.nan
    ci_str = str(ci_text).strip().strip("[]() ")
    # Try to parse 'lower-upper' pattern (handle negative signs carefully)
    match = re.match(r"([\d.]+)\s*[-–]\s*([\d.]+)", ci_str)
    if not match:
        return np.nan
    try:
        lower = float(match.group(1))
        upper = float(match.group(2))
        if lower <= 0 or upper <= 0:
            return np.nan
        # CI is on OR scale; convert to log scale
        se = (np.log(upper) - np.log(lower)) / (2 * 1.96)
        return se if se > 0 else np.nan
    except (ValueError, ZeroDivisionError):
        return np.nan


def parse_sample_size(text):
    """Extract total N from 'INITIAL SAMPLE SIZE' text via regex.

    Sums all numbers found in patterns like '2,277 African American cases'.
    """
    if pd.isna(text):
        return np.nan
    # Find all numbers (with optional commas) that are followed by a word
    numbers = re.findall(r"(\d[\d,]*)\s+\w+", str(text))
    if not numbers:
        return np.nan
    total = sum(int(n.replace(",", "")) for n in numbers)
    return total if total > 0 else np.nan


def make_var_id(chrom, pos, allele):
    """Create a minimal VAR_ID. With only one allele, use allele + '.' as sorted pair."""
    c = str(chrom)
    p = str(pos)
    a = str(allele).upper() if pd.notna(allele) else "."
    # Single-allele: can't sort a pair, use allele_. as placeholder
    a1, a2 = sorted([".", a])
    return f"{c}_{p}_{a1}_{a2}"


def process_catalog_file(filename):
    """Process a single GWAS Catalog TSV into harmonized format."""
    filepath = os.path.join(SUMSTATS, filename)
    print(f"\nProcessing: {filename}", flush=True)

    df = pd.read_csv(filepath, sep="\t")
    print(f"  Rows: {len(df)}", flush=True)

    # Parse columns
    out = pd.DataFrame()
    out["CHR"] = df["CHR_ID"].astype(str)
    out["POS"] = pd.to_numeric(df["CHR_POS"], errors="coerce")
    out["RSID"] = df["SNPS"].fillna(".")
    out["Effect_Allele"] = df["STRONGEST SNP-RISK ALLELE"].apply(extract_allele)
    out["P_VALUE"] = pd.to_numeric(df["P-VALUE"], errors="coerce")
    out["BETA"] = df["OR or BETA"].apply(parse_or_beta)

    # SE from 95% CI
    if "95% CI (TEXT)" in df.columns:
        out["SE"] = [
            parse_se_from_ci(ci, beta)
            for ci, beta in zip(df["95% CI (TEXT)"], out["BETA"])
        ]
    else:
        out["SE"] = np.nan

    # EAF
    if "RISK ALLELE FREQUENCY" in df.columns:
        out["EAF"] = pd.to_numeric(df["RISK ALLELE FREQUENCY"], errors="coerce")
    else:
        out["EAF"] = np.nan

    # N from INITIAL SAMPLE SIZE
    if "INITIAL SAMPLE SIZE" in df.columns:
        out["N"] = df["INITIAL SAMPLE SIZE"].apply(parse_sample_size)
    else:
        out["N"] = np.nan

    # MAF
    out["MAF"] = np.where(
        out["EAF"].notna(),
        np.minimum(out["EAF"].astype(float), 1 - out["EAF"].astype(float)),
        np.nan,
    )

    # VAR_ID
    out["VAR_ID"] = [
        make_var_id(c, p, a)
        for c, p, a in zip(out["CHR"], out["POS"], out["Effect_Allele"])
    ]

    # Drop rows missing critical values
    before = len(out)
    out = out.dropna(subset=["CHR", "POS", "BETA", "P_VALUE"])
    dropped = before - len(out)
    if dropped > 0:
        print(f"  Dropped {dropped} rows with missing critical values", flush=True)

    # Select output columns
    out_cols = ["VAR_ID", "RSID", "Effect_Allele", "P_VALUE", "BETA", "SE", "N", "MAF", "EAF"]
    out = out[out_cols]

    # Save
    base = filename.replace(".tsv", "")
    out_file = os.path.join(HARMONIZED, f"{base}.processed.txt.gz")
    os.makedirs(HARMONIZED, exist_ok=True)
    out.to_csv(out_file, sep="\t", index=False, compression="gzip")
    print(f"  Saved {len(out)} variants to {out_file}", flush=True)

    return out_file


def main():
    print("=== GWAS Catalog Top-Hit Extraction ===", flush=True)
    for filename in CATALOG_FILES:
        filepath = os.path.join(SUMSTATS, filename)
        if not os.path.exists(filepath):
            print(f"\nSKIP (not found): {filename}", flush=True)
            continue
        process_catalog_file(filename)
    print("\n=== Done ===", flush=True)


if __name__ == "__main__":
    main()
