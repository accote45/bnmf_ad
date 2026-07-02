#!/usr/bin/env python3
"""
harmonize_sumstats.py
Harmonize a single GWAS summary statistics file for the bNMF pipeline.

Reads a raw GWAS file, auto-detects column mappings, runs QC, performs
conditional liftover based on build info, converts OR to BETA if needed,
optionally aligns effect alleles to a reference GWAS (e.g., T2D),
and saves a standardized output file with a QC metrics JSON.

Output columns: VAR_ID, RSID, Effect_Allele, P_VALUE, BETA, SE, N, MAF, EAF

Usage:
    OPENBLAS_NUM_THREADS=1 python harmonize_sumstats.py \\
        --input-file sumstats/Tcheandjieu_NatureMed_2023.CAD.AFR.GRCh37.txt.gz \\
        --build-info sumstats/build_checks/Tcheandjieu_NatureMed_2023.CAD.AFR.GRCh37.build.json \\
        --output-file sumstats/harmonized/Tcheandjieu_NatureMed_2023.CAD.AFR.GRCh37.processed.txt.gz \\
        --preferred-build 19 \\
        --n-cores 4 \\
        --ref-gwas sumstats/harmonized/Suzuki_Nature_2024.t2d.AFR.GRCh37.txt.gz \\
        --qc-json sumstats/qc/Tcheandjieu_NatureMed_2023.CAD.AFR.GRCh37.qc.json
"""

import argparse
import json
import os
import time

import gwaslab as gl
import numpy as np
import pandas as pd
import scipy.stats as ss

from column_map import read_header, detect_params, split_compound_id


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_var_id(df):
    """Create VAR_ID with sorted alleles for cross-file matching.

    Uses vectorized string operations for memory efficiency on large DataFrames
    (50M+ rows). Alleles are sorted alphabetically so the same variant matches
    across files regardless of which allele is labeled as effect/non-effect.
    """
    chr_s = df["CHR"].astype(str)
    pos_s = df["POS"].astype(str)
    ea_s = df["EA"].astype(str)
    nea_s = df["NEA"].astype(str)

    # Sort alleles: a1 = min, a2 = max (alphabetical)
    mask = ea_s <= nea_s
    a1 = ea_s.where(mask, nea_s)
    a2 = nea_s.where(mask, ea_s)

    return chr_s + "_" + pos_s + "_" + a1 + "_" + a2


def derive_se_from_ci(sumstats_data, input_file, columns, sep):
    """Derive SE from confidence interval bounds when SE is missing.

    When a file has OR with CI bounds (ci_upper, ci_lower) but SE is all NA,
    compute SE on the log(OR) scale: SE = (ln(CI_upper) - ln(CI_lower)) / 3.92.
    This is preferred over P-value-based derivation as it preserves precision
    for extreme p-values.

    Args:
        sumstats_data: gwaslab Sumstats .data DataFrame (modified in-place).
        input_file: path to raw GWAS file (to read CI columns not loaded by gwaslab).
        columns: list of column names from file header.
        sep: file separator string.

    Returns:
        bool: True if SE was derived from CI, False otherwise.
    """
    if "SE" in sumstats_data.columns:
        return False
    if "OR" not in sumstats_data.columns:
        return False
    if "ci_upper" not in columns or "ci_lower" not in columns:
        return False

    # Read CI columns directly from raw file (not loaded by gwaslab)
    read_sep = sep if sep == "\t" else None
    ci_df = pd.read_csv(input_file, sep=read_sep, usecols=["ci_upper", "ci_lower"])

    if len(ci_df) != len(sumstats_data):
        print(f"  WARNING: CI row count ({len(ci_df)}) != data row count "
              f"({len(sumstats_data)}), skipping CI-based SE derivation", flush=True)
        return False

    ci_upper = pd.to_numeric(ci_df["ci_upper"], errors="coerce")
    ci_lower = pd.to_numeric(ci_df["ci_lower"], errors="coerce")

    # SE(log(OR)) = (ln(CI_upper) - ln(CI_lower)) / (2 * 1.96)
    se = (np.log(ci_upper) - np.log(ci_lower)) / 3.92
    sumstats_data["SE"] = se.values
    n_valid = sumstats_data["SE"].notna().sum()
    print(f"  Derived SE from CI bounds (ci_upper, ci_lower) for {n_valid} variants", flush=True)
    return True


def convert_or_to_beta(df):
    """Convert OR to BETA via log(OR) and derive SE from BETA and P.

    Only runs when the DataFrame has an OR column but no BETA column.
    SE is derived as |BETA / Z| where Z = sign(beta) * |qnorm(p/2)|.
    Returns True if conversion was performed, False otherwise.
    """
    if "OR" not in df.columns or "BETA" in df.columns:
        return False

    df["BETA"] = np.log(df["OR"].astype(float))
    print(f"  Converted OR to BETA via log(OR)", flush=True)

    if "SE" not in df.columns and "P" in df.columns:
        beta = df["BETA"].values.astype(float)
        p = df["P"].values.astype(float)
        z = np.sign(beta) * np.abs(ss.norm.ppf(p / 2))
        z_safe = np.where(np.abs(z) > 1e-300, z, np.nan)
        df["SE"] = np.abs(beta / z_safe)
        print(f"  Derived SE from BETA and P for {df['SE'].notna().sum()} variants", flush=True)

    return True


def postprocess(sumstats):
    """Standardize columns, create VAR_ID. Returns processed DataFrame."""
    df = sumstats.data.copy()
    print(f"  gwaslab columns: {list(df.columns)}", flush=True)
    print(f"  Rows after QC: {len(df)}", flush=True)

    # Create VAR_ID
    print("  Creating VAR_ID (vectorized)...", flush=True)
    df["VAR_ID"] = make_var_id(df)

    # Effect allele
    df["Effect_Allele"] = df["EA"]

    # Rename to pipeline format
    rename_map = {}
    if "P" in df.columns:
        rename_map["P"] = "P_VALUE"
    if "rsID" in df.columns:
        rename_map["rsID"] = "RSID"
    elif "SNPID" in df.columns:
        rename_map["SNPID"] = "RSID"
    df = df.rename(columns=rename_map)

    if "RSID" not in df.columns:
        df["RSID"] = "."

    # Compute MAF
    if "EAF" in df.columns:
        df["MAF"] = np.minimum(df["EAF"].astype(float), 1 - df["EAF"].astype(float))
    else:
        df["MAF"] = np.nan

    # Select output columns
    out_cols = ["VAR_ID", "RSID", "Effect_Allele", "P_VALUE", "BETA", "SE", "N", "MAF", "EAF"]
    for col in out_cols:
        if col not in df.columns:
            print(f"  WARNING: missing column {col}, filling with NA", flush=True)
            df[col] = np.nan
    df = df[out_cols].copy()

    return df


def align_alleles_to_reference(df, ref_file):
    """Align effect alleles to a reference GWAS (e.g., T2D).

    For overlapping variants (matched by VAR_ID), if the Effect_Allele
    differs from the reference, flip BETA sign, swap Effect_Allele,
    and invert EAF.

    Args:
        df: DataFrame with VAR_ID, Effect_Allele, BETA, EAF columns.
        ref_file: Path to harmonized reference GWAS (tab-separated,
            with VAR_ID and Effect_Allele columns).

    Returns:
        int: Number of alleles flipped.
    """
    ref = pd.read_csv(ref_file, sep="\t", usecols=["VAR_ID", "Effect_Allele"])
    ref = ref.drop_duplicates(subset="VAR_ID").set_index("VAR_ID")

    overlap_mask = df["VAR_ID"].isin(ref.index)
    overlapping = df.loc[overlap_mask]

    if len(overlapping) == 0:
        print(f"  Allele alignment: 0 overlapping variants with reference", flush=True)
        return 0

    ref_alleles = ref.loc[overlapping["VAR_ID"].values, "Effect_Allele"].values
    study_alleles = overlapping["Effect_Allele"].values

    flip_mask = study_alleles != ref_alleles
    n_flipped = int(flip_mask.sum())

    if n_flipped > 0:
        flip_idx = overlapping.index[flip_mask]
        df.loc[flip_idx, "BETA"] = -df.loc[flip_idx, "BETA"]
        df.loc[flip_idx, "Effect_Allele"] = ref.loc[
            df.loc[flip_idx, "VAR_ID"].values, "Effect_Allele"
        ].values
        # Flip EAF: frequency of new effect allele = 1 - old frequency
        if "EAF" in df.columns:
            df.loc[flip_idx, "EAF"] = 1.0 - df.loc[flip_idx, "EAF"].astype(float)
            # Recompute MAF
            df.loc[flip_idx, "MAF"] = np.minimum(
                df.loc[flip_idx, "EAF"].astype(float),
                1 - df.loc[flip_idx, "EAF"].astype(float),
            )

    print(
        f"  Allele alignment: {len(overlapping)} overlapping variants, "
        f"{n_flipped} alleles flipped to match reference",
        flush=True,
    )
    return n_flipped


# ---------------------------------------------------------------------------
# Main harmonization logic
# ---------------------------------------------------------------------------

def harmonize(input_file, build_info_file, output_file, preferred_build="19",
              n_cores=32, n_override=None, ref_gwas=None, qc_json_file=None):
    """Harmonize a single GWAS summary statistics file.

    Args:
        input_file: path to raw GWAS file
        build_info_file: path to JSON from check_build.py (contains actual_build)
        output_file: path for harmonized output
        preferred_build: target genome build ("19" for GRCh37, "38" for GRCh38)
        n_cores: number of cores for gwaslab operations
        n_override: fixed sample size N (for files missing N column)
        ref_gwas: path to harmonized reference GWAS for allele alignment (optional)
        qc_json_file: path for QC metrics JSON output (optional)
    """
    fname = os.path.basename(input_file)
    print(f"\n=== Harmonizing: {fname} ===", flush=True)
    t0 = time.time()

    # 1. Read build info
    with open(build_info_file) as f:
        build_info = json.load(f)
    actual_build = build_info["actual_build"]
    print(f"  Actual build: {actual_build} (labeled: {build_info['labeled_build']})", flush=True)

    if actual_build == "Unknown":
        raise ValueError(f"Cannot harmonize {fname}: actual build is Unknown")

    # Initialize QC metrics
    qc = {
        "filename": fname,
        "n_variants_raw": 0,
        "build_detected": actual_build,
        "liftover_applied": False,
        "or_to_beta_converted": False,
        "n_override_applied": False,
        "n_dropped_missing": 0,
        "se_derived_from_ci": False,
        "n_alleles_flipped": 0,
        "n_variants_final": 0,
    }

    # 2. Detect columns and load
    columns, sep = read_header(input_file)
    params = detect_params(columns, build=actual_build)
    print(f"  Detected params: {params}", flush=True)
    if sep != "\t":
        print(f"  Detected non-tab separator, passing sep={sep!r} to gwaslab", flush=True)
        params["sep"] = sep

    sumstats = gl.Sumstats(input_file, verbose=False, **params)
    # If no CHR/POS columns, extract from compound SNPID (e.g. "1:15791:C:T")
    split_compound_id(sumstats.data)

    # Compute N from N_CASES + N_CONTROLS if N column is missing
    if "N" not in sumstats.data.columns:
        for cases_col, controls_col in [("N_CASES", "N_CONTROLS"), ("num_cases", "num_controls")]:
            if cases_col in sumstats.data.columns and controls_col in sumstats.data.columns:
                sumstats.data["N"] = (
                    pd.to_numeric(sumstats.data[cases_col], errors="coerce") +
                    pd.to_numeric(sumstats.data[controls_col], errors="coerce")
                )
                print(f"  Computed N = {cases_col} + {controls_col}", flush=True)
                break

    # Drop all-NA columns before basic_check to prevent gwaslab from
    # removing all variants when data is unavailable. Covers both actual
    # NaN and string "NA" values that failed numeric conversion.
    for col_name in ("EAF", "SE"):
        if col_name not in sumstats.data.columns:
            continue
        col = pd.to_numeric(sumstats.data[col_name], errors="coerce")
        if col.isna().all():
            sumstats.data = sumstats.data.drop(columns=[col_name])
            print(f"  Dropped all-NA {col_name} column to prevent basic_check from removing all variants", flush=True)
    # Derive SE from CI bounds if SE was dropped and OR + CI are available
    qc["se_derived_from_ci"] = derive_se_from_ci(
        sumstats.data, input_file, columns, sep
    )

    sumstats.basic_check(verbose=False, n_cores=n_cores)
    qc["n_variants_raw"] = len(sumstats.data)

    # 3. Liftover if actual build differs from preferred build
    if str(actual_build) != str(preferred_build):
        print(f"  Performing liftover GRCh{actual_build} -> GRCh{preferred_build}...", flush=True)
        sumstats.liftover(
            from_build=str(actual_build),
            to_build=str(preferred_build),
            n_cores=n_cores,
        )
        qc["liftover_applied"] = True

    # 4. OR -> BETA conversion
    qc["or_to_beta_converted"] = convert_or_to_beta(sumstats.data)

    # 5. N override (for files missing sample size)
    if n_override is not None:
        if "N" not in sumstats.data.columns or sumstats.data["N"].isna().all():
            sumstats.data["N"] = n_override
            qc["n_override_applied"] = True
            print(f"  Applied N override: {n_override}", flush=True)

    # 6. Post-process: standardize columns, create VAR_ID
    df = postprocess(sumstats)

    # 7. Allele alignment to reference (if provided)
    if ref_gwas is not None:
        print(f"  Aligning alleles to reference: {os.path.basename(ref_gwas)}", flush=True)
        qc["n_alleles_flipped"] = align_alleles_to_reference(df, ref_gwas)

    # 8. Drop rows with missing critical values
    before = len(df)
    df = df.dropna(subset=["VAR_ID", "BETA", "SE", "P_VALUE"])
    qc["n_dropped_missing"] = before - len(df)
    qc["n_variants_final"] = len(df)
    print(f"  Dropped {qc['n_dropped_missing']} rows with missing BETA/SE/P", flush=True)

    # 9. Save harmonized output
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(
        output_file, sep="\t", index=False,
        compression="gzip" if output_file.endswith(".gz") else None,
    )
    print(f"  Saved {len(df)} variants to {output_file}", flush=True)

    # 10. Save QC JSON
    if qc_json_file:
        os.makedirs(os.path.dirname(qc_json_file), exist_ok=True)
        with open(qc_json_file, "w") as f:
            json.dump(qc, f, indent=2)
        print(f"  QC metrics saved to {qc_json_file}", flush=True)

    elapsed = time.time() - t0
    print(f"  Completed in {elapsed:.1f}s", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Harmonize a GWAS summary statistics file"
    )
    parser.add_argument(
        "--input-file", required=True,
        help="Path to raw GWAS summary statistics file"
    )
    parser.add_argument(
        "--build-info", required=True,
        help="Path to build info JSON from check_build.py"
    )
    parser.add_argument(
        "--output-file", required=True,
        help="Path for harmonized output file"
    )
    parser.add_argument(
        "--preferred-build", default="19",
        help="Target genome build: '19' for GRCh37, '38' for GRCh38 (default: 19)"
    )
    parser.add_argument(
        "--n-cores", type=int, default=4,
        help="Number of cores for gwaslab operations (default: 4)"
    )
    parser.add_argument(
        "--n-override", type=int, default=None,
        help="Fixed sample size N for files missing an N column"
    )
    parser.add_argument(
        "--ref-gwas", default=None,
        help="Path to harmonized reference GWAS for allele alignment (e.g., T2D)"
    )
    parser.add_argument(
        "--qc-json", default=None,
        help="Path for QC metrics JSON output"
    )
    args = parser.parse_args()

    harmonize(
        input_file=args.input_file,
        build_info_file=args.build_info,
        output_file=args.output_file,
        preferred_build=args.preferred_build,
        n_cores=args.n_cores,
        n_override=args.n_override,
        ref_gwas=args.ref_gwas,
        qc_json_file=args.qc_json,
    )


if __name__ == "__main__":
    main()
