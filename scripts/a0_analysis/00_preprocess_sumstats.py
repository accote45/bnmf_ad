#!/usr/bin/env python3
"""
00_preprocess_sumstats.py — Preprocess N harmonized EUR sumstats for LDSC.

Handles an arbitrary number of traits read from a YAML config.

Steps per trait:
  1. Load harmonized sumstats
  2. Drop rows with missing core fields
  3. Exclude multiallelic SNPs (deduplicate on VAR_ID)
  4. Filter MAF > threshold
  5. Filter to HapMap3 SNPs (inner join on VAR_ID)
  6. Derive A2 from VAR_ID
  7. Remove strand-ambiguous SNPs
  8. Resolve RSIDs (prefer HapMap3 rsid over ".")

Cross-trait steps:
  9. Align each non-reference trait's effect allele to the reference
 10. Restrict all traits to the shared SNP universe (intersection)
 11. Output LDSC-ready files: SNP, A1, A2, BETA, SE, P, N

Usage:
  python scripts/a0_analysis/00_preprocess_sumstats.py --config config/a0_5trait_config.yaml
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

STRAND_AMBIGUOUS = {frozenset({"A", "T"}), frozenset({"C", "G"})}


def is_strand_ambiguous(a1, a2):
    return frozenset({a1, a2}) in STRAND_AMBIGUOUS


def extract_a2(var_id, effect_allele):
    """Return the non-effect allele from VAR_ID (CHR_POS_A1_A2 sorted)."""
    parts = var_id.split("_")
    if len(parts) != 4:
        return None
    alleles = set([parts[2], parts[3]])
    other = alleles - {effect_allele.upper()}
    if len(other) == 1:
        return other.pop()
    return None


def log_step(label, n, log_lines):
    line = f"{label}: {n:,} SNPs"
    print(line)
    log_lines.append(line)


def filter_and_prep(df, label, hm3, maf_threshold, log_lines, exclude_mhc=True):
    """Filter and prepare a single trait DataFrame for LDSC."""
    log_lines.append(f"\n--- {label} ---")
    log_step(f"{label} raw", len(df), log_lines)

    # Drop rows with missing core fields
    df = df.dropna(subset=["BETA", "SE", "P_VALUE", "MAF", "Effect_Allele"])
    log_step(f"{label} after drop NA", len(df), log_lines)

    # Exclude multiallelic: keep VAR_IDs that appear exactly once
    counts = df["VAR_ID"].value_counts()
    df = df[df["VAR_ID"].isin(counts[counts == 1].index)]
    log_step(f"{label} after multiallelic exclusion", len(df), log_lines)

    # MAF filter
    df = df[df["MAF"] > maf_threshold]
    log_step(f"{label} after MAF > {maf_threshold}", len(df), log_lines)

    # Exclude MHC region (chr6:25,000,000-35,000,000 GRCh37)
    if exclude_mhc:
        parts = df["VAR_ID"].str.split("_", n=2)
        chr_col = parts.str[0].astype(str)
        pos_col = parts.str[1].astype(int)
        in_mhc = (chr_col == "6") & (pos_col >= 25_000_000) & (pos_col <= 35_000_000)
        df = df[~in_mhc]
        log_step(f"{label} after MHC exclusion (chr6:25-35Mb)", len(df), log_lines)

    # Filter to HapMap3 SNPs
    df = df.merge(hm3, on="VAR_ID", how="inner")
    log_step(f"{label} after HapMap3 filter", len(df), log_lines)

    # Uppercase alleles for consistency
    df["Effect_Allele"] = df["Effect_Allele"].str.upper()

    # Derive A2
    df["A2"] = df.apply(
        lambda r: extract_a2(r["VAR_ID"], r["Effect_Allele"]), axis=1
    )
    df = df.dropna(subset=["A2"])
    log_step(f"{label} after A2 derivation", len(df), log_lines)

    # Drop strand-ambiguous SNPs
    df = df[~df.apply(lambda r: is_strand_ambiguous(r["Effect_Allele"], r["A2"]), axis=1)]
    log_step(f"{label} after removing strand-ambiguous", len(df), log_lines)

    # Prefer HapMap3 rsid; fall back to file RSID only if it looks like an rsID
    has_valid_rsid = (
        df["RSID"].notna()
        & (df["RSID"] != ".")
        & df["RSID"].str.startswith("rs", na=False)
    )
    df["SNP"] = np.where(has_valid_rsid, df["RSID"], df["hm3_rsid"])
    df = df.dropna(subset=["SNP"])
    df = df[df["SNP"] != "."]
    log_step(f"{label} after rsid resolution", len(df), log_lines)

    return df


def align_to_reference(query_df, ref_df, query_label, log_lines):
    """Align a query trait's effect alleles to the reference trait.

    Returns the aligned query DataFrame restricted to shared SNPs.
    """
    log_lines.append(f"\n--- Allele alignment ({query_label} → reference) ---")

    # Get reference allele info
    ref_alleles = ref_df[["VAR_ID", "Effect_Allele", "A2"]].rename(
        columns={"Effect_Allele": "REF_Effect", "A2": "REF_A2"}
    )

    # Merge query with reference alleles
    merged = query_df.merge(ref_alleles, on="VAR_ID", how="inner")
    log_step(f"{query_label} shared SNPs with reference", len(merged), log_lines)

    # Concordant: same effect allele
    concordant = merged["Effect_Allele"] == merged["REF_Effect"]

    # Discordant but flippable: query effect allele == reference A2
    flippable = merged["Effect_Allele"] == merged["REF_A2"]

    # Flip where discordant
    merged.loc[flippable, "BETA"] = -merged.loc[flippable, "BETA"]
    merged.loc[flippable, "Effect_Allele"] = merged.loc[flippable, "REF_Effect"]
    merged.loc[flippable, "A2"] = merged.loc[flippable, "REF_A2"]

    n_concordant = concordant.sum()
    n_flipped = flippable.sum()
    n_dropped = len(merged) - n_concordant - n_flipped

    log_lines.append(f"  Concordant (no flip needed): {n_concordant:,}")
    log_lines.append(f"  Flipped (BETA sign reversed): {n_flipped:,}")
    log_lines.append(f"  Dropped (unresolvable allele mismatch): {n_dropped:,}")
    print(f"  {query_label}: Concordant={n_concordant:,} | Flipped={n_flipped:,} | Dropped={n_dropped:,}")

    # Keep only concordant + flipped
    merged = merged[concordant | flippable].copy()
    merged = merged.drop(columns=["REF_Effect", "REF_A2"])
    log_step(f"{query_label} after alignment", len(merged), log_lines)

    return merged


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Preprocess N harmonized sumstats for LDSC"
    )
    parser.add_argument("--config", default="config/a0_config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    results_dir   = cfg["results_dir"]
    hapmap3_path  = cfg["reference"]["hapmap3"]
    maf_threshold = float(cfg["filters"]["maf_threshold"])
    exclude_mhc   = cfg["filters"].get("exclude_mhc", True)
    traits_cfg    = cfg["traits"]

    os.makedirs(results_dir, exist_ok=True)
    log_lines = []

    # Identify reference trait
    ref_key = None
    for name, info in traits_cfg.items():
        if info.get("is_reference", False):
            ref_key = name
            break
    if ref_key is None:
        print("ERROR: No trait marked with is_reference: true in config")
        sys.exit(1)
    print(f"Reference trait: {ref_key} ({traits_cfg[ref_key]['label']})")
    log_lines.append(f"Reference trait: {ref_key} ({traits_cfg[ref_key]['label']})")

    # -----------------------------------------------------------------------
    # Load HapMap3 SNP list
    # -----------------------------------------------------------------------
    print("\nLoading HapMap3 SNP list...")
    hm3 = pd.read_csv(hapmap3_path, sep="\t", usecols=["VAR_ID", "rsid"])
    hm3 = hm3.rename(columns={"rsid": "hm3_rsid"})
    hm3 = hm3.drop_duplicates("VAR_ID")
    print(f"  HapMap3 SNPs loaded: {len(hm3):,}")

    # -----------------------------------------------------------------------
    # Load and filter each trait
    # -----------------------------------------------------------------------
    cols = ["VAR_ID", "RSID", "Effect_Allele", "P_VALUE", "BETA", "SE", "N", "MAF"]
    trait_dfs = {}

    for name, info in traits_cfg.items():
        path = info["file"]
        label = info.get("label", name.upper())
        print(f"\nLoading {label}: {path}")
        df = pd.read_csv(path, sep="\t", usecols=cols, compression="gzip")
        # Optional per-trait total-N override: replace the harmonized N column
        # (which may be effective N, e.g. Suzuki T2D) with a fixed total sample
        # size so LDSC h2 is estimated on the observed 0/1 scale. Only applied
        # when `n_total` is set in the config for this trait.
        n_total = info.get("n_total")
        if n_total is not None:
            df["N"] = float(n_total)
            line = f"{label}: N overridden to total = {int(n_total):,}"
            print(f"  {line}")
            log_lines.append(line)
        trait_dfs[name] = filter_and_prep(df, label, hm3, maf_threshold, log_lines, exclude_mhc=exclude_mhc)

    # -----------------------------------------------------------------------
    # Align all non-reference traits to reference
    # -----------------------------------------------------------------------
    ref_df = trait_dfs[ref_key]

    for name in list(trait_dfs.keys()):
        if name != ref_key:
            label = traits_cfg[name].get("label", name.upper())
            trait_dfs[name] = align_to_reference(
                trait_dfs[name], ref_df, label, log_lines
            )

    # -----------------------------------------------------------------------
    # Restrict all traits to shared SNP universe (intersection)
    # -----------------------------------------------------------------------
    log_lines.append("\n--- Shared SNP universe ---")
    shared_ids = None
    for name, df in trait_dfs.items():
        ids = set(df["VAR_ID"])
        if shared_ids is None:
            shared_ids = ids
        else:
            shared_ids = shared_ids & ids
    log_step("Shared SNP universe (intersection of all traits)", len(shared_ids), log_lines)

    for name in trait_dfs:
        trait_dfs[name] = trait_dfs[name][
            trait_dfs[name]["VAR_ID"].isin(shared_ids)
        ].copy()
        label = traits_cfg[name].get("label", name.upper())
        log_step(f"{label} final", len(trait_dfs[name]), log_lines)

    # -----------------------------------------------------------------------
    # Output
    # -----------------------------------------------------------------------
    out_cols = ["SNP", "Effect_Allele", "A2", "BETA", "SE", "P_VALUE", "N"]
    rename_map = {"Effect_Allele": "A1", "P_VALUE": "P"}

    log_lines.append("\n--- Output files ---")
    for name, df in trait_dfs.items():
        out_file = os.path.join(results_dir, f"{name}_for_munge.txt")
        df[out_cols].rename(columns=rename_map).to_csv(
            out_file, sep="\t", index=False
        )
        log_lines.append(f"  {name}: {out_file}  ({len(df):,} SNPs)")
        print(f"Output: {out_file}  ({len(df):,} SNPs)")

    summary_path = os.path.join(results_dir, "preprocess_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")

    print(f"\nSummary written to: {summary_path}")


if __name__ == "__main__":
    main()
