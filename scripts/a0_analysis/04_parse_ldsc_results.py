#!/usr/bin/env python3
"""
04_parse_ldsc_results.py — Parse LDSC --rg log files into a tidy CSV.

Reads each rg_{comparison}.log file, extracts the summary table, and
produces a single tidy CSV with one row per pairwise comparison.

Usage:
  python scripts/a0_analysis/04_parse_ldsc_results.py \
      --config config/a0_config.yaml \
      --results-dir results/a0_analysis \
      --output results/a0_analysis/rg_results.csv
"""

import argparse
import os
import re

import pandas as pd
import yaml


def parse_ldsc_log(log_path):
    """Parse an LDSC --rg log file and return a dict of results."""
    with open(log_path) as f:
        text = f.read()

    results = {}

    # Parse heritability of phenotype 1 (reference trait).
    # The trailing "Ratio:" line varies ("Ratio: x (y)" vs
    # "Ratio < 0 (usually indicates GC correction).") so it is not matched.
    m = re.search(r"Heritability of phenotype 1\n-+\n"
                  r"Total Observed scale h2: ([\d.e+-]+) \(([\d.e+-]+)\)\n"
                  r"Lambda GC: ([\d.e+-]+)\n"
                  r"Mean Chi\^2: ([\d.e+-]+)\n"
                  r"Intercept: ([\d.e+-]+) \(([\d.e+-]+)\)", text)
    if m:
        results["h2_obs_trait1"] = float(m.group(1))
        results["h2_obs_se_trait1"] = float(m.group(2))
        results["h2_int_trait1"] = float(m.group(5))
        results["h2_int_se_trait1"] = float(m.group(6))

    # Parse heritability of phenotype 2 (comparison trait)
    m = re.search(r"Heritability of phenotype 2/2\n-+\n"
                  r"Total Observed scale h2: ([\d.e+-]+) \(([\d.e+-]+)\)\n"
                  r"Lambda GC: ([\d.e+-]+)\n"
                  r"Mean Chi\^2: ([\d.e+-]+)\n"
                  r"Intercept: ([\d.e+-]+) \(([\d.e+-]+)\)", text)
    if m:
        results["h2_obs_trait2"] = float(m.group(1))
        results["h2_obs_se_trait2"] = float(m.group(2))
        results["h2_int_trait2"] = float(m.group(5))
        results["h2_int_se_trait2"] = float(m.group(6))

    # Parse genetic correlation
    m = re.search(r"Genetic Correlation: ([\d.e+-]+) \(([\d.e+-]+)\)\n"
                  r"Z-score: ([\d.e+-]+)\n"
                  r"P: ([\d.e+-]+)", text)
    if m:
        results["rg"] = float(m.group(1))
        results["se"] = float(m.group(2))
        results["z"] = float(m.group(3))
        results["p"] = float(m.group(4))

    # Parse genetic covariance intercept
    m = re.search(r"Genetic Covariance\n-+\n"
                  r"Total Observed scale gencov: ([\d.e+-]+) \(([\d.e+-]+)\)\n"
                  r"Mean z1\*z2: ([\d.e+-]+)\n"
                  r"Intercept: ([\d.e+-]+) \(([\d.e+-]+)\)", text)
    if m:
        results["gcov_int"] = float(m.group(4))
        results["gcov_int_se"] = float(m.group(5))

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Parse LDSC --rg logs into tidy CSV"
    )
    parser.add_argument("--config", default="config/a0_config.yaml")
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    traits_cfg = cfg["traits"]

    # Identify reference trait
    ref_key = None
    for name, info in traits_cfg.items():
        if info.get("is_reference", False):
            ref_key = name
            break

    ref_label = traits_cfg[ref_key]["label"]
    comparisons = [t for t in traits_cfg if t != ref_key]

    rows = []
    for comp in comparisons:
        log_path = os.path.join(args.results_dir, f"rg_{comp}.log")
        if not os.path.exists(log_path):
            print(f"WARNING: {log_path} not found, skipping")
            continue

        comp_label = traits_cfg[comp]["label"]
        results = parse_ldsc_log(log_path)
        results["trait1"] = ref_label
        results["trait2"] = comp_label
        rows.append(results)

    # Define tidy column order
    col_order = [
        "trait1", "trait2",
        "rg", "se", "z", "p",
        "h2_obs_trait1", "h2_obs_se_trait1",
        "h2_obs_trait2", "h2_obs_se_trait2",
        "h2_int_trait1", "h2_int_se_trait1",
        "h2_int_trait2", "h2_int_se_trait2",
        "gcov_int", "gcov_int_se",
    ]

    df = pd.DataFrame(rows)
    # Reorder columns (only include those present)
    df = df[[c for c in col_order if c in df.columns]]
    df.to_csv(args.output, index=False)
    print(f"Tidy results written to: {args.output}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
