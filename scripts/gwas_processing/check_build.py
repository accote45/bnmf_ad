#!/usr/bin/env python3
"""
Check the actual genome build of GWAS files using gwaslab's infer_build.
Loads each file with auto-detected column mappings and uses HapMap3 reference
SNPs to determine whether coordinates are GRCh37 (hg19) or GRCh38 (hg38).

Usage:
    # Multi-file mode (standalone)
    python check_build.py                        # check all available ancestries
    python check_build.py --ancestry AFR         # check AFR files only
    python check_build.py --ancestry EUR AFR     # check EUR and AFR

    # Single-file mode (for Snakemake pipeline)
    python check_build.py --input-file sumstats/file.txt.gz --output-file build.json
"""

import argparse
import glob
import json
import os

import gwaslab as gl
import numpy as np
import pandas as pd

from column_map import COLUMN_MAP, read_header, detect_params, split_compound_id

SUMSTATS_DIR = "/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/sumstats"
ANCESTRIES = ["EUR", "AFR", "HIS", "SAS", "EAS", "META"]


def discover_files(sumstats_dir, ancestry):
    """Find GWAS data files for a given ancestry in the sumstats directory."""
    pattern = os.path.join(sumstats_dir, f"*.{ancestry}.*")
    data_exts = (".txt", ".txt.gz", ".tsv", ".tsv.gz")
    files = [f for f in glob.glob(pattern) if any(f.endswith(ext) for ext in data_exts)]
    return sorted(files)


def parse_trait(filename):
    """Extract trait name from filename: Author_Journal_Year.TRAIT.Ancestry.GRCh##.ext"""
    parts = os.path.basename(filename).split(".")
    return parts[1] if len(parts) >= 3 else os.path.basename(filename)


def get_labeled_build(filename):
    """Extract the labeled genome build from the filename."""
    return "GRCh38" if "GRCh38" in filename else "GRCh37"


def infer_build_for_file(filepath):
    """Load a GWAS file with auto-detected columns and infer its genome build.

    Returns:
        str: inferred build ("19", "38", or "Unknown")
    """
    columns, sep = read_header(filepath)
    params = detect_params(columns)
    print(f"  Detected params: {params}", flush=True)
    if sep != "\t":
        params["sep"] = sep

    sumstats = gl.Sumstats(filepath, verbose=False, **params)

    # If no CHR/POS columns, extract from compound SNPID (e.g. "1:15791:C:T" or "1_10616_C_G")
    split_compound_id(sumstats.data)

    sumstats.infer_build(verbose=True)
    return sumstats.meta.get("gwaslab", {}).get("genome_build", "Unknown")


def check_single_file(input_file, output_file):
    """Check build of a single file and write JSON result (for Snakemake)."""
    fname = os.path.basename(input_file)
    print(f"\n--- Checking build: {fname} ---", flush=True)

    try:
        build = infer_build_for_file(input_file)
        result = {
            "filename": fname,
            "labeled_build": get_labeled_build(fname),
            "actual_build": str(build),
        }
        print(f"  >>> INFERRED BUILD: {build}", flush=True)
    except Exception as e:
        print(f"  ERROR: {e}", flush=True)
        result = {
            "filename": fname,
            "labeled_build": get_labeled_build(fname),
            "actual_build": "Unknown",
            "error": str(e),
        }

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Wrote build info to {output_file}", flush=True)


def check_multi_file(ancestries, sumstats_dir):
    """Check builds for all files matching the given ancestries (standalone mode)."""
    all_results = {}

    for ancestry in ancestries:
        files = discover_files(sumstats_dir, ancestry)
        if not files:
            continue

        print(f"\n{'='*70}")
        print(f"GENOME BUILD VERIFICATION FOR {ancestry} GWAS FILES")
        print(f"{'='*70}")

        for filepath in files:
            fname = os.path.basename(filepath)
            trait = parse_trait(fname)
            print(f"\n--- {trait} ({ancestry}): {fname} ---", flush=True)

            try:
                build = infer_build_for_file(filepath)
                all_results[(ancestry, trait)] = (fname, build)
                print(f"  >>> INFERRED BUILD: {build}", flush=True)
            except Exception as e:
                print(f"  ERROR: {e}", flush=True)
                all_results[(ancestry, trait)] = (fname, f"ERROR: {e}")

    # Summary
    if all_results:
        print(f"\n\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        print(f"{'Ancestry':<10} {'Trait':<10} {'Filename':<60} {'Label':>8} {'Actual':>8} {'Match':>8}")
        print("-" * 110)
        for (ancestry, trait), (fname, build) in sorted(all_results.items()):
            label = get_labeled_build(fname)
            actual_str = f"GRCh{build}" if build in ("19", "38") else str(build)
            match = (
                "OK" if (label == "GRCh37" and build == "19") or (label == "GRCh38" and build == "38")
                else "MISMATCH"
            )
            print(f"{ancestry:<10} {trait:<10} {fname:<60} {label:>8} {actual_str:>8} {match:>8}")
    else:
        print("\nNo GWAS files found.")


def main():
    parser = argparse.ArgumentParser(
        description="Check genome build of GWAS summary statistics files"
    )

    # Single-file mode (for Snakemake)
    parser.add_argument(
        "--input-file", default=None,
        help="Single GWAS file to check (enables single-file mode)"
    )
    parser.add_argument(
        "--output-file", default=None,
        help="Output JSON file for build result (requires --input-file)"
    )

    # Multi-file mode (standalone)
    parser.add_argument(
        "--ancestry", nargs="+", choices=ANCESTRIES, default=None,
        help="Ancestry/ancestries to check (default: all available)"
    )
    parser.add_argument(
        "--sumstats-dir", default=SUMSTATS_DIR,
        help="Directory containing GWAS summary statistics files"
    )
    args = parser.parse_args()

    if args.input_file:
        if not args.output_file:
            parser.error("--output-file is required when using --input-file")
        check_single_file(args.input_file, args.output_file)
    else:
        ancestries = args.ancestry if args.ancestry else ANCESTRIES
        check_multi_file(ancestries, args.sumstats_dir)


if __name__ == "__main__":
    main()
