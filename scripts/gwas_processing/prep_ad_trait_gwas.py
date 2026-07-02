#!/usr/bin/env python3
"""prep_ad_trait_gwas.py

Normalize AD trait GWAS that the standard harmonizer can't ingest as-is, so the
output can go straight into check_build.py + harmonize_sumstats.py.

Modes:
  --kind neuroticism : Nagel 2018 CTG format has Z but no BETA/SE. Derive
                       BETA/SE from Z, MAF_UKB, N via Zhu et al. 2016
                       (BETA/SE == Z, lossless for the bNMF Z-matrix).
  --kind pgc_vcf     : PGC sumstats-VCF (MDD 2025, SCZ wave3) has '##' comment
                       lines and a '#'-prefixed header. Strip them.

Usage:
  python prep_ad_trait_gwas.py --input <raw> --output <clean> --kind <mode>
"""
import argparse
import gzip

import numpy as np
import pandas as pd


def open_any(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def prep_neuroticism(inp, out):
    df = pd.read_csv(inp, sep="\t")
    for col in ("Z", "MAF_UKB", "N"):
        if col not in df.columns:
            raise SystemExit(f"neuroticism: expected column '{col}' not found in {inp}")
    z = pd.to_numeric(df["Z"], errors="coerce")
    maf = pd.to_numeric(df["MAF_UKB"], errors="coerce")
    n = pd.to_numeric(df["N"], errors="coerce")
    # Zhu et al. (Nat Genet 2016): with these, BETA/SE reproduces Z exactly.
    denom = 2.0 * maf * (1.0 - maf) * (n + z ** 2)
    denom = denom.where(denom > 0)          # guard div-by-zero / MAF=0
    root = np.sqrt(denom)
    df["BETA"] = z / root
    df["SE"] = 1.0 / root
    n_bad = int(df["BETA"].isna().sum())
    df.to_csv(out, sep="\t", index=False,
              compression="gzip" if out.endswith(".gz") else None)
    print(f"neuroticism: wrote {len(df)} rows ({n_bad} with unusable MAF/Z) -> {out}",
          flush=True)


def prep_pgc_vcf(inp, out):
    n_meta = 0
    n_data = 0
    with open_any(inp) as fin, open_any(out, "wt") as fout:
        header_written = False
        for line in fin:
            if line.startswith("##"):
                n_meta += 1
                continue
            if not header_written:
                if line.startswith("#"):      # '#CHROM...' -> 'CHROM...'
                    line = line[1:]
                fout.write(line)
                header_written = True
            else:
                fout.write(line)
                n_data += 1
    print(f"pgc_vcf: skipped {n_meta} '##' lines, wrote header + {n_data} rows -> {out}",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--kind", required=True, choices=["neuroticism", "pgc_vcf"])
    a = ap.parse_args()
    if a.kind == "neuroticism":
        prep_neuroticism(a.input, a.output)
    else:
        prep_pgc_vcf(a.input, a.output)


if __name__ == "__main__":
    main()
