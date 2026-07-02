#!/bin/bash
set -euo pipefail

SUMSTATS="/sc/arion/projects/paul_oreilly/lab/lioul01/multiancestry_polygenic/sumstats"
mkdir -p "${SUMSTATS}"

# Array of "URL|TARGET_FILENAME" pairs
FILES=(
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90693001-GCST90694000/GCST90693007/harmonised/GCST90693007.h.tsv.gz|Karczewski_NatGenet_2025.atherosclerosis.META.GRCh37.txt.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90691001-GCST90692000/GCST90691962/harmonised/GCST90691962.h.tsv.gz|Karczewski_NatGenet_2025.atherosclerosis.EUR.GRCh37.txt.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90692001-GCST90693000/GCST90692424/harmonised/GCST90692424.h.tsv.gz|Karczewski_NatGenet_2025.atherosclerosis.SAS.GRCh37.txt.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480132/harmonised/GCST90480132.h.tsv.gz|Verma_Science_2024.atherosclerosis.META.GRCh37.tsv.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475936/harmonised/GCST90475936.h.tsv.gz|Verma_Science_2024.atherosclerosis.EUR.GRCh37.tsv.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475935/harmonised/GCST90475935.h.tsv.gz|Verma_Science_2024.atherosclerosis.AFR.GRCh37.tsv.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90477001-GCST90478000/GCST90477868/harmonised/GCST90477868.h.tsv.gz|Verma_Science_2024.atherosclerosis.EAS.GRCh37.tsv.gz"
    "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475934/harmonised/GCST90475934.h.tsv.gz|Verma_Science_2024.atherosclerosis.AMR.GRCh37.tsv.gz"
)

SUCCESS=0
SKIPPED=0
FAILED=0

for ENTRY in "${FILES[@]}"; do
    IFS='|' read -r URL FNAME <<< "${ENTRY}"
    OUTPATH="${SUMSTATS}/${FNAME}"
    if [ -f "${OUTPATH}" ]; then
        echo "SKIP (exists): ${FNAME}"
        ((SKIPPED++))
        continue
    fi
    echo "Downloading: ${FNAME}"
    if wget -q -O "${OUTPATH}" "${URL}"; then
        echo "  OK: ${FNAME}"
        ((SUCCESS++))
    else
        echo "  FAILED: ${FNAME}"
        rm -f "${OUTPATH}"
        ((FAILED++))
    fi
done

echo ""
echo "Done. Success: ${SUCCESS}, Skipped: ${SKIPPED}, Failed: ${FAILED}"
