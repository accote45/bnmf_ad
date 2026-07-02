#!/usr/bin/env bash
# extract_acaf_weights.sh
# Pull the B1 cluster-weight loci out of the All of Us srWGS ACAF-threshold
# plink_bed callset into a small LOCAL plink fileset, ready for scoring by
# compute_prs_aou.sh.
#
# The AoU analysis container has NO FUSE (/dev/fuse absent), so gcsfuse can't
# mount the dataset bucket. Instead we exploit the plink1 .bed layout: a 3-byte
# magic header followed by one FIXED-SIZE block per variant (ceil(N_samples/4)
# bytes), in .bim order. So the byte offset of variant i (0-based) is
#   3 + i * B      where B = ceil(N_samples / 4)
# We stream each fileset's .bim once to find our variants' row indices, then
# fetch ONLY those byte ranges from the .bed with `gsutil -u <proj> cat -r`
# (requester-pays) and reassemble a small local .bed/.bim/.fam. ~885 small range
# reads (~100 KB each), no whole-file download, no FUSE.
#
# Run from the AoU terminal:
#   bash scripts/b1_analysis/aou/extract_acaf_weights.sh
# then score against the extracted set:
#   GENO_DIR=$OUTDIR WEIGHTS=$WEIGHTS bash scripts/b1_analysis/aou/compute_prs_aou.sh
set -euo pipefail

# ============================ PARAMETERS ====================================
ACAF_GS="${ACAF_GS:-gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/plink_bed}"
WEIGHTS="${WEIGHTS:-/home/jupyter/workspaces/duplicateofstatinresponse/multiancestry_polygenic/results/prs/b1_cluster_weights.tsv}"
OUTDIR="${OUTDIR:-/home/jupyter/workspaces/duplicateofstatinresponse/acaf_extracted}"
PROJECT="${BILLING_PROJECT:-${GOOGLE_PROJECT:-}}"     # requester-pays payer (gsutil -u)
PAR="${PAR:-16}"                                      # parallel byte-range fetches
TARGETS="${TARGETS:-$OUTDIR/targets.txt}"            # chr<TAB>pos list from weights

export PROJECT                                        # used by the fetch worker

mkdir -p "$OUTDIR"

echo "=== extract_acaf_weights (FUSE-free / gsutil byte-range): parameters ==="
echo "  ACAF_GS : $ACAF_GS"
echo "  WEIGHTS : $WEIGHTS"
echo "  OUTDIR  : $OUTDIR"
echo "  PROJECT : $PROJECT  (gsutil -u, requester-pays)"
echo "  PARALLEL: $PAR byte-range fetches"
echo "  TARGETS : $TARGETS"
echo

[ -f "$WEIGHTS" ] || { echo "ERROR: weights file not found: $WEIGHTS" >&2; exit 1; }
[ -n "$PROJECT" ] || { echo "ERROR: PROJECT empty (\$GOOGLE_PROJECT unset); requester-pays bucket needs a payer" >&2; exit 1; }
command -v gsutil >/dev/null || { echo "ERROR: gsutil not on PATH" >&2; exit 1; }

# ===== Step 1: build chr<TAB>pos target list from weights ====================
# weights cols: 1 VAR_ID37  2 CHR37  3 POS37  4 CHR38  5 POS38  ...
# Match by GRCh38 position only; alleles are resolved later by compute_prs_aou.sh.
awk -F'\t' 'NR>1 { c=$4; sub(/^chr/,"",c); print "chr"c"\t"$5 }' "$WEIGHTS" \
  | sort -u > "$TARGETS"
echo "Step 1: $(wc -l < "$TARGETS") unique target loci -> $TARGETS"
echo

# ===== Step 2: list the ACAF filesets in the bucket ==========================
echo "Step 2: listing filesets under $ACAF_GS ..."
mapfile -t BEDS < <(gsutil -u "$PROJECT" ls "$ACAF_GS/*.bed" 2>/dev/null | sort)
[ "${#BEDS[@]}" -gt 0 ] || { echo "ERROR: no .bed objects under $ACAF_GS (check access/path)" >&2; exit 1; }
echo "  found ${#BEDS[@]} fileset(s); e.g.:"
printf '    %s\n' "${BEDS[@]:0:3}"
echo

# fetch worker: $1=ord $2=start $3=end ; reads BED_GS, TMPD, PROJECT, BLK from env
# Retries until the block is exactly BLK bytes (guards against short/transient
# gsutil reads that would silently corrupt the reassembled .bed).
fetch_block() {
  local ord="$1" start="$2" end="$3" out="$TMPD/blk_$1.bin" sz=0 t=0
  while [ "$t" -lt 4 ]; do
    gsutil -u "$PROJECT" cat -r "${start}-${end}" "$BED_GS" > "$out" 2>/dev/null || true
    sz=$(stat -c %s "$out" 2>/dev/null || echo 0)
    [ "$sz" -eq "$BLK" ] && return 0
    t=$((t + 1)); sleep 1
  done
  echo "  FETCH-FAIL ord=$ord range=${start}-${end} want=$BLK got=$sz" >&2
  return 1
}
export -f fetch_block

# ===== Step 3: per-fileset byte-range extraction =============================
echo "--- Step 3: extracting loci by byte range ---"
total=0
for bed_gs in "${BEDS[@]}"; do
  prefix_gs="${bed_gs%.bed}"
  base="$(basename "$prefix_gs")"
  out="$OUTDIR/$base"

  # 3a. sample count N from .fam (download it; it's the output .fam too) -> B
  gsutil -u "$PROJECT" cp "${prefix_gs}.fam" "$out.fam" 2>/dev/null \
    || { echo "  $base: no .fam, skipping"; continue; }
  N=$(wc -l < "$out.fam")
  B=$(( (N + 3) / 4 ))                       # bytes per variant block = ceil(N/4)

  # 3b. stream .bim once; keep variants whose chr:pos is a target, with row index
  #     hits columns: <0-based_index>\t<full .bim line>
  gsutil -u "$PROJECT" cat "${prefix_gs}.bim" \
    | awk -v T="$TARGETS" '
        BEGIN { while ((getline l < T) > 0) { split(l, a, "\t"); keep[a[1] "_" a[2]] = 1 } }
        { if (($1 "_" $4) in keep) print (NR-1) "\t" $0 }
      ' > "$out.hits"

  nhits=$(wc -l < "$out.hits")
  if [ "$nhits" -eq 0 ]; then
    echo "  $base: 0 target loci present"
    rm -f "$out.fam" "$out.hits"
    continue
  fi

  # 3c. subset .bim = the .bim lines from hits (cols 2..end), in .bim order
  cut -f2- "$out.hits" > "$out.bim"

  # 3d. compute byte ranges and fetch each block in parallel (ord preserves order)
  TMPD="$out.blocks"; rm -rf "$TMPD"; mkdir -p "$TMPD"
  export BED_GS="$bed_gs" TMPD BLK="$B"
  # %.0f (not %d): byte offsets exceed 2^31, and mawk's %d saturates at INT_MAX.
  awk -v B="$B" '{ i=$1; s=3+i*B; e=s+B-1; printf "%06d %.0f %.0f\n", NR, s, e }' "$out.hits" \
    | xargs -P "$PAR" -n 3 bash -c 'fetch_block "$@"' _

  # 3e. assemble .bed = real 3-byte magic header from source + blocks in order
  gsutil -u "$PROJECT" cat -r 0-2 "$bed_gs" > "$out.bed"      # 0x6c 0x1b 0x01 (SNP-major)
  cat "$TMPD"/blk_*.bin >> "$out.bed"                          # glob = sorted = .bim order

  # verify: block count AND exact .bed size = 3 + nhits*B (plink1 .bed invariant)
  nblk=$(ls "$TMPD"/blk_*.bin 2>/dev/null | wc -l)
  [ "$nblk" -eq "$nhits" ] || { echo "  $base: ERROR fetched $nblk/$nhits blocks" >&2; exit 1; }
  expected=$(( 3 + nhits * B )); actual=$(stat -c %s "$out.bed")
  [ "$actual" -eq "$expected" ] || { echo "  $base: ERROR .bed size $actual != expected $expected" >&2; exit 1; }
  rm -rf "$TMPD" "$out.hits"

  total=$((total + nhits))
  echo "  $base: $nhits variants (N=$N, block=$B bytes)"
done

echo
echo "=== Done. Extracted $total variants total into: $OUTDIR ==="
echo "Score against them with:"
echo "  GENO_DIR=$OUTDIR WEIGHTS=$WEIGHTS bash scripts/b1_analysis/aou/compute_prs_aou.sh"
