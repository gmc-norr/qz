#!/usr/bin/env bash
# BAM benchmark: bz (self-contained, reference-free) vs samtools CRAM 3.0 / 3.1.
#
# bz roundtrip is lossless at the RECORD level but re-frames BGZF, so the
# decompressed BAM is NOT byte-identical to the input; we verify losslessness by
# diffing `samtools view` (records) of input vs roundtrip. CRAM is the ratio/speed
# competitor — it needs a reference whose contigs MATCH the BAM header (see the
# REF_FULL/REF_CHR20 note below); bench_cram fails the row if samtools falls back to
# embed_ref. CRAM's roundtrip is not re-verified here.
#
# bz requires a coordinate-sorted BAM (@HD SO:coordinate). If the input is labeled
# SO:unsorted but is actually sorted (has a .bai), we reheader it into a cached copy.
#
# Default runs only the chr20 tier; set DO_FULL=1 for the deferred 51 GB WGS tier.
# Override via env: BAM_FULL, BAM_CHR20, REF_FULL, REF_CHR20, DO_FULL, BZ, SAM, T,
#                   WORK, RES, LOG, MIN_FREE_GB.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

BAM_FULL=${BAM_FULL:-$HERE/../real_data/HG002.novaseq.pcr-free.35x.bwamem2.dedup.grch38_no_alt.bam}
BAM_CHR20=${BAM_CHR20:-$HERE/../real_data/HG002.GRCh38.chr20.bam}
# CRAM needs a reference whose contigs match the BAM header — INCLUDING @SQ contigs
# no read aligns to, because samtools resolves every @SQ at header-write time. A
# single miss makes it silently enable embed_ref=2 (it embeds the reference INTO the
# CRAM), which bloats the file and yields a bogus ratio. So the reference is per-BAM
# (bench_cram aborts the row if embed_ref is triggered anyway):
#   * the WGS BAM is aligned to grch38_no_alt (no decoys)            → REF_FULL.
#   * the chr20 BAM is aligned to GRCh38_full_plus_hs38d1 (2580 SQ,
#     incl. decoys)                                                  → REF_CHR20 (hs38DH).
REF_FULL=${REF_FULL:-${REF:-$HERE/../real_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}}
REF_CHR20=${REF_CHR20:-$HERE/../real_data/reference/Homo_sapiens_assembly38.fasta}
DO_FULL=${DO_FULL:-0}             # set DO_FULL=1 to also run the deferred 51 GB WGS tier
BZ=${BZ:-$HERE/../target/release/bz}
SAM=${SAM:-samtools}
WORK=${WORK:-/tmp/qzbench_bam}
RES=${RES:-$HERE/results/bam.tsv}
LOG=${LOG:-$HERE/results/bam.log}
MIN_FREE_GB=${MIN_FREE_GB:-200}

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"

free_gb=$(df -BG --output=avail "$WORK" | tail -1 | tr -dc '0-9')
if [ "${free_gb:-0}" -lt "$MIN_FREE_GB" ]; then
  echo "WARN: only ${free_gb}G free at $WORK (<$MIN_FREE_GB G); the 51GB BAM rows may fail on disk." | tee -a "$LOG" >&2
fi

if ! need_tool "$SAM" samtools; then echo "samtools required for BAM bench" >&2; exit 1; fi

# fingerprints (computed once per BAM): full record stream + QUAL-stripped
prep_fingerprints() { # bam ref -> sets FP_FULL, FP_NOQUAL, CUR_BYTES, CUR_BAM, CUR_REF
  CUR_BAM="$1"; CUR_REF="$2"
  [ -f "$CUR_REF.fai" ] || "$SAM" faidx "$CUR_REF" 2>>"$LOG" || echo "WARN: faidx $CUR_REF failed" >&2
  # bz requires @HD SO:coordinate; reheader header-only if mislabeled. The chr20
  # BAM is genuinely coordinate-sorted at the record level (verified: 0 out-of-order
  # POS over 18.3M records) but carries a stale @HD SO:unsorted label, so this
  # header-only rewrite runs for it. The records are NOT touched/re-sorted.
  local so; so=$("$SAM" view -H "$CUR_BAM" 2>/dev/null | grep -m1 -oP 'SO:\K[a-z]+' || echo "")
  if [ "$so" != "coordinate" ]; then
    echo "input @HD SO:$so != coordinate; reheadering to SO:coordinate (header-only rewrite)" | tee -a "$LOG"
    local pre="$WORK/$(basename "$CUR_BAM").sorted.bam"
    "$SAM" view -H "$CUR_BAM" 2>/dev/null | sed 's/SO:unsorted/SO:coordinate/' > "$WORK/hdr.sam"
    "$SAM" reheader "$WORK/hdr.sam" "$CUR_BAM" > "$pre" 2>>"$LOG"
    CUR_BAM="$pre"
  fi
  CUR_BYTES=$(stat -c %s "$CUR_BAM")
  FP_FULL=$("$SAM" view -@ "$T" "$CUR_BAM" 2>>"$LOG" | md5sum | awk '{print $1}')
  FP_NOQUAL=$("$SAM" view -@ "$T" "$CUR_BAM" 2>>"$LOG" | awk 'BEGIN{OFS="\t"}{$11="*";print}' | md5sum | awk '{print $1}')
}

bench_bz() { # label lossless(yes|no) <extra bz flags...>
  local label="$1" lossless="$2"; shift 2
  local arc="$WORK/$label.bz" rt="$WORK/$label.rt.bam"
  section "bz $label"
  run_timed "$LOG" -- "$BZ" compress -i "$CUR_BAM" -o "$arc" -t "$T" -f --numa off "$@"
  local crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  if [ "$crc" -ne 0 ] || [ "$cbytes" -eq 0 ]; then
    emit_row "$RES" "$label" "$lossless" "0" "$CUR_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"; rm -f "$arc"; return; fi
  run_timed "$LOG" -- "$BZ" decompress -i "$arc" -o "$rt" -t "$T" -f
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU rtchk m
  if [ "$lossless" = yes ]; then
    m=$("$SAM" view -@ "$T" "$rt" 2>>"$LOG" | md5sum | awk '{print $1}')
    [ "$m" = "$FP_FULL" ] && rtchk=OK || rtchk=MISMATCH
  else
    m=$("$SAM" view -@ "$T" "$rt" 2>>"$LOG" | awk 'BEGIN{OFS="\t"}{$11="*";print}' | md5sum | awk '{print $1}')
    [ "$m" = "$FP_NOQUAL" ] && rtchk=QUAL_ONLY_DIFF || rtchk=NONQUAL_MISMATCH
  fi
  emit_row "$RES" "$label" "$lossless" "$cbytes" "$CUR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rtchk"
  rm -f "$arc" "$rt"
}

bench_cram() { # label version(3.0|3.1) [extra --output-fmt-option args...]
  local label="$1" ver="$2"; shift 2
  local cram="$WORK/$label.cram" rt="$WORK/$label.rt.bam"
  section "samtools CRAM $ver ($label)  ref=$(basename "$CUR_REF")"
  local log0; log0=$(wc -l < "$LOG")
  run_timed "$LOG" -- "$SAM" view -C -T "$CUR_REF" --output-fmt-option version=$ver "$@" -@ "$T" -o "$cram" "$CUR_BAM"
  local crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$cram" 2>/dev/null || echo 0)
  # Reference/contig mismatch → samtools silently embeds the reference (embed_ref=2)
  # and the CRAM bloats; that ratio is meaningless, so fail the row loudly instead
  # of publishing it (this is the bug that made CRAM look like 1.29× — it was 1.83×).
  if tail -n +$((log0+1)) "$LOG" | grep -q "Enabling embed_ref"; then
    echo "FATAL: CRAM $label fell back to embed_ref — $(basename "$CUR_REF") does not match $(basename "$CUR_BAM") contigs; ratio invalid" | tee -a "$LOG" >&2
    emit_row "$RES" "$label" "yes" "0" "$CUR_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "EMBED_REF_FAIL"; rm -f "$cram"; return; fi
  if [ "$crc" -ne 0 ] || [ "$cbytes" -eq 0 ]; then
    emit_row "$RES" "$label" "yes" "0" "$CUR_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"; rm -f "$cram"; return; fi
  run_timed "$LOG" -- "$SAM" view -b -T "$CUR_REF" -@ "$T" -o "$rt" "$cram"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  emit_row "$RES" "$label" "yes" "$cbytes" "$CUR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "n-a"
  rm -f "$cram" "$rt"
}

bench_pigz() { # label  — generic re-gzip baseline of the (already-BGZF) BAM
  command -v pigz >/dev/null 2>&1 || { echo "pigz not found; skipping $1" | tee -a "$LOG"; return; }
  local label="$1"
  local gz="$WORK/$label.gz" rt="$WORK/$label.rt"
  section "pigz -9 ($label)"
  run_timed "$LOG" -- bash -c "pigz -9 -p $T -c '$CUR_BAM' > '$gz'"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$gz" 2>/dev/null || echo 0)
  run_timed "$LOG" -- bash -c "pigz -d -p $T -c '$gz' > '$rt'"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU rtchk=MISMATCH
  [ "$(md5_of "$rt")" = "$(md5_of "$CUR_BAM")" ] && rtchk=OK
  emit_row "$RES" "$label" "yes" "$cbytes" "$CUR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rtchk"
  rm -f "$gz" "$rt"
}

# CRAM 3.1 with its strong codecs on (fqzcomp qualities, name tokeniser, range
# coder) — the fairest competitor; samtools does NOT enable these by default.
CRAM31_BEST=(--output-fmt-option use_fqz=1 --output-fmt-option use_tok=1 --output-fmt-option use_arith=1)

# (a) full WGS BAM — deferred headline (51 GB; needs ~200 GB scratch). DO_FULL=1.
if [ "$DO_FULL" = 1 ]; then
  prep_fingerprints "$BAM_FULL" "$REF_FULL"
  echo "bam full: $CUR_BAM bytes=$CUR_BYTES ref=$CUR_REF" | tee -a "$LOG"
  bench_bz  bz_full_default  yes
  bench_bz  bz_full_rq1      no   --reduce-quality 1
  bench_bz  bz_full_rq2      no   --reduce-quality 2
  bench_bz  bz_full_rq3      no   --reduce-quality 3
  bench_pigz pigz_full
  bench_cram cram30_full      3.0
  bench_cram cram31_full      3.1
  bench_cram cram31_best_full 3.1 "${CRAM31_BEST[@]}"
fi

# (b) chr20 BAM — the published tier: bz default + level sweep + lossy + CRAM + pigz.
prep_fingerprints "$BAM_CHR20" "$REF_CHR20"
echo "bam chr20: $CUR_BAM bytes=$CUR_BYTES ref=$CUR_REF" | tee -a "$LOG"
bench_bz  bz_chr20_default yes
bench_bz  bz_chr20_level1  yes  -l 1
bench_bz  bz_chr20_level2  yes  -l 2
bench_bz  bz_chr20_level3  yes  -l 3
bench_bz  bz_chr20_rq1     no   --reduce-quality 1
bench_bz  bz_chr20_rq2     no   --reduce-quality 2
bench_bz  bz_chr20_rq3     no   --reduce-quality 3
bench_pigz pigz_chr20
bench_cram cram30_chr20      3.0
bench_cram cram31_chr20      3.1
bench_cram cram31_best_chr20 3.1 "${CRAM31_BEST[@]}"

echo "DONE $(date)" >> "$LOG"
show_results "$RES"
rm -rf "$WORK"
