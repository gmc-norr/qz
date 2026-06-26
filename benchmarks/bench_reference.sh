#!/usr/bin/env bash
# Reference-based FASTQ benchmark (deep tier): single-end + paired, qz reference/
# reference-fast + cluster, vs SPRING. Reference is compress-only.
#
# Override via env: INPUT_R1, INPUT_R2, REF, QZ, SPRING, SAM, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

INPUT_R1=${INPUT_R1:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R1.fastq}
INPUT_R2=${INPUT_R2:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R2.fastq}
REF=${REF:-$HERE/../real_data/chr20.fa}
QZ=${QZ:-$HERE/../target/release/qz}
SPRING=${SPRING:-spring}
SAM=${SAM:-samtools}
EXPECT_M5=${EXPECT_M5:-b18e6c531b0bd70e949a7fc20859cb01}
WORK=${WORK:-/tmp/qzbench_reference}
RES=${RES:-$HERE/results/reference.tsv}
LOG=${LOG:-$HERE/results/reference.log}

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"

R1_BYTES=$(stat -c %s "$INPUT_R1"); R2_BYTES=$(stat -c %s "$INPUT_R2")
ORIG_BYTES=$((R1_BYTES + R2_BYTES))
R1_MD5=$(md5sum "$INPUT_R1" | awk '{print $1}'); R2_MD5=$(md5sum "$INPUT_R2" | awk '{print $1}')
echo "reference: R1=$INPUT_R1 R2=$INPUT_R2 REF=$REF total_bytes=$ORIG_BYTES threads=$T $(date)" | tee -a "$LOG"

# Reference must match the reads' origin or mapping is meaningless.
REF_M5=$("$SAM" dict "$REF" 2>/dev/null | grep -m1 -oP 'SN:chr20\b.*M5:\K[0-9a-f]+' || echo "")
if [ "$REF_M5" != "$EXPECT_M5" ]; then
  echo "FATAL: $REF chr20 M5=$REF_M5 != expected $EXPECT_M5 (reference does not match the reads)" | tee -a "$LOG" >&2
  exit 2
fi
echo "reference identity OK: chr20 M5=$REF_M5" | tee -a "$LOG"

# Build the canonical index for BOTH variants, keyed to the actual read length
# via --like (151 -> canonical r150). require-index-by-default then guarantees
# compress finds it; we never hardcode the sidecar name.
section "index preflight"
"$QZ" index "$REF" --like "$INPUT_R1" -t "$T" 2>>"$LOG" || { echo "FATAL: index build failed" >&2; exit 2; }
"$QZ" index "$REF" --like "$INPUT_R1" --reference-fast -t "$T" 2>>"$LOG" || { echo "FATAL: index --reference-fast build failed" >&2; exit 2; }

pair_md5_ok() { [ "$(md5_of "$1")" = "$R1_MD5" ] && [ "$(md5_of "$2")" = "$R2_MD5" ] && echo OK || echo MISMATCH; }

# --- paired qz (denominator R1+R2) ---
bench_pair() { # name <extra flags...>
  local name="$1"; shift
  local arc="$WORK/$name.qz" pre="$WORK/$name.out"
  section "qz $name (paired)"
  run_timed "$LOG" -- "$QZ" compress -i "$INPUT_R1" -i "$INPUT_R2" -o "$arc" -w "$WORK" -t "$T" -f --numa off "$@"
  local crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  if [ "$crc" -ne 0 ] || [ "$cbytes" -eq 0 ]; then
    emit_row "$RES" "qz_${name}_paired" "yes" "0" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"; rm -f "$arc"; return; fi
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$pre" -w "$WORK" -t "$T" -f
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  emit_row "$RES" "qz_${name}_paired" "yes" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(pair_md5_ok "${pre}_R1.fastq" "${pre}_R2.fastq")"
  rm -f "$arc" "${pre}_R1.fastq" "${pre}_R2.fastq"
}

# --- single-end qz on R1 (denominator R1 only) ---
bench_single() { # name <extra flags...>
  local name="$1"; shift
  local arc="$WORK/$name.qz" dec="$WORK/$name.se.fastq"
  section "qz $name (single R1)"
  run_timed "$LOG" -- "$QZ" compress -i "$INPUT_R1" -o "$arc" -w "$WORK" -t "$T" -f --numa off "$@"
  local crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  if [ "$crc" -ne 0 ] || [ "$cbytes" -eq 0 ]; then
    emit_row "$RES" "qz_${name}_single" "yes" "0" "$R1_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"; rm -f "$arc"; return; fi
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(md5_of "$dec")" = "$R1_MD5" ] && local rt=OK || local rt=MISMATCH
  emit_row "$RES" "qz_${name}_single" "yes" "$cbytes" "$R1_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$arc" "$dec"
}

# single-end deep rows (R1)
section "qz cluster (single R1)"
arc="$WORK/cl_se.qz" dec="$WORK/cl_se.fastq"
run_timed "$LOG" -- "$QZ" compress -i "$INPUT_R1" -o "$arc" -w "$WORK" -t "$T" -f --numa off --cluster
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
emit_row "$RES" "qz_cluster_single" "set" "$cbytes" "$R1_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(cluster_verify "$QZ" "$arc" "$T")"
rm -f "$arc" "$dec"
bench_single reference       --reference "$REF"
bench_single reference_fast  --reference "$REF" --reference-fast

# paired deep rows (R1+R2)
bench_pair reference       --reference "$REF"
bench_pair reference_fast  --reference "$REF" --reference-fast

# SPRING baselines: single (R1) and paired
if need_tool "$SPRING" SPRING; then
  section "spring single (R1)"
  sp="$WORK/sp_se.spring" dec="$WORK/sp_se.fastq"
  run_timed "$LOG" -- "$SPRING" -c -i "$INPUT_R1" -o "$sp" -t "$T" -w "$WORK"
  cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$sp" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$SPRING" -d -i "$sp" -o "$dec" -t "$T" -w "$WORK"
  ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(md5_of "$dec")" = "$R1_MD5" ] && rt=OK || rt=MISMATCH
  emit_row "$RES" "spring_single" "yes" "$cbytes" "$R1_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$sp" "$dec"

  section "spring paired"
  sp="$WORK/sp.spring" d1="$WORK/sp_R1.fastq" d2="$WORK/sp_R2.fastq"
  run_timed "$LOG" -- "$SPRING" -c -i "$INPUT_R1" "$INPUT_R2" -o "$sp" -t "$T" -w "$WORK"
  cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$sp" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$SPRING" -d -i "$sp" -o "$d1" "$d2" -t "$T" -w "$WORK"
  ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  emit_row "$RES" "spring_paired" "yes" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(pair_md5_ok "$d1" "$d2")"
  rm -f "$sp" "$d1" "$d2"
fi

echo "DONE $(date)" >> "$LOG"
show_results "$RES"
rm -rf "$WORK"
