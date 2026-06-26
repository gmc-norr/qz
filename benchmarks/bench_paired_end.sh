#!/usr/bin/env bash
# Paired-end FASTQ benchmark: qz (paired) vs SPRING (paired).
# qz paired decompress: -o <prefix> emits <prefix>_R1.fastq + <prefix>_R2.fastq.
# Roundtrip = md5 of BOTH mates vs the originals.
#
# Override via env: INPUT_R1, INPUT_R2, QZ, SPRING, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

INPUT_R1=${INPUT_R1:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R1.fastq}
INPUT_R2=${INPUT_R2:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R2.fastq}
QZ=${QZ:-$HERE/../target/release/qz}
SPRING=${SPRING:-spring}
WORK=${WORK:-/tmp/qzbench_paired}
RES=${RES:-$HERE/results/paired.tsv}
LOG=${LOG:-$HERE/results/paired.log}

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"

R1_BYTES=$(stat -c %s "$INPUT_R1"); R2_BYTES=$(stat -c %s "$INPUT_R2")
ORIG_BYTES=$((R1_BYTES + R2_BYTES))
R1_MD5=$(md5sum "$INPUT_R1" | awk '{print $1}'); R2_MD5=$(md5sum "$INPUT_R2" | awk '{print $1}')
echo "paired: R1=$INPUT_R1 R2=$INPUT_R2 total_bytes=$ORIG_BYTES threads=$T $(date)" | tee -a "$LOG"

pair_md5_ok() { # decR1 decR2  -> echoes OK|MISMATCH
  [ "$(md5_of "$1")" = "$R1_MD5" ] && [ "$(md5_of "$2")" = "$R2_MD5" ] && echo OK || echo MISMATCH
}

# bench_qz <variant> <extra qz compress flags...>   (always lossless here)
bench_qz() {
  local name="$1"; shift
  local arc="$WORK/$name.qz" pre="$WORK/$name.out"
  section "qz $name"
  run_timed "$LOG" -- "$QZ" compress -i "$INPUT_R1" -i "$INPUT_R2" -o "$arc" -w "$WORK" -t "$T" -f "$@"
  local crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  if [ "$crc" -ne 0 ] || [ "$cbytes" -eq 0 ]; then
    emit_row "$RES" "qz_$name" "yes" "0" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"
    rm -f "$arc"; return
  fi
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$pre" -w "$WORK" -t "$T" -f
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  local rt; rt=$(pair_md5_ok "${pre}_R1.fastq" "${pre}_R2.fastq")
  emit_row "$RES" "qz_$name" "yes" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$arc" "${pre}_R1.fastq" "${pre}_R2.fastq"
}

# Deep paired: cluster (order-drop). Output is a permutation -> qz verify oracle.
# Paired-cluster decode takes TWO exact -o paths (not the prefix contract).
section "qz cluster (paired)"
arc="$WORK/cluster.qz"
run_timed "$LOG" -- "$QZ" compress -i "$INPUT_R1" -i "$INPUT_R2" -o "$arc" -w "$WORK" -t "$T" -f --numa off --cluster
crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
if [ "$crc" -eq 0 ] && [ "$cbytes" -gt 0 ]; then
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$WORK/cl_R1.fastq" -o "$WORK/cl_R2.fastq" -w "$WORK" -t "$T" -f
  ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  emit_row "$RES" "qz_cluster" "set" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(cluster_verify "$QZ" "$arc" "$T")"
else
  emit_row "$RES" "qz_cluster" "set" "0" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"
fi
rm -f "$arc" "$WORK/cl_R1.fastq" "$WORK/cl_R2.fastq"

# ---- SPRING paired (lossless, order-preserving) ----
if need_tool "$SPRING" SPRING; then
  section "spring paired"
  sp="$WORK/spring.spring" d1="$WORK/sp_R1.fastq" d2="$WORK/sp_R2.fastq"
  run_timed "$LOG" -- "$SPRING" -c -i "$INPUT_R1" "$INPUT_R2" -o "$sp" -t "$T" -w "$WORK"
  crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$sp" 2>/dev/null || echo 0)
  if [ "$crc" -eq 0 ] && [ "$cbytes" -gt 0 ]; then
    run_timed "$LOG" -- "$SPRING" -d -i "$sp" -o "$d1" "$d2" -t "$T" -w "$WORK"
    ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
    rt=$(pair_md5_ok "$d1" "$d2")
    emit_row "$RES" "spring_paired" "yes" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  else
    emit_row "$RES" "spring_paired" "yes" "0" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"
  fi
  rm -f "$sp" "$d1" "$d2"
fi

echo "DONE $(date)" >> "$LOG"
show_results "$RES"
rm -rf "$WORK"
