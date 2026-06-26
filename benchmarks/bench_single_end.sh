#!/usr/bin/env bash
# Single-end FASTQ benchmark: qz (all production settings) vs SPRING.
# Metrics per variant: compressed size, ratio, compress/decompress wall + peak
# RSS + CPU%, and md5 roundtrip (lossless variants only).
#
# Override via env: INPUT, QZ, SPRING, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

INPUT=${INPUT:-$HERE/../real_data/ERR3239334_1.30m.fastq}
QZ=${QZ:-$HERE/../target/release/qz}
SPRING=${SPRING:-spring}
WORK=${WORK:-/tmp/qzbench_single}
RES=${RES:-$HERE/results/single.tsv}
LOG=${LOG:-$HERE/results/single.log}

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"

ORIG_BYTES=$(stat -c %s "$INPUT")
ORIG_MD5=$(md5sum "$INPUT" | awk '{print $1}')
echo "single-end: INPUT=$INPUT bytes=$ORIG_BYTES md5=$ORIG_MD5 threads=$T $(date)" | tee -a "$LOG"

# bench_qz <variant> <lossless yes|no> <extra qz compress flags...>
bench_qz() {
  local name="$1" lossless="$2"; shift 2
  local arc="$WORK/$name.qz" dec="$WORK/$name.out.fastq"
  section "qz $name"
  run_timed "$LOG" -- "$QZ" compress -i "$INPUT" -o "$arc" -w "$WORK" -t "$T" -f --numa off "$@"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU rt="n-a"
  if [ "$lossless" = "yes" ]; then
    [ "$(md5_of "$dec")" = "$ORIG_MD5" ] && rt="OK" || rt="MISMATCH"
  fi
  emit_row "$RES" "qz_$name" "$lossless" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$arc" "$dec"
}

bench_qz default  yes
bench_qz fast     yes  --fast
bench_qz ultra1   yes  --ultra 1
bench_qz ultra2   yes  --ultra 2
bench_qz ultra3   yes  --ultra 3

# cluster: output is a permutation -> qz verify oracle (not md5)
section "qz cluster"
arc="$WORK/cluster.qz"
run_timed "$LOG" -- "$QZ" compress -i "$INPUT" -o "$arc" -w "$WORK" -t "$T" -f --numa off --cluster
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
dec="$WORK/cluster.out.fastq"
run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
emit_row "$RES" "qz_cluster" "set" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(cluster_verify "$QZ" "$arc" "$T")"
rm -f "$arc" "$dec"

# NUMA on-vs-off (default codec): the sharding wall lever
section "qz default --numa auto"
arc="$WORK/numa_auto.qz" dec="$WORK/numa_auto.out.fastq"
run_timed "$LOG" -- "$QZ" compress -i "$INPUT" -o "$arc" -w "$WORK" -t "$T" -f --numa auto
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f --numa auto
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
[ "$(md5_of "$dec")" = "$ORIG_MD5" ] && rt=OK || rt=MISMATCH
emit_row "$RES" "qz_default_numa_auto" "yes" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
rm -f "$arc" "$dec"

# ---- SPRING (lossless, order-preserving = no -r) ----
if need_tool "$SPRING" SPRING; then
  section "spring lossless"
  sp="$WORK/spring.spring" dec="$WORK/spring.out.fastq"
  run_timed "$LOG" -- "$SPRING" -c -i "$INPUT" -o "$sp" -t "$T" -w "$WORK"
  cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$sp" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$SPRING" -d -i "$sp" -o "$dec" -t "$T" -w "$WORK"
  ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(md5_of "$dec")" = "$ORIG_MD5" ] && rt=OK || rt=MISMATCH
  emit_row "$RES" "spring_lossless" "yes" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$sp" "$dec"
fi

# ---- pigz -9 (general-purpose baseline) ----
if need_tool pigz pigz; then
  section "pigz -9"
  gz="$WORK/pigz.gz" dec="$WORK/pigz.out"
  run_timed "$LOG" -- bash -c "pigz -9 -p $T -c '$INPUT' > '$gz'"
  cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$gz" 2>/dev/null || echo 0)
  run_timed "$LOG" -- bash -c "pigz -d -p $T -c '$gz' > '$dec'"
  ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(md5_of "$dec")" = "$ORIG_MD5" ] && rt=OK || rt=MISMATCH
  emit_row "$RES" "pigz_9" "yes" "$cbytes" "$ORIG_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$gz" "$dec"
fi

echo "DONE $(date)" >> "$LOG"
show_results "$RES"
rm -rf "$WORK"
