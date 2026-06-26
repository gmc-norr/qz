#!/usr/bin/env bash
# Fill-in rows for a like-for-like qz vs SPRING head-to-head. The main per-mode
# scripts already cover order-preserving qz vs SPRING (no -r) for single/reference
# and qz --cluster vs nothing for the reorder axis. This adds the missing halves:
#   - qz default PAIRED (order-preserving) to face SPRING-paired (no -r)
#   - SPRING -r (reorder; mates stay paired) for single 30M, single chr20 R1, and
#     paired chr20 — the reorder-axis competitor for qz --cluster.
# SPRING -r output is a permutation, so its roundtrip is a MULTISET check (sorted
# records), mirroring qz --cluster's qz-verify oracle (set-losslessness).
#
# Override via env: SE, R1, R2, QZ, SPRING, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

SE=${SE:-$HERE/../real_data/ERR3239334_1.30m.fastq}
R1=${R1:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R1.fastq}
R2=${R2:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R2.fastq}
QZ=${QZ:-$HERE/../target/release/qz}
SPRING=${SPRING:-spring}
WORK=${WORK:-/tmp/qzbench_springcmp}
RES=${RES:-$HERE/results/spring_compare.tsv}
LOG=${LOG:-$HERE/results/spring_compare.log}

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"

SE_BYTES=$(stat -c %s "$SE"); R1_BYTES=$(stat -c %s "$R1"); R2_BYTES=$(stat -c %s "$R2")
PAIR_BYTES=$((R1_BYTES + R2_BYTES))
R1_MD5=$(md5sum "$R1" | awk '{print $1}'); R2_MD5=$(md5sum "$R2" | awk '{print $1}')
echo "spring-compare: SE=$SE R1=$R1 R2=$R2 threads=$T $(date)" | tee -a "$LOG"

# order-independent multiset fingerprint (records as tab-joined 4-line groups)
ms_single() { paste - - - - < "$1" | LC_ALL=C sort --parallel="$T" -S 40G -T "$WORK" | md5sum | awk '{print $1}'; }
ms_pair()   { paste <(paste - - - - < "$1") <(paste - - - - < "$2") | LC_ALL=C sort --parallel="$T" -S 40G -T "$WORK" | md5sum | awk '{print $1}'; }

section "precompute input multisets"
MS_SE=$(ms_single "$SE");  echo "MS_SE=$MS_SE"   | tee -a "$LOG"
MS_R1=$(ms_single "$R1");  echo "MS_R1=$MS_R1"   | tee -a "$LOG"
MS_PAIR=$(ms_pair "$R1" "$R2"); echo "MS_PAIR=$MS_PAIR" | tee -a "$LOG"

# ---- (1) qz default PAIRED (order-preserving) ----
section "qz default (paired, order-preserving)"
arc="$WORK/qzpe.qz" pre="$WORK/qzpe.out"
run_timed "$LOG" -- "$QZ" compress -i "$R1" -i "$R2" -o "$arc" -w "$WORK" -t "$T" -f --numa off
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$pre" -w "$WORK" -t "$T" -f
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
[ "$(md5_of "${pre}_R1.fastq")" = "$R1_MD5" ] && [ "$(md5_of "${pre}_R2.fastq")" = "$R2_MD5" ] && rt=OK || rt=MISMATCH
emit_row "$RES" "qz_default_paired" "yes" "$cbytes" "$PAIR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
rm -f "$arc" "${pre}_R1.fastq" "${pre}_R2.fastq"

# ---- SPRING -r helper (single): reorder, multiset roundtrip ----
spring_r_single() { # label input raw_bytes input_ms
  local label="$1" in="$2" raw="$3" ms="$4"
  section "spring -r $label"
  local sp="$WORK/$label.spring" dec="$WORK/$label.out.fastq"
  run_timed "$LOG" -- "$SPRING" -c -r -i "$in" -o "$sp" -t "$T" -w "$WORK"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$sp" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$SPRING" -d -i "$sp" -o "$dec" -t "$T" -w "$WORK"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(ms_single "$dec")" = "$ms" ] && local rt=PASS || local rt=FAIL
  emit_row "$RES" "$label" "set" "$cbytes" "$raw" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$sp" "$dec"
}

# ---- (2) SPRING -r single, 30M ----
spring_r_single spring_reorder_single_30m "$SE" "$SE_BYTES" "$MS_SE"
# ---- (3) SPRING -r single, chr20 R1 ----
spring_r_single spring_reorder_single_chr20 "$R1" "$R1_BYTES" "$MS_R1"

# ---- (4) SPRING -r paired, chr20 ----
section "spring -r paired (chr20)"
sp="$WORK/sp_r_pe.spring" d1="$WORK/sp_r_R1.fastq" d2="$WORK/sp_r_R2.fastq"
run_timed "$LOG" -- "$SPRING" -c -r -i "$R1" "$R2" -o "$sp" -t "$T" -w "$WORK"
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes=$(stat -c %s "$sp" 2>/dev/null || echo 0)
run_timed "$LOG" -- "$SPRING" -d -i "$sp" -o "$d1" "$d2" -t "$T" -w "$WORK"
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
[ "$(ms_pair "$d1" "$d2")" = "$MS_PAIR" ] && rt=PASS || rt=FAIL
emit_row "$RES" "spring_reorder_paired_chr20" "set" "$cbytes" "$PAIR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
rm -f "$sp" "$d1" "$d2"

echo "DONE $(date)" >> "$LOG"
show_results "$RES"
rm -rf "$WORK"
