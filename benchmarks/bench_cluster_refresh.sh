#!/usr/bin/env bash
# Re-measure ONLY the cluster rows affected by the zstd worker cap-lift (16→48). Ratios +
# decompress are unchanged (cap-lift is compress-only, byte-identical) — this confirms the
# compress speedup at the SUMMARY/README configs. Uses the shared harness for identical framing.
# Override via env: QZ, S30M, R1, R2, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

QZ=${QZ:-$HERE/../target/release/qz}
S30M=${S30M:-$HERE/../real_data/ERR3239334_1.30m.fastq}
R1=${R1:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R1.fastq}
R2=${R2:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R2.fastq}
T=${T:-72}
WORK=${WORK:-/tmp/qzbench_clrefresh}
RES=${RES:-$HERE/results/cluster_refresh.tsv}
LOG=${LOG:-$HERE/results/cluster_refresh.log}

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"
echo "cluster-refresh: T=$T loadavg=$(cut -d' ' -f1-3 /proc/loadavg) $(date)" | tee -a "$LOG"

# 1) single-end 30M --cluster (matches bench_single_end.sh qz_cluster)
S30M_BYTES=$(stat -c %s "$S30M")
section "single 30M cluster"
arc="$WORK/s.qz"
run_timed "$LOG" -- "$QZ" compress -i "$S30M" -o "$arc" -w "$WORK" -t "$T" -f --numa off --cluster
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$WORK/s.out.fastq" -w "$WORK" -t "$T" -f
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
emit_row "$RES" "single30m_cluster" "set" "$cbytes" "$S30M_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(cluster_verify "$QZ" "$arc" "$T")"
rm -f "$arc" "$WORK/s.out.fastq"

# 2) paired chr20 --cluster (matches bench_paired_end.sh qz_cluster; two -o paths)
PAIR_BYTES=$(( $(stat -c %s "$R1") + $(stat -c %s "$R2") ))
section "paired chr20 cluster"
arc="$WORK/p.qz"
run_timed "$LOG" -- "$QZ" compress -i "$R1" -i "$R2" -o "$arc" -w "$WORK" -t "$T" -f --numa off --cluster
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$WORK/p_R1.fastq" -o "$WORK/p_R2.fastq" -w "$WORK" -t "$T" -f
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
emit_row "$RES" "paired_chr20_cluster" "set" "$cbytes" "$PAIR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(cluster_verify "$QZ" "$arc" "$T")"
rm -f "$arc" "$WORK/p_R1.fastq" "$WORK/p_R2.fastq"

# 3) cluster single-end on chr20 R1 (matches bench_reference.sh qz_cluster_single)
R1_BYTES=$(stat -c %s "$R1")
section "cluster single R1 chr20"
arc="$WORK/cs.qz"
run_timed "$LOG" -- "$QZ" compress -i "$R1" -o "$arc" -w "$WORK" -t "$T" -f --numa off --cluster
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$WORK/cs.out.fastq" -w "$WORK" -t "$T" -f
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
emit_row "$RES" "cluster_single_chr20R1" "set" "$cbytes" "$R1_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(cluster_verify "$QZ" "$arc" "$T")"
rm -f "$arc" "$WORK/cs.out.fastq"

echo "DONE loadavg=$(cut -d' ' -f1-3 /proc/loadavg) $(date)" | tee -a "$LOG"
show_results "$RES"
rm -rf "$WORK"
