#!/usr/bin/env bash
# Measure pigz -9 general-purpose baselines on the datasets that lacked a pigz row,
# matching bench_single_end.sh's pigz methodology (`pigz -9 -p T -c`) so they slot
# directly into the Performance tables. The single-end 30M tier already carries its
# pigz row (5.03×) from the main 2026-06-24 pass and is left as-is.
#
# Rows produced:
#   pigz_paired_chr20   — R1+R2 (7.47 GB); also the reference-paired baseline (same raw reads)
#   pigz_ref_single_R1  — R1 alone (3.74 GB); the reference-single baseline
#   pigz_bam_chr20      — re-gzip of the reheadered BGZF BAM (same input bz/CRAM use)
#
# Override via env: R1, R2, BAM, SAM, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

R1=${R1:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R1.fastq}
R2=${R2:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R2.fastq}
BAM=${BAM:-$HERE/../real_data/HG002.GRCh38.chr20.bam}
SAM=${SAM:-samtools}
T=${T:-72}
WORK=${WORK:-/tmp/qzbench_pigz}
RES=${RES:-$HERE/results/pigz_all.tsv}
LOG=${LOG:-$HERE/results/pigz_all.log}

need_tool pigz pigz || { echo "pigz required" >&2; exit 1; }
rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"
echo "pigz-all: T=$T $(pigz --version 2>&1 | head -1) loadavg=$(cut -d' ' -f1-3 /proc/loadavg) $(date)" | tee -a "$LOG"

# 1) paired chr20 (R1+R2) — gzip each mate separately (how you'd store a pair), sum sizes.
section "pigz -9 paired (R1+R2)"
R1B=$(stat -c %s "$R1"); R2B=$(stat -c %s "$R2"); PAIRB=$((R1B + R2B))
m1=$(md5_of "$R1"); m2=$(md5_of "$R2")
g1="$WORK/r1.gz" g2="$WORK/r2.gz" d1="$WORK/r1.out" d2="$WORK/r2.out"
run_timed "$LOG" -- bash -c "pigz -9 -p $T -c '$R1' > '$g1' && pigz -9 -p $T -c '$R2' > '$g2'"
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU
cbytes=$(( $(stat -c %s "$g1") + $(stat -c %s "$g2") ))
run_timed "$LOG" -- bash -c "pigz -d -p $T -c '$g1' > '$d1' && pigz -d -p $T -c '$g2' > '$d2'"
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
{ [ "$(md5_of "$d1")" = "$m1" ] && [ "$(md5_of "$d2")" = "$m2" ]; } && rt=OK || rt=MISMATCH
emit_row "$RES" "pigz_paired_chr20" "yes" "$cbytes" "$PAIRB" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
rm -f "$g1" "$g2" "$d1" "$d2"

# 2) reference-single baseline = pigz on R1 alone (3.74 GB).
section "pigz -9 single R1 (3.74 GB)"
g="$WORK/r1only.gz" d="$WORK/r1only.out"
run_timed "$LOG" -- bash -c "pigz -9 -p $T -c '$R1' > '$g'"
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU; cbytes=$(stat -c %s "$g")
run_timed "$LOG" -- bash -c "pigz -d -p $T -c '$g' > '$d'"
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
[ "$(md5_of "$d")" = "$m1" ] && rt=OK || rt=MISMATCH
emit_row "$RES" "pigz_ref_single_R1" "yes" "$cbytes" "$R1B" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
rm -f "$g" "$d"

# 3) BAM chr20 — reheader header-only to SO:coordinate (identical to bench_bam.sh, so the
#    denominator matches the bz/CRAM rows), then re-gzip the BGZF bytes. A BAM is already
#    DEFLATE-compressed (BGZF), so this generic baseline barely shrinks it (~1.0×) — the point.
section "pigz -9 BAM (reheadered; re-gzip of BGZF)"
RBAM="$WORK/chr20.coord.bam"
"$SAM" view -H "$BAM" 2>>"$LOG" | sed 's/SO:unsorted/SO:coordinate/' > "$WORK/hdr.sam"
"$SAM" reheader "$WORK/hdr.sam" "$BAM" > "$RBAM" 2>>"$LOG"
BAMB=$(stat -c %s "$RBAM"); mB=$(md5_of "$RBAM")
echo "reheadered BAM bytes=$BAMB" | tee -a "$LOG"
g="$WORK/bam.gz" d="$WORK/bam.out"
run_timed "$LOG" -- bash -c "pigz -9 -p $T -c '$RBAM' > '$g'"
cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU; cbytes=$(stat -c %s "$g")
run_timed "$LOG" -- bash -c "pigz -d -p $T -c '$g' > '$d'"
ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
[ "$(md5_of "$d")" = "$mB" ] && rt=OK || rt=MISMATCH
emit_row "$RES" "pigz_bam_chr20" "yes" "$cbytes" "$BAMB" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
rm -f "$g" "$d" "$RBAM" "$WORK/hdr.sam"

echo "DONE loadavg=$(cut -d' ' -f1-3 /proc/loadavg) $(date)" | tee -a "$LOG"
show_results "$RES"
rm -rf "$WORK"
