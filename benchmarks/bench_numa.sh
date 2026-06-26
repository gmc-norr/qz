#!/usr/bin/env bash
# NUMA on-vs-off sweep. For each production mode, ONE workload is compressed AND
# decoded at --numa off, then again at --numa auto. Rows come in pairs
# (<mode>_numa_off vs <mode>_numa_auto); the delta isolates the multi-socket
# sharding win on BOTH compress and decompress. Decode is the bigger NUMA lever
# (~1.6-1.9x per the design notes), so it is measured for every mode here — the
# per-mode scripts only pinned compress and let decode default to --numa auto.
#
# Override via env: SE, R1, R2, REF, QZ, BZ, SAM, BAM, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

SE=${SE:-$HERE/../real_data/ERR3239334_1.30m.fastq}
R1=${R1:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R1.fastq}
R2=${R2:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R2.fastq}
REF=${REF:-$HERE/../real_data/chr20.fa}
QZ=${QZ:-$HERE/../target/release/qz}
BZ=${BZ:-$HERE/../target/release/bz}
SAM=${SAM:-samtools}
BAM=${BAM:-$HERE/../real_data/HG002.GRCh38.chr20.bam}
WORK=${WORK:-/tmp/qzbench_numa}
RES=${RES:-$HERE/results/numa.tsv}
LOG=${LOG:-$HERE/results/numa.log}

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"

SE_BYTES=$(stat -c %s "$SE"); SE_MD5=$(md5sum "$SE" | awk '{print $1}')
R1_BYTES=$(stat -c %s "$R1"); R2_BYTES=$(stat -c %s "$R2"); PAIR_BYTES=$((R1_BYTES + R2_BYTES))
R1_MD5=$(md5sum "$R1" | awk '{print $1}'); R2_MD5=$(md5sum "$R2" | awk '{print $1}')
pair_md5_ok() { [ "$(md5_of "$1")" = "$R1_MD5" ] && [ "$(md5_of "$2")" = "$R2_MD5" ] && echo OK || echo MISMATCH; }
echo "numa sweep: threads=$T $(date)" | tee -a "$LOG"

# reference index (untimed) so reference rows find it
"$QZ" index "$REF" --like "$R1" -t "$T" 2>>"$LOG" || { echo "FATAL: index build failed" >&2; exit 2; }

# coordinate-labeled BAM for bz
SO=$("$SAM" view -H "$BAM" 2>/dev/null | grep -m1 -oP 'SO:\K[a-z]+' || echo "")
BAM_C="$BAM"
if [ "$SO" != coordinate ]; then
  BAM_C="$WORK/bam.coord.bam"
  "$SAM" view -H "$BAM" 2>/dev/null | sed 's/SO:unsorted/SO:coordinate/' > "$WORK/hdr.sam"
  "$SAM" reheader "$WORK/hdr.sam" "$BAM" > "$BAM_C" 2>>"$LOG"
fi
BAM_BYTES=$(stat -c %s "$BAM_C")
BAM_FP=$("$SAM" view -@ "$T" "$BAM_C" 2>>"$LOG" | md5sum | awk '{print $1}')

# qz single-end (30M): compress+decode at numa mode $1
numa_single() { # mode
  local M="$1" arc="$WORK/se.qz" dec="$WORK/se.out.fastq"
  section "single numa=$M"
  run_timed "$LOG" -- "$QZ" compress -i "$SE" -o "$arc" -w "$WORK" -t "$T" -f --numa "$M"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f --numa "$M"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(md5_of "$dec")" = "$SE_MD5" ] && local rt=OK || local rt=MISMATCH
  emit_row "$RES" "single_numa_$M" "yes" "$cbytes" "$SE_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$arc" "$dec"
}

# qz paired (chr20): compress+decode at numa mode $1
numa_paired() { # mode
  local M="$1" arc="$WORK/pe.qz" pre="$WORK/pe.out"
  section "paired numa=$M"
  run_timed "$LOG" -- "$QZ" compress -i "$R1" -i "$R2" -o "$arc" -w "$WORK" -t "$T" -f --numa "$M"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$pre" -w "$WORK" -t "$T" -f --numa "$M"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  emit_row "$RES" "paired_numa_$M" "yes" "$cbytes" "$PAIR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(pair_md5_ok "${pre}_R1.fastq" "${pre}_R2.fastq")"
  rm -f "$arc" "${pre}_R1.fastq" "${pre}_R2.fastq"
}

# qz reference single (type 4, chr20 R1): compress+decode at numa mode $1
numa_ref_single() { # mode
  local M="$1" arc="$WORK/rs.qz" dec="$WORK/rs.out.fastq"
  section "reference_single numa=$M"
  run_timed "$LOG" -- "$QZ" compress -i "$R1" --reference "$REF" -o "$arc" -w "$WORK" -t "$T" -f --numa "$M"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f --numa "$M"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(md5_of "$dec")" = "$R1_MD5" ] && local rt=OK || local rt=MISMATCH
  emit_row "$RES" "reference_single_numa_$M" "yes" "$cbytes" "$R1_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$arc" "$dec"
}

# qz reference paired (type 2, chr20): compress+decode at numa mode $1
numa_ref_paired() { # mode
  local M="$1" arc="$WORK/rp.qz" pre="$WORK/rp.out"
  section "reference_paired numa=$M"
  run_timed "$LOG" -- "$QZ" compress -i "$R1" -i "$R2" --reference "$REF" -o "$arc" -w "$WORK" -t "$T" -f --numa "$M"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$QZ" decompress -i "$arc" -o "$pre" -w "$WORK" -t "$T" -f --numa "$M"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  emit_row "$RES" "reference_paired_numa_$M" "yes" "$cbytes" "$PAIR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(pair_md5_ok "${pre}_R1.fastq" "${pre}_R2.fastq")"
  rm -f "$arc" "${pre}_R1.fastq" "${pre}_R2.fastq"
}

# bz (chr20 BAM): compress+decode at numa mode $1
numa_bz() { # mode
  local M="$1" arc="$WORK/bz.bz" rt="$WORK/bz.rt.bam"
  section "bz numa=$M"
  run_timed "$LOG" -- "$BZ" compress -i "$BAM_C" -o "$arc" -t "$T" -f --numa "$M"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$BZ" decompress -i "$arc" -o "$rt" -t "$T" -f --numa "$M"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU m
  m=$("$SAM" view -@ "$T" "$rt" 2>>"$LOG" | md5sum | awk '{print $1}')
  [ "$m" = "$BAM_FP" ] && local rt2=OK || local rt2=MISMATCH
  emit_row "$RES" "bz_numa_$M" "yes" "$cbytes" "$BAM_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt2"
  rm -f "$arc" "$rt"
}

for fn in numa_single numa_paired numa_ref_single numa_ref_paired numa_bz; do
  "$fn" off
  "$fn" auto
done

echo "DONE $(date)" >> "$LOG"
show_results "$RES"
rm -rf "$WORK"
