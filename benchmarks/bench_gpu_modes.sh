#!/usr/bin/env bash
# GPU × NUMA across the OTHER modes: bz (BAM), paired --cluster, and reference.
# Two binaries each (CPU build vs --features cuda build). Configs per workload:
# CPU (--numa off), NUMA (--numa auto), GPU (--numa off), GPU+NUMA (--numa auto).
#
# Note: --cluster's sequence stream is zstd-long (not BSC), quality is fqz, headers
# columnar — only a tiny StrandBits stream is BSC. So the GPU (which only accelerates
# libbsc's BWT) is expected to do ~nothing for cluster; this measures that.
# Reference DOES use BSC for its big streams (consensus backing, fallback literals,
# R2 delta), so the GPU engages there — but reference compress is serialization-bound.
#
# Override via env: R1, R2, REF, BAM, QZ, GPU_QZ, BZ, GPU_BZ, SAM, CUDA_LIB, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

R1=${R1:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R1.fastq}
R2=${R2:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R2.fastq}
REF=${REF:-$HERE/../real_data/chr20.fa}
BAM=${BAM:-$HERE/../real_data/HG002.GRCh38.chr20.bam}
QZ=${QZ:-$HERE/../target/release/qz}
GPU_QZ=${GPU_QZ:-$HERE/../target-cuda/release/qz}
BZ=${BZ:-$HERE/../target/release/bz}
GPU_BZ=${GPU_BZ:-$HERE/../target-cuda/release/bz}
SAM=${SAM:-samtools}
CUDA_LIB=${CUDA_LIB:-/usr/local/cuda/lib64}
WORK=${WORK:-/tmp/qzbench_gpumodes}
RES=${RES:-$HERE/results/gpu_modes.tsv}
LOG=${LOG:-$HERE/results/gpu_modes.log}

export LD_LIBRARY_PATH="$CUDA_LIB:${LD_LIBRARY_PATH:-}"

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"

R1_BYTES=$(stat -c %s "$R1"); R2_BYTES=$(stat -c %s "$R2"); PAIR_BYTES=$((R1_BYTES + R2_BYTES))
R1_MD5=$(md5sum "$R1" | awk '{print $1}'); R2_MD5=$(md5sum "$R2" | awk '{print $1}')
pair_md5_ok() { [ "$(md5_of "$1")" = "$R1_MD5" ] && [ "$(md5_of "$2")" = "$R2_MD5" ] && echo OK || echo MISMATCH; }
echo "gpu-modes: T=$T $(date)" | tee -a "$LOG"

# reference index (untimed)
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

# ---- bz (BAM): label binary numa ----
bz_cfg() { # label bin numa
  local label="$1" bin="$2" numa="$3"
  local arc="$WORK/$label.bz" rt="$WORK/$label.rt.bam"
  section "$label"
  run_timed "$LOG" -- "$bin" compress -i "$BAM_C" -o "$arc" -t "$T" -f --numa "$numa"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$bin" decompress -i "$arc" -o "$rt" -t "$T" -f --numa "$numa"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU m
  m=$("$SAM" view -@ "$T" "$rt" 2>>"$LOG" | md5sum | awk '{print $1}')
  [ "$m" = "$BAM_FP" ] && local rt2=OK || local rt2=MISMATCH
  emit_row "$RES" "$label" "yes" "$cbytes" "$BAM_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt2"
  rm -f "$arc" "$rt"
}
bz_cfg bz_cpu      "$BZ"     off
bz_cfg bz_numa     "$BZ"     auto
bz_cfg bz_gpu      "$GPU_BZ" off
bz_cfg bz_gpu_numa "$GPU_BZ" auto

# ---- paired --cluster: label bin numa (qz verify oracle) ----
cl_cfg() { # label bin numa
  local label="$1" bin="$2" numa="$3"
  local arc="$WORK/$label.qz"
  section "$label"
  run_timed "$LOG" -- "$bin" compress -i "$R1" -i "$R2" -o "$arc" -w "$WORK" -t "$T" -f --numa "$numa" --cluster
  local crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  if [ "$crc" -ne 0 ] || [ "$cbytes" -eq 0 ]; then
    emit_row "$RES" "$label" "set" "0" "$PAIR_BYTES" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"; rm -f "$arc"; return; fi
  run_timed "$LOG" -- "$bin" decompress -i "$arc" -o "$WORK/${label}_R1.fastq" -o "$WORK/${label}_R2.fastq" -w "$WORK" -t "$T" -f --numa "$numa"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  emit_row "$RES" "$label" "set" "$cbytes" "$PAIR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(cluster_verify "$bin" "$arc" "$T")"
  rm -f "$arc" "$WORK/${label}_R1.fastq" "$WORK/${label}_R2.fastq"
}
cl_cfg cluster_paired_cpu      "$QZ"     off
cl_cfg cluster_paired_numa     "$QZ"     auto
cl_cfg cluster_paired_gpu      "$GPU_QZ" off
cl_cfg cluster_paired_gpu_numa "$GPU_QZ" auto

# ---- reference single (R1) + paired (R1+R2): label bin numa paired? ----
ref_single_cfg() { # label bin numa
  local label="$1" bin="$2" numa="$3"
  local arc="$WORK/$label.qz" dec="$WORK/$label.out.fastq"
  section "$label"
  run_timed "$LOG" -- "$bin" compress -i "$R1" --reference "$REF" -o "$arc" -w "$WORK" -t "$T" -f --numa "$numa"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$bin" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f --numa "$numa"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(md5_of "$dec")" = "$R1_MD5" ] && local rt=OK || local rt=MISMATCH
  emit_row "$RES" "$label" "yes" "$cbytes" "$R1_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$arc" "$dec"
}
ref_paired_cfg() { # label bin numa
  local label="$1" bin="$2" numa="$3"
  local arc="$WORK/$label.qz" pre="$WORK/$label.out"
  section "$label"
  run_timed "$LOG" -- "$bin" compress -i "$R1" -i "$R2" --reference "$REF" -o "$arc" -w "$WORK" -t "$T" -f --numa "$numa"
  local cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  run_timed "$LOG" -- "$bin" decompress -i "$arc" -o "$pre" -w "$WORK" -t "$T" -f --numa "$numa"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  emit_row "$RES" "$label" "yes" "$cbytes" "$PAIR_BYTES" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$(pair_md5_ok "${pre}_R1.fastq" "${pre}_R2.fastq")"
  rm -f "$arc" "${pre}_R1.fastq" "${pre}_R2.fastq"
}
ref_single_cfg ref_single_cpu      "$QZ"     off
ref_single_cfg ref_single_numa     "$QZ"     auto
ref_single_cfg ref_single_gpu      "$GPU_QZ" off
ref_single_cfg ref_single_gpu_numa "$GPU_QZ" auto
ref_paired_cfg ref_paired_cpu      "$QZ"     off
ref_paired_cfg ref_paired_numa     "$QZ"     auto
ref_paired_cfg ref_paired_gpu      "$GPU_QZ" off
ref_paired_cfg ref_paired_gpu_numa "$GPU_QZ" auto

echo "DONE $(date)" >> "$LOG"
show_results "$RES"
rm -rf "$WORK"
