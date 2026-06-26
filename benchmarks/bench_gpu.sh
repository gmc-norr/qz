#!/usr/bin/env bash
# GPU vs CPU vs NUMA vs GPU+NUMA matrix (single-end FASTQ).
# GPU is compile-time (`--features cuda`), so this uses TWO binaries:
#   QZ      = CPU build      (target/release/qz)
#   GPU_QZ  = --features cuda (target-cuda/release/qz, links libcudart)
# Configs per (dataset, variant): CPU (--numa off), NUMA (--numa auto),
# GPU (--numa off), GPU+NUMA (--numa auto). Each compress+decompress, md5
# roundtrip. Output is byte-identical CPU vs GPU (GPU only changes the BWT
# backend), so ratios match; this measures the speed/RSS/CPU deltas.
#
# Override via env: SE5, SE30, QZ, GPU_QZ, CUDA_LIB, T, WORK, RES, LOG.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/lib/bench_common.sh"

SE5=${SE5:-$HERE/../real_data/ERR3239334_1.5m.fastq}
SE30=${SE30:-$HERE/../real_data/ERR3239334_1.30m.fastq}
QZ=${QZ:-$HERE/../target/release/qz}
GPU_QZ=${GPU_QZ:-$HERE/../target-cuda/release/qz}
CUDA_LIB=${CUDA_LIB:-/usr/local/cuda/lib64}
WORK=${WORK:-/tmp/qzbench_gpu}
RES=${RES:-$HERE/results/gpu.tsv}
LOG=${LOG:-$HERE/results/gpu.log}

export LD_LIBRARY_PATH="$CUDA_LIB:${LD_LIBRARY_PATH:-}"   # so the cuda binary + its numa re-execs find libcudart

rm -rf "$WORK"; mkdir -p "$WORK" "$(dirname "$RES")"
: > "$LOG"
tsv_header "$RES"

if [ ! -x "$GPU_QZ" ]; then echo "FATAL: GPU binary not found at $GPU_QZ (build: CUDA_PATH=/usr/local/cuda NVCC_APPEND_FLAGS=-arch=sm_75 CARGO_TARGET_DIR=target-cuda cargo build --release --features cuda -p qz-cli)" >&2; exit 2; fi
"$GPU_QZ" --version >/dev/null 2>&1 || { echo "FATAL: GPU binary won't run (libcudart missing? set CUDA_LIB)" >&2; exit 2; }

echo "gpu matrix: 5M=$SE5 30M=$SE30 T=$T $(date)" | tee -a "$LOG"
echo "GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null)" | tee -a "$LOG"

# run_cfg <bin> <numa off|auto> <variant-label> <input> <input_md5> <raw_bytes> <extra qz flags...>
run_cfg() {
  local bin="$1" numa="$2" label="$3" in="$4" md5="$5" raw="$6"; shift 6
  local arc="$WORK/$label.qz" dec="$WORK/$label.out.fastq"
  section "$label"
  run_timed "$LOG" -- "$bin" compress -i "$in" -o "$arc" -w "$WORK" -t "$T" -f --numa "$numa" "$@"
  local crc=$T_RC cs=$T_SECS crss=$T_RSS_KB ccpu=$T_CPU cbytes; cbytes=$(stat -c %s "$arc" 2>/dev/null || echo 0)
  if [ "$crc" -ne 0 ] || [ "$cbytes" -eq 0 ]; then
    emit_row "$RES" "$label" "yes" "0" "$raw" "$cs" "$crss" "$ccpu" "NA" "0" "NA" "COMPRESS_FAIL"; rm -f "$arc"; return; fi
  run_timed "$LOG" -- "$bin" decompress -i "$arc" -o "$dec" -w "$WORK" -t "$T" -f --numa "$numa"
  local ds=$T_SECS drss=$T_RSS_KB dcpu=$T_CPU
  [ "$(md5_of "$dec")" = "$md5" ] && local rt=OK || local rt=MISMATCH
  emit_row "$RES" "$label" "yes" "$cbytes" "$raw" "$cs" "$crss" "$ccpu" "$ds" "$drss" "$dcpu" "$rt"
  rm -f "$arc" "$dec"
}

# 4 configs for one (dataset,variant): cpu / numa / gpu / gpu+numa
matrix() { # tag input <extra flags...>
  local tag="$1" in="$2"; shift 2
  local md5 raw; md5=$(md5sum "$in" | awk '{print $1}'); raw=$(stat -c %s "$in")
  run_cfg "$QZ"     off  "${tag}_cpu"      "$in" "$md5" "$raw" "$@"
  run_cfg "$QZ"     auto "${tag}_numa"     "$in" "$md5" "$raw" "$@"
  run_cfg "$GPU_QZ" off  "${tag}_gpu"      "$in" "$md5" "$raw" "$@"
  run_cfg "$GPU_QZ" auto "${tag}_gpu_numa" "$in" "$md5" "$raw" "$@"
}

# 5M (shallow) — default + ultra-2
matrix 5m_default  "$SE5"
matrix 5m_ultra2   "$SE5"  --ultra 2
# 30M (deep single) — default + ultra-2
matrix 30m_default "$SE30"
matrix 30m_ultra2  "$SE30" --ultra 2

echo "DONE $(date)" >> "$LOG"
show_results "$RES"
rm -rf "$WORK"
