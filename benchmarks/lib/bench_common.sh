#!/usr/bin/env bash
# Shared helpers for the qz/bz end-to-end mode benchmarks.
# Source this from a per-mode script:  source "$(dirname "$0")/lib/bench_common.sh"
#
# Provides:
#   run_timed <logfile> -- cmd...   -> sets globals T_SECS, T_RSS_KB, T_RSS_GB, T_CPU, T_RC
#   ratio <orig_bytes> <comp_bytes> -> prints "X.XX" (or "NA")
#   pct_of <comp_bytes> <orig_bytes>-> prints "NN.N" (comp as % of orig)
#   md5_of <file>                   -> prints md5 hex (or "MISSING")
#   tsv_header <resultfile>         -> writes the standard column header
#   emit_row <resultfile> <variant> <lossless> <comp_bytes> <orig_bytes> \
#            <comp_s> <comp_rss_kb> <comp_cpu> <decomp_s> <decomp_rss_kb> \
#            <decomp_cpu> <roundtrip>
#   need_tool <path> <label>        -> warns + returns 1 if missing/not executable
#   section <title>                 -> log banner
#
# Conventions (verified against the binaries):
#   * qz/bz refuse to overwrite an existing -o; always pass -f or clean first.
#   * paired/reference decompress: -o <prefix> emits <prefix>_R1.fastq + _R2.fastq.
#   * bz roundtrip is byte-identical at the RECORD level, not the BGZF-byte level;
#     compare with `samtools view` diff, not md5 (see bench_bam.sh).
set -u

: "${T:=72}"                      # threads (override: T=NN)
: "${TIME_BIN:=/usr/bin/time}"    # GNU time (the `-v` form)

# run_timed <logfile> -- cmd...
# Appends the command's stdout+the time report to <logfile>.
# Sets: T_SECS (wall, float seconds), T_RSS_KB, T_RSS_GB, T_CPU ("NNN%"), T_RC (exit code).
run_timed() {
  local logf="$1"; shift
  [ "${1:-}" = "--" ] && shift
  local tf
  tf="$(mktemp)"
  "$TIME_BIN" -v "$@" >>"$logf" 2>"$tf"
  T_RC=$?
  T_SECS=$(grep "Elapsed (wall clock)" "$tf" | sed -E 's/.*: //' \
            | awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1 }')
  T_RSS_KB=$(grep "Maximum resident set size" "$tf" | sed -E 's/.*: //')
  T_CPU=$(grep "Percent of CPU" "$tf" | sed -E 's/.*: //')
  : "${T_SECS:=NA}" "${T_RSS_KB:=0}" "${T_CPU:=NA}"
  T_RSS_GB=$(awk -v k="$T_RSS_KB" 'BEGIN{ printf "%.2f", k/1048576 }')
  cat "$tf" >>"$logf"
  rm -f "$tf"
  return "$T_RC"
}

ratio() { # orig comp
  awk -v o="$1" -v c="$2" 'BEGIN{ if (c>0) printf "%.2f", o/c; else print "NA" }'
}

pct_of() { # comp orig
  awk -v c="$1" -v o="$2" 'BEGIN{ if (o>0) printf "%.1f", 100*c/o; else print "NA" }'
}

md5_of() { # file
  if [ -f "$1" ]; then md5sum "$1" | awk '{print $1}'; else echo "MISSING"; fi
}

tsv_header() { # resultfile
  printf 'variant\tlossless\tcomp_bytes\traw_bytes\tratio\tcomp_s\tcomp_rss_gb\tcomp_cpu\tdecomp_s\tdecomp_rss_gb\tdecomp_cpu\troundtrip\n' > "$1"
}

# emit_row <rf> <variant> <lossless> <comp_bytes> <orig_bytes> <comp_s> <comp_rss_kb> <comp_cpu> <decomp_s> <decomp_rss_kb> <decomp_cpu> <roundtrip>
emit_row() {
  local rf="$1" variant="$2" lossless="$3" cbytes="$4" obytes="$5"
  local cs="$6" crss="$7" ccpu="$8" ds="$9" drss="${10}" dcpu="${11}" rt="${12}"
  local cgb dgb
  cgb=$(awk -v k="$crss" 'BEGIN{printf "%.2f", k/1048576}')
  dgb=$(awk -v k="$drss" 'BEGIN{printf "%.2f", k/1048576}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$variant" "$lossless" "$cbytes" "$obytes" "$(ratio "$obytes" "$cbytes")" \
    "$cs" "$cgb" "$ccpu" "$ds" "$dgb" "$dcpu" "$rt" >> "$rf"
}

need_tool() { # path label
  if [ ! -x "$1" ] && ! command -v "$1" >/dev/null 2>&1; then
    echo "WARN: $2 not found/executable at '$1' — skipping its rows." >&2
    return 1
  fi
  return 0
}

# cluster_verify <qz> <archive> <threads> -> echoes PASS|FAIL
# Cluster output is a permutation of the input, so md5 does NOT apply; qz verify
# checks the embedded multiset checksum (set-losslessness).
cluster_verify() {
  if "$1" verify -i "$2" -t "$3" >/dev/null 2>&1; then echo PASS; else echo FAIL; fi
}

section() { echo; echo "=== $* ==="; }

show_results() { # resultfile
  echo; echo "===== RESULTS: $1 ====="; column -t -s $'\t' "$1"
}
