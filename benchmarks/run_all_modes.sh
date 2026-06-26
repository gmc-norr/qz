#!/usr/bin/env bash
# Run the end-to-end mode benchmarks (single, paired, reference, bam) and print
# a combined report. Each sub-script writes results/<mode>.tsv + results/<mode>.log.
#
# Usage:
#   ./run_all_modes.sh                 # all four modes
#   ./run_all_modes.sh single bam      # only the named modes
#   T=36 ./run_all_modes.sh            # override thread count for all
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
MODES=("$@")
[ "${#MODES[@]}" -eq 0 ] && MODES=(single paired reference bam)

declare -A SCRIPT=(
  [single]="$HERE/bench_single_end.sh"
  [paired]="$HERE/bench_paired_end.sh"
  [reference]="$HERE/bench_reference.sh"
  [bam]="$HERE/bench_bam.sh"
)

echo "############################################################"
echo "# qz/bz mode benchmarks — $(date)"
echo "# threads=${T:-72}  modes=${MODES[*]}"
echo "############################################################"

for m in "${MODES[@]}"; do
  s="${SCRIPT[$m]:-}"
  if [ -z "$s" ]; then echo "unknown mode: $m (valid: single paired reference bam)"; continue; fi
  echo; echo "################## MODE: $m ##################"
  bash "$s" || echo "WARN: mode '$m' exited non-zero (partial results may exist)"
done

echo; echo "############################################################"
echo "# COMBINED RESULTS"
echo "############################################################"
for m in "${MODES[@]}"; do
  tsv="$HERE/results/${m}.tsv"
  [ -f "$tsv" ] || continue
  echo; echo "## $m"
  column -t -s $'\t' "$tsv"
done
echo
echo "TSVs + per-tool logs in: $HERE/results/"
