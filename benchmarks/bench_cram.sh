#!/usr/bin/env bash
# bz (lossless v12) vs samtools CRAM — size, time, peak RSS, CPU.
# Threads pinned equal for fairness. CRAM is reference-based (needs -T ref);
# bz is self-contained (no reference at compress or decompress).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

IN=${IN:-$HERE/../real_data/HG002.GRCh38.2x151.chr20_only.bwamem2.bam}
REF=${REF:-$HERE/../real_data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}
BZ=${BZ:-$HERE/../target/release/bz}
T=${T:-16}
D=/tmp/cram_bench; mkdir -p "$D"
BAM=$(stat -c %s "$IN")
TIME=/usr/bin/time; LOG=$D/t.log
echo "input BAM = $BAM bytes, threads = $T"

run() {  # run <label> <cmd...>
  local label="$1"; shift
  $TIME -v "$@" >/dev/null 2>"$LOG"; local rc=$?
  python3 - "$label" "$rc" "$LOG" <<'PY'
import sys,re
label,rc,log=sys.argv[1],sys.argv[2],sys.argv[3]
t=open(log).read()
g=lambda p:(re.search(p,t) or [None,"?"])[1]
rss=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)',t).group(1))
user=float(re.search(r'User time \(seconds\): ([\d.]+)',t).group(1))
sysv=float(re.search(r'System time \(seconds\): ([\d.]+)',t).group(1))
print(f"{label:26s} rc={rc} wall={g(r'Elapsed.*?: (.*)'):9s} "
      f"RSS={rss/1048576:.2f}GB CPU={g(r'Percent of CPU this job got: (.*)'):7s} cpu_s={user+sysv:.0f}")
PY
}
ratio() { # ratio <file> <name>
  python3 - "$2" "$(stat -c %s "$1")" "$BAM" <<'PY'
import sys; n,s,b=sys.argv[1],int(sys.argv[2]),int(sys.argv[3])
print(f"  {n:9s} {s:>13,} bytes  ratio {b/s:.3f}x  ({100*s/b:.1f}% of BAM)")
PY
}

echo "==== COMPRESS ===="
run "bz v12 (lossless)" "$BZ" compress -i "$IN" -o "$D/out.bz" -t $T -f
run "samtools CRAM 3.0" samtools view -C -T "$REF" -@ $T -o "$D/out30.cram" "$IN"
run "samtools CRAM 3.1" samtools view -C -T "$REF" --output-fmt-option version=3.1 -@ $T -o "$D/out31.cram" "$IN"
echo "==== ratios (vs original BAM) ===="
ratio "$D/out.bz" bz_v12; ratio "$D/out30.cram" CRAM3.0; ratio "$D/out31.cram" CRAM3.1
echo "==== DECOMPRESS ===="
run "bz decompress"         "$BZ" decompress -i "$D/out.bz" -o "$D/rt.bam" -t $T -f
run "samtools CRAM3.1->BAM" samtools view -b -T "$REF" -@ $T -o "$D/rt31.bam" "$D/out31.cram"
echo DONE
