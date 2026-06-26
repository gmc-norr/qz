#!/usr/bin/env bash
# Manual ratio + losslessness validation for single-end reference mode (archive_type 4).
# Compress HG002 chr20 R1 ALONE against chr20.fa, compare vs single-end DEFAULT,
# assert lossless roundtrip + clean verify. RUN ON AN IDLE BOX (no concurrent
# wall benchmarks) — peak-RSS / wall numbers are otherwise contaminated.
#
# Override via env: R1, REF, QZ, WORK, T.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
R1=${R1:-$HERE/../real_data/HG002.chr20.novaseq.2x151.R1.fastq}
REF=${REF:-$HERE/../real_data/chr20.fa}
QZ=${QZ:-$HERE/../target/release/qz}
WORK=${WORK:-/tmp/qz_se_ref_validate}
T=${T:-72}

rm -rf "$WORK"; mkdir -p "$WORK"
ORIG_BYTES=$(stat -c %s "$R1")
echo "R1=$R1  ORIG_BYTES=$ORIG_BYTES  REF=$REF  threads=$T  $(date)"

# 1. single-end reference
/usr/bin/time -v "$QZ" compress -i "$R1" -o "$WORK/ref.qz" --reference "$REF" -w "$WORK" -t "$T" -f \
  2> "$WORK/ref.compress.time"
REF_BYTES=$(stat -c %s "$WORK/ref.qz")
AT=$(od -An -tu1 -N4 "$WORK/ref.qz" | awk '{print $4}')
echo "archive_type byte = $AT (expect 4)"

# 2. single-end default (no --reference)
/usr/bin/time -v "$QZ" compress -i "$R1" -o "$WORK/def.qz" -w "$WORK" -t "$T" -f \
  2> "$WORK/def.compress.time"
DEF_BYTES=$(stat -c %s "$WORK/def.qz")

# 3. lossless roundtrip of the reference archive
/usr/bin/time -v "$QZ" decompress -i "$WORK/ref.qz" -o "$WORK/ref.out.fastq" -w "$WORK" -t "$T" -f \
  2> "$WORK/ref.decompress.time"
if cmp -s "$R1" "$WORK/ref.out.fastq"; then RT=OK; else RT=MISMATCH; fi

# 4. verify deep + fast
"$QZ" verify -i "$WORK/ref.qz" -t "$T"        && VD=OK || VD=FAIL
"$QZ" verify -i "$WORK/ref.qz" -t "$T" --fast && VF=OK || VF=FAIL

ref_ratio=$(echo "scale=3; $ORIG_BYTES/$REF_BYTES" | bc)
def_ratio=$(echo "scale=3; $ORIG_BYTES/$DEF_BYTES" | bc)
lift=$(echo "scale=3; $DEF_BYTES/$REF_BYTES" | bc)
c_rss=$(awk '/Maximum resident/{print $6}' "$WORK/ref.compress.time")
d_rss=$(awk '/Maximum resident/{print $6}' "$WORK/ref.decompress.time")

echo "=== single-end reference validation ==="
echo "orig_bytes      : $ORIG_BYTES"
echo "reference_bytes : $REF_BYTES  (ratio ${ref_ratio}x)"
echo "default_bytes   : $DEF_BYTES  (ratio ${def_ratio}x)"
echo "lift (def/ref)  : ${lift}x  (>1 => reference smaller, expected PASS)"
echo "compress RSS KB : $c_rss"
echo "decompress RSS KB: $d_rss"
echo "roundtrip       : $RT  (expect OK)"
echo "verify deep     : $VD  (expect OK)"
echo "verify fast     : $VF  (expect OK)"

[ "$AT" = "4" ] && [ "$RT" = "OK" ] && [ "$VD" = "OK" ] && [ "$VF" = "OK" ] \
  && echo "VALIDATION: PASS" || { echo "VALIDATION: FAIL"; exit 1; }
