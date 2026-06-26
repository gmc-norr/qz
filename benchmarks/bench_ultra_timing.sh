#!/usr/bin/env bash
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SM=$HERE/../real_data/ERR3239334_1.10m.fastq
W=$HERE/../real_data/_vt; mkdir -p "$W"
QZ=$HERE/../target/release/qz
OMD5=$(md5sum "$SM"|awk '{print $1}')
RES=$HERE/results_ultra.tsv
echo -e "mode\tcomp_bytes\tcomp_s\tcomp_rss_kb\tdecomp_s\tdecomp_rss_kb\tmd5" > "$RES"
gett(){ grep "Elapsed" "$1"|sed -E 's/.*: //'|awk -F: '{if(NF==3)print $1*3600+$2*60+$3;else print $1*60+$2}'; }
getr(){ grep "Maximum resident" "$1"|sed -E 's/.*: //'; }
for L in 1 2 3 4 5; do
  /usr/bin/time -v $QZ compress -i "$SM" -o "$W/u.qz" -w "$W" -t 72 -f --ultra $L >/dev/null 2>"$W/ct"
  cs=$(gett "$W/ct"); crss=$(getr "$W/ct"); cb=$(stat -c %s "$W/u.qz")
  /usr/bin/time -v $QZ decompress -i "$W/u.qz" -o "$W/u.out" -w "$W" -t 72 -f >/dev/null 2>"$W/dt"
  ds=$(gett "$W/dt"); drss=$(getr "$W/dt")
  [ "$(md5sum "$W/u.out"|awk '{print $1}')" = "$OMD5" ] && m=OK || m=MISMATCH
  echo -e "qz_ultra$L\t$cb\t$cs\t$crss\t$ds\t$drss\t$m" >> "$RES"
  rm -f "$W/u.qz" "$W/u.out"
done
rm -rf "$W"; echo DONE >> "$RES"
