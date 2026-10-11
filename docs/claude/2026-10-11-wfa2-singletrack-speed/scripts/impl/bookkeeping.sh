#!/usr/bin/env bash
# bookkeeping.sh: self instructions per function of the new build's replay (callgrind of replay.sh), and per source line
# of the per-step bookkeeping (wavefront_compute.c, wavefront_compute_affine.c's driver, wavefront_slab.c, wavefront.c),
# lines with at least 0.3% of wavefront_align's instructions.
set -uo pipefail
W=$HOME/wfaspeed; R=$W/replay; H=$W/new/lib/wfa2-lib
for s in pe pb ont; do
  for t in singletrack classic; do
    o=$R/cg.$s.$t.new.out
    tot=$(awk -v s=$s -v t=$t '$1==s && $2==t && $3=="new" {print $4}' $R/cg.tsv)
    echo "== $s $t: wavefront_align $tot instructions; self shares"
    callgrind_annotate --threshold=100 $o 2>/dev/null | sed -nE 's/^ *([0-9,]+) +\(.*\) +[^ ]*:([a-zA-Z_0-9.]+) .*/\1 \2/p; s/^ *([0-9,]+) +[^ ]*:([a-zA-Z_0-9.]+) .*/\1 \2/p' | tr -d , | \
      awk -v T=$tot '$1/T >= 0.01 {printf "   %5.1f%%  %s\n", 100*$1/T, $2}' | head -24
    if [ $t = singletrack ]; then
      callgrind_annotate --auto=yes --show-percs=no --context=0 --include=$H/wavefront $o 2>/dev/null > $R/lines.$s.$t.txt
      echo "   -- lines (>= 0.3%) of the bookkeeping"
      awk -v T=$tot '
        /^-- Auto-annotated source:/ {f=$NF; sub(".*/","",f); next}
        f ~ /^(wavefront_compute\.c|wavefront_compute_affine\.c|wavefront_slab\.c|wavefront\.c)$/ {
          if (match($0, /^ *[0-9,]+ /)) { n=substr($0, RSTART, RLENGTH); gsub(/[ ,]/,"",n); rest=substr($0, RLENGTH+1);
            if (n/T >= 0.003) printf "   %5.2f%%  %-26s %s\n", 100*n/T, f, substr(rest,1,110) } }' $R/lines.$s.$t.txt | head -60
    fi
  done
done
