#!/bin/bash
# Exclusive ("self") instruction counts per function of a callgrind.sh run, top N, merged over
# inlining as callgrind sees it: cg_flat.sh LABEL [N]
source "$(dirname "$0")/env.sh"
out=$PERF_DIR/cg/$1; n=${2:-40}
f=$out/flat.txt
[ -s $f ] || callgrind_annotate --auto=no --tree=none --threshold=100 $out/callgrind.out 2>/dev/null | c++filt -n > $f
grep "PROGRAM TOTALS" $f
awk -v n=$n '/^ *[0-9,]+ .*\(/ && !/PROGRAM/ { if (++c <= n) print }' $f | sed -E 's/ +/ /g' | cut -c1-200
