#!/bin/bash
# Per-pair exclusive instruction counts: the difference between a long and a short callgrind.sh run
# of the same reads, so index load, model parsing and other fixed costs cancel.
#   cg_diff.sh BIG_LABEL BIG_PAIRS SMALL_LABEL SMALL_PAIRS [N]
# Prints the top N "file:function" entries by Ir per pair, and (second table) the same summed by source file.
source "$(dirname "$0")/env.sh"
big=$PERF_DIR/cg/$1; bp=$2; small=$PERF_DIR/cg/$3; sp=$4; n=${5:-40}
flat() { [ -s $1/flat.txt ] || callgrind_annotate --auto=no --tree=none --threshold=100 $1/callgrind.out 2>/dev/null | c++filt -n > $1/flat.txt; }
flat $big; flat $small
parse() { awk '/^ *[0-9][0-9,]* +\( *[0-9.]+%\)/ && !/PROGRAM TOTALS/ { v=$1; gsub(",","",v); line=$0; sub(/^ *[0-9,]+ +\( *[0-9.]+%\) +/, "", line); sub(/ \[[^]]*\]$/, "", line); print v "\t" line }' $1/flat.txt; }
parse $big > $big/flat.tsv; parse $small > $small/flat.tsv
awk -F'\t' -v bp=$bp -v sp=$sp 'NR==FNR { s[$2] += $1; next } { d[$2] += $1 } END { for (k in d) { v = (d[k] - s[k]) / (bp - sp); if (v >= 1) printf "%.0f\t%s\n", v, k } }' $small/flat.tsv $big/flat.tsv | sort -rn > $big/perpair.tsv
tot=$(awk -F'\t' '{t += $1} END {print t}' $big/perpair.tsv)
echo "per-pair total (listed entries): $tot Ir"
head -$n $big/perpair.tsv | awk -F'\t' -v t=$tot '{ printf "%7d %5.1f%%  %s\n", $1, 100*$1/t, substr($2, 1, 190) }'
echo; echo "by source file:"
awk -F'\t' '{ f = $2; sub(/:.*/, "", f); sum[f] += $1 } END { for (f in sum) printf "%d\t%s\n", sum[f], f }' $big/perpair.tsv | sort -rn | head -25 | awk -F'\t' -v t=$tot '{ printf "%7d %5.1f%%  %s\n", $1, 100*$1/t, $2 }'
