#!/bin/bash
# Simulated last-level data-read misses and first-level read misses per pair, by function: the difference between
# a long and a short cachesim.sh run.  cg_cachediff.sh BIG_LABEL BIG_PAIRS SMALL_LABEL SMALL_PAIRS [N]
source "$(dirname "$0")/env.sh"
big=$PERF_DIR/cg/$1; bp=$2; small=$PERF_DIR/cg/$3; sp=$4; n=${5:-20}
parse() {
  [ -s $1/flatc.txt ] || callgrind_annotate --auto=no --tree=none --threshold=100 --show=Ir,D1mr,DLmr $1/callgrind.out 2>/dev/null | c++filt -n > $1/flatc.txt
  awk '/^ *[0-9][0-9,]* +\( *[0-9.]+%\)/ && !/PROGRAM TOTALS/ {
         line = $0; gsub(/,/, "", line)
         match(line, /^ *[0-9]+ +\( *[0-9.]+%\) +[0-9]+ +\( *[0-9.]+%\) +[0-9]+ +\( *[0-9.]+%\) +/)
         head = substr(line, 1, RLENGTH); name = substr(line, RLENGTH + 1); sub(/ \[[^]]*\]$/, "", name)
         split(head, f, /[ ()%]+/); ir = f[2]; d1 = f[4]; dl = f[6]
         print name "\t" ir "\t" d1 "\t" dl }' $1/flatc.txt
}
parse $big > $big/c.tsv; parse $small > $small/c.tsv
echo "per pair: Ir, D1 read misses, LL read misses (simulated 12 MB LL)"
awk -F'\t' -v bp=$bp -v sp=$sp -v n=$n -v tmp=/tmp/cachediff.$$ 'NR==FNR { ir[$1] = $2; d1[$1] = $3; dl[$1] = $4; next }
  { i = ($2 - ir[$1]) / (bp - sp); d = ($3 - d1[$1]) / (bp - sp); l = ($4 - dl[$1]) / (bp - sp); ti += i; td += d; tl += l
    if (l >= 0.3) printf "%8.0f %7.1f %7.1f  %s\n", i, d, l, substr($1, 1, 150) > tmp
  } END { printf "%8.0f %7.1f %7.1f  TOTAL\n", ti, td, tl }' $small/c.tsv $big/c.tsv
sort -k3 -rn /tmp/cachediff.$$ | head -$n; rm -f /tmp/cachediff.$$
