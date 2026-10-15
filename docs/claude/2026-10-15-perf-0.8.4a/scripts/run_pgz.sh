#!/usr/bin/env bash
# Trial A: ParallelGzip on single-member long-read gzip files, 1-4 threads, 3 alternated rounds on cores 0-3.
# Usage: run_pgz.sh BIN OUT [files...]
set -uo pipefail
BIN=$1; OUT=$2; shift 2
cd ~/perf085/data
FILES=${*:-ont.g6.gz ont.g1.gz ont.p6.gz ont.ena.gz hifi.ena.gz}
cat $FILES > /dev/null   # into the page cache
: > $OUT
for round in 1 2 3; do
  for f in $FILES; do
    for t in 1 2 3 4; do
      line=$(taskset -c 0-3 $BIN $f $t 2>&1 | tail -1)
      printf 'round %s\tfile %s\t%s\n' $round $f "$line" >> $OUT
    done
  done
done
awk -F'\t' '{split($2,f," "); split($3,t," "); for(i=1;i<=NF;i++){split($i,kv," "); v[kv[1]]=kv[2]} key=f[2]" t"v["threads"]; g[key]=g[key]" "v["GB/s"]; c[key]="ahead="v["ahead"]" in_order="v["in_order"]" streamed="v["streamed"]" cpu="v["cpu"]} END{for(k in g) print k"\tGB/s"g[k]"\t"c[k]}' $OUT | sort
