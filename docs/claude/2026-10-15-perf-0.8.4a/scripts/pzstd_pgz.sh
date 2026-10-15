#!/usr/bin/env bash
# Trial A4: the ONT FASTQ by pzstd (independent frames with their sizes), read on 1-4 threads, against gzip -6.
set -uo pipefail
W=~/perf085
cd $W/data
[ -s ont.pz.zst ] || taskset -c 0-3 nice -n 5 pzstd -q -3 -p 4 -c ont.fq > ont.pz.zst
ls -la ont.pz.zst
cat ont.g6.gz ont.pz.zst > /dev/null
for round in 1 2 3; do
  for f in ont.g6.gz ont.pz.zst; do
    for t in 1 2 4; do
      line=$(taskset -c 0-3 $W/bin/bench_pgz $f $t 2>&1 | tail -1)
      printf 'round %s\tfile %s\t%s\n' $round $f "$line"
    done
  done
done | awk -F'\t' '{for(i=1;i<=NF;i++){split($i,kv," "); v[kv[1]]=kv[2]} print v["file"]" t"v["threads"]" GB/s "v["GB/s"]" cpu "v["cpu"]" rss "v["maxrss_mb"]}' | sort | awk '{k=$1" "$2; g[k]=g[k]" "$4} END{for(k in g) print k": GB/s"g[k]}' | sort
