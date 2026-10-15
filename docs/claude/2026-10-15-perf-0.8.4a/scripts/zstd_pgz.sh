#!/usr/bin/env bash
# Trial A3: the ONT FASTQ as zstd (one frame from zstd -T4; frames of 4 MB from split + zstd, as a multi-frame file)
# read through ThreadedGzIstream on 1-4 threads, against gzip -6 (one member).
set -uo pipefail
W=~/perf085
cd $W/data
command -v pzstd || true
[ -s ont.z3.zst ] || taskset -c 0-3 nice -n 5 zstd -q -3 -T4 -c ont.fq > ont.z3.zst
[ -s ont.zf.zst ] || taskset -c 0-3 nice -n 5 split -b 4M --filter='zstd -q -3 -c' ont.fq > ont.zf.zst
ls -la ont.g6.gz ont.z3.zst ont.zf.zst
cat ont.g6.gz ont.z3.zst ont.zf.zst > /dev/null
for round in 1 2; do
  for f in ont.g6.gz ont.z3.zst ont.zf.zst; do
    for t in 1 2 4; do
      line=$(taskset -c 0-3 $W/bin/bench_pgz $f $t 2>&1 | tail -1)
      printf 'round %s\tfile %s\t%s\n' $round $f "$line"
    done
  done
done | awk -F'\t' '{for(i=1;i<=NF;i++){split($i,kv," "); v[kv[1]]=kv[2]} print v["file"]" t"v["threads"]" GB/s "v["GB/s"]" cpu "v["cpu"]" rss "v["maxrss_mb"]}'
