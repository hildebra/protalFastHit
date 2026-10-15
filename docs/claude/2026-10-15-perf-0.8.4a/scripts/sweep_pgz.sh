#!/usr/bin/env bash
# Trial A2: ParallelGzip's piece size and gap (the lead's in-order lead before the workers' first piece) on the
# long-read files, 3 and 4 threads, 2 alternated rounds on cores 0-3.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/95e4f938-fd80-46a7-b960-30fecd6a2871/scratchpad
W=~/perf085
cd $W
cp $SP/bench_pgz.cpp $W/
taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -DPGZ_NEW -I src/src -I src/src/Utilities -I src/lib \
  bench_pgz.cpp src/src/IO/ParallelGzip.cpp -lisal -lzstd -o bin/bench_pgz 2> bench_compile.log || { cat bench_compile.log; exit 1; }
cd data
OUT=$W/pgz_sweep.txt
: > $OUT
cat ont.g6.gz hifi.ena.gz > /dev/null
for round in 1 2; do
  for f in ont.g6.gz hifi.ena.gz; do
    for t in 3 4; do
      for piece in 524288 1048576 2097152; do
        for gap in 1 2 3 0; do
          line=$(PGZ_PIECE=$piece PGZ_GAP=$gap taskset -c 0-3 $W/bin/bench_pgz $f $t 2>&1 | tail -1)
          printf 'round %s\tfile %s\tpiece %s\tgap %s\t%s\n' $round $f $piece $gap "$line" >> $OUT
        done
      done
    done
  done
done
awk -F'\t' '{for(i=1;i<=NF;i++){split($i,kv," "); v[kv[1]]=kv[2]} key=v["file"]" t"v["threads"]" piece"v["piece"]" gap"v["gap"]; g[key]=g[key]" "v["GB/s"]; c[key]="in_order="v["in_order"]" ahead="v["ahead"]" cpu="v["cpu"]" rss="v["maxrss_mb"]} END{for(k in g) print k"\tGB/s"g[k]"\t"c[k]}' $OUT | sort
