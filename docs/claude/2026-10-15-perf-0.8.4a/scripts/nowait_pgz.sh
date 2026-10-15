#!/usr/bin/env bash
# Trial A5: the lead decoding a piece a worker is still on (PGZ_NOWAIT) against waiting for it (protal's), with the
# lead's wait times (PGZ_STATS); the experiment copy of ParallelGzip.cpp, 3 alternated rounds, cores 0-3.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/95e4f938-fd80-46a7-b960-30fecd6a2871/scratchpad
W=~/perf085
cd $W
mkdir -p exp
cp $SP/exp/ParallelGzip.cpp exp/ParallelGzip.cpp
taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -DPGZ_NEW -I src/src -I src/src/IO -I src/src/Utilities -I src/lib \
  bench_pgz.cpp exp/ParallelGzip.cpp -lisal -lzstd -o bin/bench_pgz_exp 2> exp_compile.log || { cat exp_compile.log; exit 1; }
cd data
OUT=$W/pgz_nowait.txt
: > $OUT
for round in 1 2 3; do
  for f in ont.g6.gz hifi.ena.gz; do
    for t in 2 3 4; do
      for mode in wait nowait; do
        if [ $mode = nowait ]; then export PGZ_NOWAIT=1; else unset PGZ_NOWAIT; fi
        res=$(PGZ_STATS=1 taskset -c 0-3 $W/bin/bench_pgz_exp $f $t 2>&1)
        printf 'round %s\tfile %s\tmode %s\t%s\t%s\n' $round $f $mode "$(echo "$res" | grep threads | tail -1)" "$(echo "$res" | grep pgz_stats | tail -1)" >> $OUT
      done
    done
  done
done
unset PGZ_NOWAIT
awk -F'\t' '{for(i=1;i<=NF;i++){n=split($i,kv," "); for(j=1;j+1<=n;j+=2) v[kv[j]]=kv[j+1]} key=v["file"]" t"v["threads"]" "v["mode"]; g[key]=g[key]" "v["GB/s"]; c[key]="waits="v["waits"]" wait_s="v["wait_s"]" not_waited="v["not_waited"]" in_order="v["in_order"]" in_order_s="v["in_order_s"]" resolve_wait_s="v["resolve_wait_s"]" cpu="v["cpu"]} END{for(k in g) print k"\tGB/s"g[k]"\t"c[k]}' $OUT | sort
