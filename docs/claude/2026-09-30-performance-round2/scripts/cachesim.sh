#!/bin/bash
# Simulated cache misses (callgrind --cache-sim) of the first PAIRS pairs, 1 thread: cachesim.sh LABEL DB READS_DIR PAIRS
# The last-level cache is simulated as 12 MB, 12-way, 64-byte lines (this CPU's L3). Simulated only: it has no
# prefetcher, no TLB and no timing; it says how many distinct lines a stage touches that a 12 MB cache cannot hold.
source "$(dirname "$0")/env.sh"
label=$1; db=$2; reads=$3; pairs=$4
out=$PERF_DIR/cg/$label; rm -rf $out; mkdir -p $out/reads
head -n $((pairs*4)) $(ls $reads/*_R1.fq* | head -1) > $out/reads/sub_R1.fq
head -n $((pairs*4)) $(ls $reads/*_R2.fq* | head -1) > $out/reads/sub_R2.fq
valgrind --tool=callgrind --cache-sim=yes --LL=12582912,12,64 --callgrind-out-file=$out/callgrind.out \
  $PROF_BIN --db $db -1 $out/reads/sub_R1.fq -2 $out/reads/sub_R2.fq -o $out/out -t 1 --no_qcmsa --no_profile > $out/stdout.log 2> $out/stderr.log
tail -1 $out/stderr.log
