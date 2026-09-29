#!/bin/bash
# Read throughput of SeqReaderPE alone (threads take pairs and do nothing else), igzstream against
# ThreadedGzIstream: bench_reader.sh R1.fq.gz R2.fq.gz [REPS]
source "$(dirname "$0")/env.sh"
r1=$1; r2=$2; reps=${3:-3}
g++ -std=c++20 -O3 -fopenmp -I$PROTAL_SRC/src -I$PROTAL_SRC/src/IO -I$PROTAL_SRC/lib -I$PROTAL_SRC/lib/gzstream \
  $here/bench_reader.cpp $PROTAL_SRC/src/IO/FastxReader.cpp $PROTAL_SRC/lib/gzstream/gzstream.C -lz -o $PERF_DIR/bench_reader || exit 1
for rep in $(seq 1 $reps); do
  for t in 1 8; do
    $PERF_DIR/bench_reader gz $r1 $r2 $t
    $PERF_DIR/bench_reader threaded $r1 $r2 $t
  done
done
