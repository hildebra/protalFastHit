#!/bin/bash
# The syncmer scan alone (bench_syncmers.cpp), HEAD's against the tree's, on the reads of a FASTQ, pinned and alternated:
#   [PIN=cpu] bench_syncmers.sh TREE_A TREE_B READS.fq [READS] [ROUNDS] [REPS]   (trees: $PERF_DIR/pair/<name>/src)
# Each tree is built twice: as it is, and stopping after the per-window codes are filled (FILL_ONLY: the k-mers, cores
# and s-mers, before the AVX2 evaluation of the windows), so the fill and the evaluation can be told apart.
source "$(dirname "$0")/env.sh"
A=$1; B=$2; fq=$3; reads=${4:-200000}; rounds=${5:-15}; reps=${6:-3}
W=$PERF_DIR/bsync; mkdir -p $W
for t in $A $B; do
  T=$PERF_DIR/pair/$t/src
  inc=$(find $T/src -type d | sed "s/^/-I/" | tr "\n" " ")
  for v in full fill; do
    mkdir -p $W/$t-$v/SequenceUtils
    cp $T/src/SequenceUtils/KmerIterator.h $W/$t-$v/SequenceUtils/
    if [ $v = fill ]; then
      # stop where the window evaluation starts, in the AVX2 scan only (the first match is ScanWindows')
      awk '/void ScanWindowsAvx2/ { in_scan = 1 } in_scan && /uint32_t const count = m_minimizer.SmerCount\(\), t = m_minimizer.T\(\);/ { print "            return;  // FILL_ONLY"; in_scan = 0 } { print }' \
        $T/src/SequenceUtils/KmerIterator.h > $W/$t-$v/SequenceUtils/KmerIterator.h
      grep -q "FILL_ONLY" $W/$t-$v/SequenceUtils/KmerIterator.h || { echo "no stop inserted in $t"; exit 1; }
    fi
    g++ -O3 -DNDEBUG -std=c++20 -march=x86-64-v3 -mtune=generic -pthread -I$W/$t-$v $inc "$(dirname "$0")/bench_syncmers.cpp" -o $W/bench_$t-$v || exit 1
  done
done
for i in $(seq 1 $reps); do
  for v in full fill; do
    for t in $A $B; do echo "$t $v: $(${PIN:+taskset -c $PIN} $W/bench_$t-$v $fq $reads $rounds)"; done
  done
done
