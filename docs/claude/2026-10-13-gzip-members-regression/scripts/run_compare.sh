#!/bin/bash
# Before (9c27756) against after (the working tree with the member reader and the released pages): bench_pgz built
# against each, on the files of bench_data.sh, at 1, 2 and 4 gzip threads, ROUNDS rounds alternating old and new;
# cores 0-3, niced. old.tar: `git archive -o scripts/old.tar 9c27756 src lib/gzstream` from Git Bash (not kept); the
# new sources are ~/pgzfix/src (the working tree synced there).
set -uo pipefail
SP=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
B=~/pgzbench
mkdir -p $B/old/src
tar -xf $SP/old.tar -C $B/old/src
CXX="taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG"
$CXX -I $B/old/src/src -I $B/old/src/src/Utilities -I $B/old/src/lib $SP/bench_pgz.cpp $B/old/src/src/IO/ParallelGzip.cpp \
  -lisal -lzstd -o $B/bench_old || exit 1
$CXX -DPGZ_NEW -I ~/pgzfix/src/src -I ~/pgzfix/src/src/Utilities -I ~/pgzfix/src/lib $SP/bench_pgz.cpp \
  ~/pgzfix/src/src/IO/ParallelGzip.cpp -lisal -lzstd -o $B/bench_new || exit 1
cd $B/data
for f in ${FILES:-one m8m m1m m64k mixed}; do cat $f.gz > /dev/null; done
uptime
for round in $(seq 1 ${ROUNDS:-3}); do
  for f in ${FILES:-one m8m m1m m64k mixed}; do
    for t in ${THREADS:-1 2 4}; do
      for v in old new; do
        printf 'r%s\t%s\t%s\t' $round $f $v
        taskset -c 0-3 nice -n 5 $B/bench_$v $f.gz $t
      done
    done
  done
done
uptime
