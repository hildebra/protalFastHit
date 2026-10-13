#!/bin/bash
# Time bench_pgz on each test file at 1, 2 and 4 threads (cores 0-3, niced), ROUNDS alternated rounds.
set -uo pipefail
SP=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
W=~/pgzdiag
cd $W
taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -I src/src -I src/src/Utilities -I src/lib \
  $SP/bench_pgz.cpp src/src/IO/ParallelGzip.cpp -lisal -lzstd -o bench_pgz 2> bench_pgz_compile.log
echo "compile exit $?"; head -5 bench_pgz_compile.log
uptime
for f in ${FILES:-one m8m m1m m64k}; do cat data/$f.gz > /dev/null; done
for round in $(seq 1 ${ROUNDS:-2}); do
  for f in ${FILES:-one m8m m1m m64k}; do
    for t in ${THREADS:-1 2 4}; do
      printf 'r%s\t%s\t' $round $f
      taskset -c 0-3 nice -n 5 ./bench_pgz data/$f.gz $t
    done
  done
done
uptime
