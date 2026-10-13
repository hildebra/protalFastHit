#!/bin/bash
# Sync the working tree, compile bench_gz2 (bench_gz + PGZ_PIECE) with the current ParallelGzip, time THREADS threads
# at PIECES piece sizes (0: the default), with PROTAL_PGZ_TIMING, ROUNDS rounds.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
P=/mnt/c/Users/hildebra/Documents/locDev/protal
W=~/pgzbuild
rsync -a --checksum --no-times --delete \
  --exclude='/data' --exclude='/build*/' --exclude='/local' --exclude='/.git' --exclude='/.claude' --exclude='/misc' \
  --exclude='/tests/misc' --exclude='/docs/claude' $P/ $W/src/
cd $W
g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -I src/src -I src/src/Utilities -I src/lib \
  $SP/bench_gz2.cpp src/src/IO/ParallelGzip.cpp -lisal -lzstd -o bench_gz 2> bench_compile.log
echo "compile exit $?"; head -5 bench_compile.log
F=~/gzpar/data/R1.member.gz
cat $F > /dev/null
for round in $(seq 1 ${ROUNDS:-2}); do
  for t in ${THREADS:-1 2 3 4}; do
    for piece in ${PIECES:-0}; do
      PGZ_PIECE=$piece PROTAL_PGZ_TIMING=1 taskset -c 0-3 ./bench_gz $F $t 2>&1 | sed "s/^/r$round p$piece\t/" | cut -c1-280
    done
  done
done
