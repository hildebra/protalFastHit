#!/usr/bin/env bash
# Builds the two experiment binaries of this report in ~/mt-audit (WSL, Release, Ninja, niced):
#   bin/protal-instr  1c11a00 with mtaudit.patch (wall-clock counters around the locks and hand-offs,
#                     printed as MTAUDIT lines on stderr when protal exits)
#   readbench/readbench  the input path alone (readbench.cpp), against the patched tree
# and copies the plain 0.7.1 build of the benchmark of the same day (1c11a00, the same flags) as
# bin/protal-base. Neither the patch nor readbench is meant for the tree.
set -euo pipefail
M=${MT_AUDIT:-$HOME/mt-audit}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${REPO:-/mnt/c/Users/hildebra/Documents/locDev/protal}
mkdir -p $M/bin $M/logs $M/readbench
if [ ! -d $M/instr-src ]; then
  mkdir -p $M/instr-src
  git -C $REPO archive 1c11a00 | tar -x -C $M/instr-src
  (cd $M/instr-src && patch -p1 -s < $HERE/mtaudit.patch)
fi
cd $M/instr-src
[ -f build/build.ninja ] || nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $M/logs/configure_instr.log 2>&1
nice cmake --build build --target protal -j 4 > $M/logs/build_instr.log 2>&1 || { tail -30 $M/logs/build_instr.log; exit 1; }
cp build/protal $M/bin/protal-instr
cp $HOME/bench071/bin/protal-0.7.1 $M/bin/protal-base

T=$M/instr-src
cp $HERE/readbench.cpp $M/readbench/
cd $M/readbench
g++ -std=c++20 -O3 -march=x86-64-v2 -fopenmp -I$T/src -I$T/src/IO -I$T/src/SequenceUtils -I$T/src/Utilities \
    -I$T/build/zlib-ng -I$T/lib/zlib-ng readbench.cpp $T/src/IO/FastxReader.cpp \
    $T/build/zlib-ng/libz-ng.a -ldeflate -lpthread -o readbench
ls -la $M/bin $M/readbench/readbench
