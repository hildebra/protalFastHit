#!/bin/bash
# rapidgzip directly: baseline vs -march=native, chunk and read sizes, with its statistics once.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/rgzbuild
B=$W/build
ZNG=$B/zlib-ng
cd $W
for arch in x86-64 native; do
  g++ -O3 -std=c++20 -march=$arch -pthread -DNDEBUG -DLIBRAPIDARCHIVE_WITH_ISAL \
    -isystem src/lib/rapidgzip/src -isystem src/lib/rapidgzip-zlib -isystem $B/isa-l-install/include/isa-l -I $ZNG \
    $SP/bench_rgz.cpp $B/isa-l-install/lib/libisal.a $(find $B -name "libz-ng.a" | head -1) -o bench_rgz_$arch 2> bench_rgz_$arch.log
  echo "compile $arch exit $?"
done
F=~/gzpar/data/R1.member.gz
taskset -c 0-3 ./bench_rgz_x86-64 $F 3 4 4 1 2>&1 | tail -60
for round in 1 2; do
  for arch in x86-64 native; do
    for cfg in "1 4 4" "3 4 1" "3 4 4" "3 4 16" "3 8 8" "3 1 4" "4 4 4"; do
      echo -n "$arch: "; taskset -c 0-3 ./bench_rgz_$arch $F $cfg 0
    done
  done
done
