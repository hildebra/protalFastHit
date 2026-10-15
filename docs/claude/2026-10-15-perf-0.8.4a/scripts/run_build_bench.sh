#!/usr/bin/env bash
# The hashed key map's build at r226's size: 4d89d86's HashedKeys.h (one shared counter) against the working tree's
# (per-partition tallies), on 1, 2 and 4 threads (cores 0-3).
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/95e4f938-fd80-46a7-b960-30fecd6a2871/scratchpad/hkeys
W=~/hkbuild
mkdir -p $W/new/Hash $W/old/Hash $W/common
cd $W
tar -xf $SP/src_now.tar -C $W/common
cp $W/common/src/Hash/HashedKeys.h $W/new/Hash/HashedKeys.h
cp $SP/HashedKeys_4d89d86.h $W/old/Hash/HashedKeys.h
for v in old new; do
  taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -fopenmp -DNDEBUG -I $W/$v -I $W/common/src -I $W/common/src/Hash \
    -I $W/common/src/Utilities $SP/hk_build_bench.cpp -o $W/bench_$v 2> $W/compile_$v.log || { cat $W/compile_$v.log; exit 1; }
done
free -g | head -2
for v in old new; do
  echo "== $v"
  taskset -c 0-3 $W/bench_$v 1 2 4
done
