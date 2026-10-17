#!/usr/bin/env bash
# Sorted SAM records of ref and work (last round), and the tally's cost per GB of SAM text (tally_bench.cpp, -O3,
# x86-64-v3 as protal's hot code), on core 0.
ROOT=$HOME/spikecount
here=$(dirname "$0")
runs=$ROOT/runs
for s in dohlman hku; do
    a=$(ls "$runs/ref.$s.3/$s".sam* | head -1); b=$(ls "$runs/work.$s.3/$s".sam* | head -1)
    if cmp -s <(zstd -dcq "$a" | grep -v '^@' | sort) <(zstd -dcq "$b" | grep -v '^@' | sort); then
        echo "$s: SAM records identical (sorted), $(zstd -dcq "$b" | grep -vc '^@') records"
    else
        echo "$s: SAM records DIFFER (sorted)"
    fi
done
g++ -std=c++20 -O3 -march=x86-64-v3 -I "$ROOT/work/src" -o "$ROOT/tally_bench" "$here/tally_bench.cpp" || exit 1
for s in hku dohlman; do
    b=$(ls "$runs/work.$s.3/$s".sam* | head -1)
    zstd -dcq "$b" | taskset -c 0 "$ROOT/tally_bench"
done
