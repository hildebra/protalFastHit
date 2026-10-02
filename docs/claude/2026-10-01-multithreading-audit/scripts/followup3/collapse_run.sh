#!/usr/bin/env bash
# FALLOC_FL_COLLAPSE_RANGE on a freshly written file of 100 KB, 40 MB and 400 MB (dirty, or synced first), three
# passes each, on ext4 in the home folder. After the real-run check.
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
W=$HOME/mt-work/copy; mkdir -p $W
g++ -std=c++20 -O2 $SP/collapse.cpp -o $W/collapse || exit 1
cat /proc/loadavg
for rep in 1 2 3; do
  for b in 100000 40000000 400000000; do
    $W/collapse $W/c.bin $b
    $W/collapse $W/c.bin $b sync
  done
done
echo COLLAPSE DONE
