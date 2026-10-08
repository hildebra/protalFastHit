#!/bin/bash
# Builds and runs the deflate probe on block 6 (content 294045, 3580 bytes) of the failing run's s2.
set -eu
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/99a2c9fe-d722-41e1-bc32-3fc6baabdeb2/scratchpad
W=~/lrdet/fail39
cd $W
cp $SP/probe_deflate.cpp .
taskset -c 0-3 g++ -O2 -g -std=c++17 -I ~/pipefix/src/src probe_deflate.cpp -o probe_deflate -lisal
echo "one-thread block: $((1905 - 26)) bytes of deflate; four-thread block: $((1901 - 26))"
taskset -c 0-3 ./probe_deflate one_s2.fq 294045 3580
echo "--- valgrind (fresh, exact-size buffer only matters)"
taskset -c 0-3 valgrind --error-exitcode=9 -q ./probe_deflate one_s2.fq 294045 3580 2>&1 | head -60 || true
