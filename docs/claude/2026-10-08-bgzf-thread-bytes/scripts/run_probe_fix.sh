#!/bin/bash
# Builds the probe of the fix against the working tree's Bgzf.h (read from /mnt/c, a header only) and runs it.
set -eu
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/99a2c9fe-d722-41e1-bc32-3fc6baabdeb2/scratchpad
W=~/lrdet/fix
mkdir -p $W/IO
cd $W
cp /mnt/c/Users/hildebra/Documents/locDev/protal/src/IO/Bgzf.h IO/
cp $SP/probe_fix.cpp .
taskset -c 0-3 nice -n 5 g++ -O2 -g -std=c++20 -Wall -Wextra -msse4.2 -I . probe_fix.cpp -o probe_fix -lisal
for run in 1 2 3; do
    echo "== process $run"
    taskset -c 0-3 nice -n 5 ./probe_fix ~/lrdet/fail39/one_s2.fq 294045 3580
done
