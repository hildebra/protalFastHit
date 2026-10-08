#!/bin/bash
set -eu
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/99a2c9fe-d722-41e1-bc32-3fc6baabdeb2/scratchpad
W=~/lrdet/fix
cd $W
cp /mnt/c/Users/hildebra/Documents/locDev/protal/src/IO/Bgzf.h IO/
cp $SP/probe_buckets.cpp .
taskset -c 0-3 nice -n 5 g++ -O2 -g -std=c++20 -msse4.2 -I . probe_buckets.cpp -o probe_buckets -lisal
taskset -c 0-3 nice -n 5 ./probe_buckets ~/lrdet/fail39/one_s2.fq 294045 3580
uname -r
