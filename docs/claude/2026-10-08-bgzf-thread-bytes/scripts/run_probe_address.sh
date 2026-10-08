#!/bin/bash
# Builds and runs the address probe on block 6 (content 294045, 3580 bytes; BGZF offset 141815) of the failing s2.
set -eu
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/99a2c9fe-d722-41e1-bc32-3fc6baabdeb2/scratchpad
cd ~/lrdet/fail39
cp $SP/probe_address.cpp .
taskset -c 0-3 g++ -O2 -g -std=c++17 probe_address.cpp -o probe_address -lisal
taskset -c 0-3 nice -n 5 ./probe_address one_s2.fq 294045 3580 one_s2.fq.gz four_s2.fq.gz 141815 8192
