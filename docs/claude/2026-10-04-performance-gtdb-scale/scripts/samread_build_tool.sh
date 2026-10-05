#!/usr/bin/env bash
set -euo pipefail
W=$HOME/perf-gtdb; X=$W/samread; S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
g++ -O2 -std=c++20 -I$W/work/src -I$W/work/include $S/samread/seekable.cpp -o $X/seekable -lzstd -lpthread 2>&1 | head -20
for f in mapped mixed; do $X/seekable $X/$f.sam $X/$f.sam.zst; done
ls -la $X
