#!/usr/bin/env bash
# The profiler's reader path alone (readsam: SamInput + the chunk reader, nothing done with the chunks), built
# against the committed src (final2 = f359d49 + 85e53e6 + a212559), on the 5M-pair SAM, with 0-4 decompression
# threads, three passes, niced. After the callgrind runs. Rows into ~/mt-work/input/readsam.tsv.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
T=$HOME/mt-work/final2/src; W=$HOME/mt-work/input
mkdir -p $W/readsam
g++ -std=c++20 -O3 -march=x86-64-v2 -I$T/src -I$T/src/IO -I$T/src/Utilities -I$T/src/Profiling -I$T/build/zlib-ng -I$T/lib/zlib-ng \
    $SP/readsam.cpp $T/build/zlib-ng/libz-ng.a -lzstd -ldeflate -lpthread -o $W/readsam/readsam || { echo "FAIL build"; exit 1; }
SAM=$HOME/mt-audit/runs/b6/s1.sam.zst
cat $SAM > /dev/null
[ -f $W/readsam.tsv ] || printf "rep\tload\tline\n" > $W/readsam.tsv
for rep in 1 2 3; do for d in 0 1 2 3 4; do
  printf "%s\t%s\t%s\n" $rep $(cut -d' ' -f1 /proc/loadavg) "$(nice $W/readsam/readsam $SAM $d)" | tee -a $W/readsam.tsv
done; done
echo READSAM DONE
