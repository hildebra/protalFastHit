#!/bin/bash
# Truncated one-member gzip files of several sizes and levels, read with 1 and 3 inflate threads (bench_gz):
# is the truncation always reported?
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
cd ~/rgzbuild
B=build
g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -DLIBRAPIDARCHIVE_WITH_ISAL \
  -I src/src -I src/lib -I $B/isa-l-install/include -I $B/zlib-ng \
  -isystem src/lib/rapidgzip/src -isystem src/lib/rapidgzip-zlib -isystem $B/isa-l-install/include/isa-l \
  $SP/bench_gz.cpp src/src/IO/ParallelGzip.cpp \
  $B/isa-l-install/lib/libisal.a $B/zlib-ng/libz-ng.a -lzstd -o bench_gz 2> bench_compile.log
echo "compile exit $?"
T=~/rgzbuild/trunc
rm -rf $T; mkdir -p $T
for mb in 0.1 0.5 1 2 5 20; do
  bytes=$(awk -v m=$mb 'BEGIN { printf "%d", m * 1048576 }')
  head -c $bytes ~/gzpar/data/R1.fq > $T/plain_$mb.txt
  for level in 1 6 9; do
    gzip -$level -c $T/plain_$mb.txt > $T/g_${mb}_$level.gz
    size=$(stat -c %s $T/g_${mb}_$level.gz)
    for frac in 0.25 0.5 0.9 0.999; do
      cut=$(awk -v s=$size -v f=$frac 'BEGIN { printf "%d", s * f }')
      head -c $cut $T/g_${mb}_$level.gz > $T/cut.gz
      one=$(PROTAL_INFLATE_THREADS=1 ./bench_gz $T/cut.gz 1 | awk -F'\t' '{print $9, $4}')
      three=$(PROTAL_INFLATE_THREADS=3 ./bench_gz $T/cut.gz 3 2>/dev/null | awk -F'\t' '{print $3, $9, $4}')
      echo "MB $mb level $level cut $frac: one [$one] three [$three]"
    done
  done
done
