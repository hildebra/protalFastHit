#!/bin/bash
# Inflate speed of one-member gzip FASTQ: ISA-L on one thread vs rapidgzip on 2-4 (cores 0-3), alternated, 3 rounds.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/rgzbuild
B=$W/build
ZNG=$(dirname $(find $B -name zlib-ng.h -path "*zlib-ng*" | head -1))
echo "zlib-ng.h in $ZNG"
cd $W
g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -DLIBRAPIDARCHIVE_WITH_ISAL \
  -I src/src -I src/lib -I $B/isa-l-install/include -I $ZNG \
  -isystem src/lib/rapidgzip/src -isystem src/lib/rapidgzip-zlib -isystem $B/isa-l-install/include/isa-l \
  $SP/bench_gz.cpp src/src/IO/ParallelGzip.cpp \
  $B/isa-l-install/lib/libisal.a $(find $B -name "libz-ng.a" | head -1) -lzstd -o bench_gz 2> bench_compile.log
echo "compile exit $? ($(wc -l < bench_compile.log) lines of messages)"
F=~/gzpar/data/R1.member.gz
ls -la $F
cat $F > /dev/null  # page cache
for round in 1 2 3; do
  for t in 1 2 3 4; do
    taskset -c 0-3 ./bench_gz $F $t
  done
done
