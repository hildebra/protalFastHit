#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
R=/mnt/c/Users/hildebra/Documents/locDev/protal/lib
cd $D
taskset -c 0-3 nice -n 5 g++ -O2 -march=native -std=c++17 -pthread -DNDEBUG -DLIBRAPIDARCHIVE_WITH_ISAL -o bin/bfinder bfinder.cpp \
  -isystem $R/rapidgzip/src -isystem $R/rapidgzip-zlib -isystem ~/rgzbuild/build/isa-l-install/include/isa-l -isystem ~/rgzbuild/build/zlib-ng \
  ~/rgzbuild/build/isa-l-install/lib/libisal.a ~/rgzbuild/build/zlib-ng/libz-ng.a 2>&1 | grep -E "error" | head -20
for s in 1 4 16; do taskset -c 0-3 nice -n 5 bin/bfinder ~/gzpar/data/R1.member.gz logs/R1.bounds.tsv $s; done 2>&1 | tee logs/bfinder.log
