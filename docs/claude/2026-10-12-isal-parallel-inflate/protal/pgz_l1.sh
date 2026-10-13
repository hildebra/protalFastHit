#!/bin/bash
# Level-1 random FASTQ (as the unit test's RandomFastq) read in 16 KiB pieces at gap 1: how pieces were decoded.
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/pgzdbg
P=/mnt/c/Users/hildebra/Documents/locDev/protal
mkdir -p $W/src/src && rsync -a --checksum --no-times $P/src/ $W/src/src/
cd $W
g++ -O2 -g -std=c++20 -march=x86-64 -pthread -I src/src -I src/src/Utilities $SP/pgz_counts.cpp src/src/IO/ParallelGzip.cpp -lisal -lzstd -o pgz_counts 2> counts_compile.log || { head counts_compile.log; exit 1; }
mkdir -p $W/cnt && cd $W/cnt
perl -e 'srand(21); for $i (0..19999) { $l = 100 + int(rand(101)); print "\@q$i/1\n"; print join("", map { (qw(A C G T))[int(rand(4))] } 1..$l), "\n+\n"; print join("", map { chr(35 + int(rand(40))) } 1..$l), "\n"; }' > rand.fq
for lv in 1 6; do gzip -$lv -c rand.fq > rand$lv.gz; done
ls -la rand1.gz rand6.gz
for f in rand1.gz rand6.gz; do for t in 2 3 5; do taskset -c 0-3 $W/pgz_counts $f $t 16384 1; done; done
