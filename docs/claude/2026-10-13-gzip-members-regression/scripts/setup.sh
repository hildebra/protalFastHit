#!/bin/bash
# Build bench_gz (ThreadedGzIstream over ParallelGzip, commit 9c27756) and make test gzip files of one Illumina-like
# FASTQ: one member (gzip -6), and members of 64 KiB, 1 MiB and 8 MiB of FASTQ each, concatenated.
set -uo pipefail
SP=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
W=~/pgzdiag
mkdir -p $W/src $W/data
cd $W
# src.tar made beforehand from Git Bash: git archive -o scripts/src.tar 9c27756 src lib/gzstream (not kept)
tar -xf $SP/src.tar -C $W/src
cp /mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-12-isal-parallel-inflate/protal/bench_gz2.cpp $W/
taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -I src/src -I src/src/Utilities -I src/lib \
  bench_gz2.cpp src/src/IO/ParallelGzip.cpp -lisal -lzstd -o bench_gz 2> bench_compile.log
echo "compile exit $?"; head -5 bench_compile.log
cd $W/data
if [ ! -s R1.fq ]; then
  taskset -c 0-3 nice -n 5 perl $SP/make_fastq.pl 4000000 7 > R1.fq
fi
ls -la R1.fq
[ -s one.gz ] || taskset -c 0-3 nice -n 5 gzip -6 -c R1.fq > one.gz
[ -s m64k.gz ] || taskset -c 0-3 nice -n 5 split -b 64K --filter='gzip -6' R1.fq > m64k.gz
[ -s m1m.gz ] || taskset -c 0-3 nice -n 5 split -b 1M --filter='gzip -6' R1.fq > m1m.gz
[ -s m8m.gz ] || taskset -c 0-3 nice -n 5 split -b 8M --filter='gzip -6' R1.fq > m8m.gz
ls -la
for f in one m64k m1m m8m; do printf '%s ' $f; head -c 18 $f.gz | od -A n -t x1; done
