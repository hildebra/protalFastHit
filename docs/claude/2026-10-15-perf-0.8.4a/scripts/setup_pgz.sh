#!/usr/bin/env bash
# Trial A setup: HEAD's source in ~/perf085/src, bench_pgz (ThreadedGzIstream over ParallelGzip), and long-read test
# files: the local ONT mock's FASTQ recompressed in one member by gzip -6, gzip -1 and pigz -6, plus ENA's own file.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/95e4f938-fd80-46a7-b960-30fecd6a2871/scratchpad
W=~/perf085
mkdir -p $W/src $W/data $W/bin
cd $W
if [ ! -f src/src/IO/ParallelGzip.cpp ]; then tar -xf $SP/src_head.tar -C $W/src; fi
cp $SP/bench_pgz.cpp $W/
taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -DPGZ_NEW -I src/src -I src/src/Utilities -I src/lib \
  bench_pgz.cpp src/src/IO/ParallelGzip.cpp -lisal -lzstd -o bin/bench_pgz 2> bench_compile.log
echo "compile exit $?"; head -20 bench_compile.log
cd $W/data
ONT=~/Mocks/reads/d6331_ont/d6331_ont_ERR15459768/ERR15459768.fq.gz
HIFI=~/Mocks/reads/d6331_hifi/d6331_hifi_SRR13128012/SRR13128012.fq.gz
[ -s ont.fq ] || taskset -c 0-3 nice -n 5 zcat $ONT > ont.fq
ls -la ont.fq
head -c 400 ont.fq; echo
awk 'NR%4==2{n++; b+=length($0)} END{printf "reads %d bases %d mean %.0f\n", n, b, b/n}' ont.fq
[ -s ont.g6.gz ] || taskset -c 0-3 nice -n 5 gzip -6 -c ont.fq > ont.g6.gz
[ -s ont.g1.gz ] || taskset -c 0-3 nice -n 5 gzip -1 -c ont.fq > ont.g1.gz
[ -s ont.p6.gz ] || taskset -c 0-3 nice -n 5 pigz -6 -p 4 -c ont.fq > ont.p6.gz
ln -sf $ONT ont.ena.gz
ln -sf $HIFI hifi.ena.gz
ls -laL
