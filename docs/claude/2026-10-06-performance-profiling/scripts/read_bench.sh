#!/usr/bin/env bash
# read_bench.sh: protal's gzip input and BGZF output before (HEAD: zlib-ng streams, libdeflate BGZF at level 6) and
# after (ISA-L) the change, on the 477 MB FASTQ of the zstd-reads report (~/zstdbench: reads.fq, reads.gzip6.fq.gz,
# a single gzip -6 member; reads.bgzf.fq.gz, libdeflate's BGZF; reads.zstd3.fq.zst). The read benchmark is that
# report's bench_read_formats.cpp (ThreadedGzStream as protal reads its input); bgzf_write.cpp times bgzf::CompressFile.
# The two builds' runs alternate. Sources: ~/perf8/ref (HEAD) and ~/perf8/isal; ISA-L in ~/isal/prefix.
set -euo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
HERE=$REPO/docs/claude/2026-10-06-performance-profiling/scripts
BENCH=$REPO/docs/claude/2026-10-05-zstd-reads/bench_read_formats.cpp
D=$HOME/zstdbench; W=$HOME/perf8/readbench; I=$HOME/isal/prefix
mkdir -p $W && cd $W
old=$HOME/perf8/ref; new=$HOME/perf8/isal
zng=$(dirname "$(find $old/build -name 'libz-ng*.a' | head -1)")
for prog in bench_read_formats:$BENCH bgzf_write:$HERE/bgzf_write.cpp; do
  name=${prog%%:*}; src=${prog#*:}
  g++ -O3 -DNDEBUG -std=c++20 -pthread -DZLIBNG_NATIVE_API -I$old/src -I$old/src/IO -I$old/build/generated -I$old/build/zlib-ng \
      -I$old/lib/zlib-ng $src -o $name.old -L$zng -lz-ng -ldeflate -lzstd
  g++ -O3 -DNDEBUG -std=c++20 -pthread -I$new/src -I$new/src/IO -I$I/include $src -o $name.new $I/lib/libisal.a -lzstd
done
echo "== BGZF written by bgzf::CompressFile ($(date +%T), load $(cut -d' ' -f1-3 /proc/loadavg))"
for t in 1 6; do
  for v in old new; do echo "-- $v, $t threads"; ./bgzf_write.$v $D/reads.fq $W/reads.$v.bgzf.fq.gz $t --repeats 3; done
done
echo "-- zstd -3, one thread (as the simulator's zstd reads), 3 runs: wall s"
TIMEFORMAT="    %R s wall, %U s user"
for r in 1 2 3; do { time zstd -3 -T1 -q -f -c $D/reads.fq > $W/reads.zstd3.fq.zst; } 2>&1; done
ls -l $W/reads.zstd3.fq.zst | awk '{print "    " $5 " bytes"}'
echo "== read by ThreadedGzStream ($(date +%T), load $(cut -d' ' -f1-3 /proc/loadavg))"
for round in 1 2; do
  for v in old new; do
    echo "-- $v, round $round"
    ./bench_read_formats.$v $D/reads.gzip6.fq.gz $D/reads.bgzf.fq.gz $W/reads.new.bgzf.fq.gz $D/reads.zstd3.fq.zst --repeats 5
  done
done
echo "== the content: $(zcat $W/reads.new.bgzf.fq.gz | md5sum | cut -c1-12) (ISA-L BGZF), $(md5sum < $D/reads.fq | cut -c1-12) (reads.fq)"
echo "READBENCHDONE $(date +%T)"
