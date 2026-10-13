#!/bin/bash
# Item 3: reading a gzip file from the disk rather than the page cache, before (0ab85f4) and with the read-ahead hints
# (the working tree): bench_pgz (../2026-10-13-gzip-members-regression/scripts/bench_pgz.cpp) on one member (pieces, a
# mapping; and on one thread, read()) and on 64 KiB members (by members, read()), each run after the file's pages were
# dropped from the page cache (dd iflag=nocache, posix_fadvise DONTNEED; fincore shows what stayed). WSL, cores 0-3, an
# ext4 file system on a virtual disk: a local disk, not a network file system. ROUNDS alternated rounds.
set -uo pipefail
R=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude
G=$R/2026-10-13-gzip-members-regression/scripts
B=~/pgzread
mkdir -p $B/old/src $B/data
tar -xf $B/old.tar -C $B/old/src   # git archive -o old.tar 0ab85f4 src lib/gzstream (Git Bash; not kept)
CXX="taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG"
$CXX -DPGZ_NEW -I $B/old/src/src -I $B/old/src/src/Utilities -I $B/old/src/lib $G/bench_pgz.cpp $B/old/src/src/IO/ParallelGzip.cpp \
  -lisal -lzstd -o $B/bench_old || exit 1
$CXX -DPGZ_NEW -I ~/pgzfix/src/src -I ~/pgzfix/src/src/Utilities -I ~/pgzfix/src/lib $G/bench_pgz.cpp \
  ~/pgzfix/src/src/IO/ParallelGzip.cpp -lisal -lzstd -o $B/bench_new || exit 1
cd $B/data
if [ ! -s one.gz ]; then
  taskset -c 0-3 nice -n 5 perl $G/make_fastq.pl 4000000 7 > R1.fq
  taskset -c 0-3 nice -n 5 gzip -6 -c R1.fq > one.gz
  taskset -c 0-3 nice -n 5 split -b 64K --filter='gzip -6' R1.fq > m64k.gz
  rm -f R1.fq
fi
uptime
for round in $(seq 1 ${ROUNDS:-3}); do
  for case in "one.gz 4" "one.gz 1" "m64k.gz 4"; do
    set -- $case
    for v in old new; do
      for cache in cold warm; do
        [ $cache = cold ] && dd if=$1 iflag=nocache count=0 status=none
        cached=$(fincore -n -b $1 2> /dev/null | awk '{ printf "%.0f", $1 / 1048576 }')
        printf 'r%s\t%s\t%s\t%s\tcached_mb %s\t' $round $1 $v $cache "${cached:-NA}"
        taskset -c 0-3 nice -n 5 $B/bench_$v $1 $2 | cut -f1-9
      done
    done
  done
done
uptime
