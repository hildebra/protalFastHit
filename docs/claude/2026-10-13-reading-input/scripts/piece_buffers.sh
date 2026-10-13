#!/bin/bash
# Item 4: ParallelGzip's output buffers against its speed. bench_pgz on one gzip member (pieces) at 2 and 4 threads,
# built with the working tree's ParallelGzip.cpp and with copies whose look-ahead is changed:
#   V0 as is: workers from gap = 2T + 4 pieces ahead of the reader's thread, up to ahead = gap + 2T
#   V1 ahead = gap + T;  V2 gap = T + 2, ahead = gap + 2T;  V3 gap = T + 2, ahead = gap + T
# Max RSS and wall time, warm page cache, cores 0-3, ROUNDS alternated rounds. Data: cold_reads.sh's one.gz.
set -uo pipefail
G=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
B=~/pgzread
S=~/pgzfix/src
mkdir -p $B/variants
GAP='if (impl.gap == 0) impl.gap = 2 \* size_t{impl.threads} + 4;'
AHEAD='impl.ahead = impl.gap + 2 \* size_t{impl.threads};'
grep -q "impl.gap = 2 \* size_t{impl.threads} + 4;" $S/src/IO/ParallelGzip.cpp || { echo "gap line not found"; exit 1; }
grep -q "impl.ahead = impl.gap + 2 \* size_t{impl.threads};" $S/src/IO/ParallelGzip.cpp || { echo "ahead line not found"; exit 1; }
cp $S/src/IO/ParallelGzip.cpp $B/variants/V0.cpp
sed "s/$AHEAD/impl.ahead = impl.gap + size_t{impl.threads};/" $S/src/IO/ParallelGzip.cpp > $B/variants/V1.cpp
sed "s/$GAP/if (impl.gap == 0) impl.gap = size_t{impl.threads} + 2;/" $S/src/IO/ParallelGzip.cpp > $B/variants/V2.cpp
sed -e "s/$GAP/if (impl.gap == 0) impl.gap = size_t{impl.threads} + 2;/" -e "s/$AHEAD/impl.ahead = impl.gap + size_t{impl.threads};/" \
  $S/src/IO/ParallelGzip.cpp > $B/variants/V3.cpp
for v in V0 V1 V2 V3; do
  cmp -s $B/variants/$v.cpp $S/src/IO/ParallelGzip.cpp && [ $v != V0 ] && { echo "$v unchanged"; exit 1; }
  # The copy compiles next to the sources it includes.
  cp $B/variants/$v.cpp $S/src/IO/ParallelGzip_$v.cpp
  taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -DPGZ_NEW -I $S/src -I $S/src/Utilities \
    -I $S/lib $G/bench_pgz.cpp $S/src/IO/ParallelGzip_$v.cpp -lisal -lzstd -o $B/bench_$v || exit 1
  rm -f $S/src/IO/ParallelGzip_$v.cpp
done
cd $B/data
cat one.gz > /dev/null
uptime
for round in $(seq 1 ${ROUNDS:-3}); do
  for t in 2 4; do
    for v in V0 V1 V2 V3; do
      printf 'r%s\t%s\t' $round $v
      taskset -c 0-3 nice -n 5 $B/bench_$v one.gz $t | cut -f1,6-9
    done
  done
done
uptime
