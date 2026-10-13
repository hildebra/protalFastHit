#!/bin/bash
# Item 4: how many output buffers ParallelGzip makes (TakeBuffer finding the pool empty), their bytes, and how many
# live at most (taken and not given back), on one gzip member at 2 and 4 threads: a copy of the working tree's
# ParallelGzip.cpp with counters printed when the reader is destroyed. Cores 0-3, warm page cache.
set -uo pipefail
G=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
B=~/pgzread
S=~/pgzfix/src
C=$S/src/IO/ParallelGzip_count.cpp
perl -0 -pe '
  s/(std::unique_ptr<Bytes> TakeBuffer\(\) \{\n)/$1            { std::lock_guard<std::mutex> l(pool_mutex); if (++live > max_live) max_live = live; }\n/ or die "take";
  s/(return std::make_unique<Bytes>\(\);)/{ std::lock_guard<std::mutex> l(pool_mutex); made++; } $1/ or die "make";
  s/(void GiveBack\(std::unique_ptr<Bytes>& bytes\) \{\n            if \(!bytes\) return;\n)/$1            { std::lock_guard<std::mutex> l(pool_mutex); live--; }\n/ or die "give";
  s/(~Impl\(\) \{\n)/$1            size_t pooled = 0; for (auto const& b : pool) pooled += b->capacity;\n            std::fprintf(stderr, "buffers made %zu, live at most %zu, pooled at the end %zu (%.0f MB of capacity), ahead %zu slots, gap %zu\\n", made, max_live, pool.size(), pooled \/ 1e6, ahead, gap);\n/ or die "dtor";
  s/(std::vector<std::unique_ptr<Bytes>> pool;)/$1\n        size_t made = 0, live = 0, max_live = 0;/ or die "fields";
' $S/src/IO/ParallelGzip.cpp > $C || exit 1
sed -i 's/#include <cstring>/#include <cstring>\n#include <cstdio>/' $C
taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -DPGZ_NEW -I $S/src -I $S/src/Utilities -I $S/lib \
  $G/bench_pgz.cpp $C -lisal -lzstd -o $B/bench_count || exit 1
rm -f $C
cd $B/data
for t in 2 4; do for round in 1 2; do
  printf 'r%s\t' $round; taskset -c 0-3 nice -n 5 $B/bench_count one.gz $t 2>&1 | tr '\n' '\t' | cut -c1-400; echo
done; done
