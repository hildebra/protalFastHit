#!/bin/bash
# Follow-up #4: GenomeLoader::LoadAllGenomes (bench_gene_arena.cpp) built against two protal
# checkouts, before and after the gene arena, on N synthetic genes of 900-1199 bytes (119 per taxon,
# a sparse reference.fna), 1 and 8 threads.
#   gene_arena.sh OLD_CHECKOUT NEW_CHECKOUT [N]
# A configured build of NEW_CHECKOUT under PERF_DIR provides protal_version.h.
source "$(dirname "$0")/env.sh"
old=$1; new=$2; n=${3:-2000000}
A=$OUT/arena; mkdir -p $A
if [ ! -s $A/genes$n.map ]; then
  awk -v n=$n 'BEGIN { pos = 0; for (i = 0; i < n; i++) { tax = int(i / 119) + 1; g = i % 119 + 1; len = 900 + (tax * 7 + g * 13) % 300
      printf "%d\t%d\t%d\t%d\n", tax, g, pos, pos + len; pos += len + 1 }
      printf "%d\n", pos > "/dev/stderr" }' > $A/genes$n.map 2> $A/genes$n.size
  rm -f $A/genes$n.fna; truncate -s $(cat $A/genes$n.size) $A/genes$n.fna
fi
echo "$(wc -l < $A/genes$n.map) genes, reference $(cat $A/genes$n.size) bytes"
gen_dir=$(ls -d $PERF_DIR/*/build*/generated $PERF_DIR/build*/generated 2>/dev/null | head -1)
for v in old new; do
  src=$old; [ $v = new ] && src=$new
  inc="$(find $src/src -type d | sed 's/^/-I/' | tr '\n' ' ') -I$src/lib -I$src/lib/tsl -I$src/lib/robin -I$src/lib/gzstream"
  g++ -std=c++20 -O3 -march=x86-64-v3 -DNDEBUG $inc ${gen_dir:+-I$gen_dir} $here/bench_gene_arena.cpp -o $A/bench_$v -lzstd -lz -ldeflate -fopenmp -pthread || exit 1
done
cat $A/genes$n.fna > /dev/null   # into the page cache
for rep in 1 2; do for t in 1 8; do for v in old new; do
  echo "  $v: $($A/bench_$v $A/genes$n.fna $A/genes$n.map $t)"
done; done; done
