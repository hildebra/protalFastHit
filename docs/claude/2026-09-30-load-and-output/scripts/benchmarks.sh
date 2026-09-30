#!/bin/bash
# The stand-alone benchmarks: compression libraries on reads and a SAM (bench_compress.cpp), protal's
# per-gene start-up steps at db900 and GTDB r226 size (bench_genes.cpp, synthetic reference.map and
# unique_kmers.tsv, a sparse reference.fna), and sizing the gene sequences (bench_prepare.cpp).
#   benchmarks.sh READS_R1.fq.gz SAM
# Needs g++, zlib, libdeflate and zstd headers (libdeflate-dev, zlib1g-dev, libzstd-dev).
source "$(dirname "$0")/env.sh"
reads=$1; sam=$2
g++ -O2 -march=x86-64-v3 $here/bench_compress.cpp -o $OUT/bench_compress -lz -ldeflate -lzstd || exit 1
g++ -std=c++20 -O3 -march=x86-64-v3 $here/bench_prepare.cpp -o $OUT/bench_prepare -pthread || exit 1
inc="$(find $PROTAL_SRC/src -type d | sed 's/^/-I/' | tr '\n' ' ') -I$PROTAL_SRC/lib -I$PROTAL_SRC/lib/tsl -I$PROTAL_SRC/lib/robin -I$PROTAL_SRC/lib/gzstream"
gen_dir=$(ls -d $PERF_DIR/build*/generated 2>/dev/null | head -1)   # protal_version.h from a configured build
g++ -std=c++20 -O3 -march=x86-64-v3 -DNDEBUG $inc ${gen_dir:+-I$gen_dir} $here/bench_genes.cpp -o $OUT/bench_genes -lzstd -lz -fopenmp || exit 1

echo "== compression"
$OUT/bench_compress inflate $reads
$OUT/bench_compress deflate $sam

echo "== per-gene steps"
G=$OUT/genes; mkdir -p $G
# taxa and genes per taxon: 900 x 112 (db900's size); GTDB r226: 136,646 bacteria x 119, 6,968 archaea x 52
gen() {
  local name=$1; shift
  awk -v spec="$*" 'BEGIN { n = split(spec, a, " "); tax = 0; pos = 0
      for (i = 1; i <= n; i += 2) for (t = 0; t < a[i]; t++) { tax++
        for (g = 1; g <= a[i+1]; g++) { len = 900 + (tax * 7 + g * 13) % 300
          printf "%d\t%d\t%d\t%d\n", tax, g, pos, pos + len > "/dev/stderr"
          printf "%d\t%d\t%d\t%.6f\t%d\t%.6f\t%d\t%.6f\t%d\n", tax, g, 12, 0.0123, 3, 0.0031, 1, 0.0010, len - 30
          pos += len + 1 } }
      printf "%d\n", pos > "/dev/fd/3" }' 2> $G/$name.map > $G/$name.uk 3> $G/$name.size
  rm -f $G/$name.fna; truncate -s $(cat $G/$name.size) $G/$name.fna
  echo "$name: $(wc -l < $G/$name.map) genes, reference $(cat $G/$name.size) bytes"
}
gen db900like 900 112
gen gtdb 136646 119 6968 52
for name in db900like gtdb; do
  for rep in 1 2; do
    /usr/bin/time -f "%e s wall, %M kB max RSS" $OUT/bench_genes $G/$name.fna $G/$name.map $G/$name.uk $G/$name.header.sam 2> $G/time.txt
    echo "  $name rep $rep: $(tail -1 $G/time.txt); SAM header $(stat -c %s $G/$name.header.sam) bytes"
  done
done
rm -f $G/gtdb.header.sam $G/gtdb.fna

echo "== sizing gene sequences, 4M genes"
for rep in 1 2; do $OUT/bench_prepare 4000000; done
