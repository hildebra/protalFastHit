#!/usr/bin/env bash
# The read-format benchmark of this report. Usage: bash run_bench.sh BUILD_DIR SOURCE_DIR READS.fq.gz WORK_DIR
# BUILD_DIR: a protal build (its vendored zlib-ng); READS: a BGZF FASTQ of simulate_metagenomes (here the first reads
# of 500,000 pairs, 150 bp HSXt), used three times over. Needs gzip and zstd.
set -euo pipefail
build=$(realpath "$1"); src=$(realpath "$2"); reads=$(realpath "$3"); work=$4
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$work" && cd "$work"
for copy in 1 2 3; do zcat "$reads"; done > reads.fq             # three copies: ~0.5 GB of FASTQ
for copy in 1 2 3; do cat "$reads"; done > reads.bgzf.fq.gz       # the simulator's BGZF (libdeflate level 6), thrice
gzip -6 -c reads.fq > reads.gzip6.fq.gz          # one gzip member, as sequencers' files are
for level in 1 3 9 19; do zstd -q -f -T0 -"$level" reads.fq -o reads.zstd"$level".fq.zst; done
zlib_ng=$(dirname "$(find "$build" -name 'libz-ng*.a' | head -1)")
g++ -O3 -DNDEBUG -std=c++20 -pthread -DZLIBNG_NATIVE_API -I"$src/src" -I"$src/src/IO" -I"$build/generated" \
    -I"$build/zlib-ng" -I"$src/lib/zlib-ng" "$here/bench_read_formats.cpp" -o bench_read_formats \
    -L"$zlib_ng" -lz-ng -ldeflate -lzstd
ls -l reads.*
./bench_read_formats reads.fq reads.bgzf.fq.gz reads.gzip6.fq.gz reads.zstd1.fq.zst reads.zstd3.fq.zst \
    reads.zstd9.fq.zst reads.zstd19.fq.zst --repeats 5 | tee bench.tsv
