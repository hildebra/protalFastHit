#!/bin/bash
# Follow-up #5: the reader alone (bench_reader.cpp) built against two protal checkouts, before and
# after taking FASTQ lines with memchr under the reader lock, on plain and gzipped read pairs, 1, 4
# and 8 threads, alternated. Also used for zlib against zlib-ng (the checkouts before and after it).
#   [ZNG_BUILD=<build dir>/zlib-ng] reader_lock.sh OLD_CHECKOUT NEW_CHECKOUT R1.fq R2.fq R1.fq.gz R2.fq.gz
# A checkout with lib/zlib-ng (vendored zlib-ng) needs ZNG_BUILD, the zlib-ng directory of a
# configured build of it (zlib-ng.h, libz-ng.a); older checkouts use the system zlib.
source "$(dirname "$0")/env.sh"
old=$1; new=$2; shift 2
for v in old new; do
  src=$old; [ $v = new ] && src=$new
  gz="-lz"
  if [ -d $src/lib/zlib-ng ]; then
    [ -f "$ZNG_BUILD/libz-ng.a" ] || { echo "$src vendors zlib-ng: set ZNG_BUILD to the zlib-ng directory of its build"; exit 2; }
    gz="-DWITH_GZFILEOP -DZLIBNG_NATIVE_API -I$ZNG_BUILD $ZNG_BUILD/libz-ng.a"
  fi
  g++ -std=c++20 -O3 -march=x86-64-v3 -DNDEBUG -fopenmp -I$src/src -I$src/src/IO -I$src/lib -I$src/lib/gzstream \
    $here/bench_reader.cpp $src/src/IO/FastxReader.cpp -o $OUT/bench_reader_$v $gz -ldeflate || exit 1
done
cat "$@" > /dev/null   # into the page cache
for rep in 1 2; do
  for kind in plain gz; do
    r1=$1; r2=$2; [ $kind = gz ] && { r1=$3; r2=$4; }
    for t in 1 4 8; do for v in old new; do
      echo "rep $rep $kind $v: $($OUT/bench_reader_$v $r1 $r2 $t)"
    done; done
  done
done
