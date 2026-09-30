#!/bin/bash
# Follow-up: gzip inflate with zlib, zlib-ng and libdeflate (bench_inflate.cpp) on gzip files as
# sequencers write them (one member, not BGZF), alternated.
#   ZNG_BUILD=<build dir>/zlib-ng inflate_libs.sh FILE.gz [FILE.gz ...]
# ZNG_BUILD is the zlib-ng directory of a configured protal build (zlib-ng.h, libz-ng.a); zlib's
# headers (zlib1g-dev) are needed for the zlib variant only.
source "$(dirname "$0")/env.sh"
[ -f "$ZNG_BUILD/libz-ng.a" ] || { echo "set ZNG_BUILD to the zlib-ng directory of a protal build"; exit 2; }
g++ -O2 -DWITH_ZLIB $here/bench_inflate.cpp -o $OUT/inflate_zlib -lz || exit 1
g++ -O2 -DWITH_ZLIB_NG -DWITH_GZFILEOP -DZLIBNG_NATIVE_API -I$ZNG_BUILD $here/bench_inflate.cpp -o $OUT/inflate_zlib_ng $ZNG_BUILD/libz-ng.a || exit 1
g++ -O2 -DWITH_LIBDEFLATE $here/bench_inflate.cpp -o $OUT/inflate_libdeflate -ldeflate || exit 1
for f in "$@"; do
  cat $f > /dev/null   # into the page cache
  echo "== $(basename $f): $(stat -c %s $f) bytes"
  for rep in 1 2; do for v in zlib zlib_ng libdeflate; do $OUT/inflate_$v $f 2; done; done
done
