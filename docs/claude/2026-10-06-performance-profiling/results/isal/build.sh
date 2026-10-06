#!/bin/bash
# Build the benchmark and test programs against ~/isal/prefix (static ISA-L and zlib-ng).
set -e
cd ~/isal/bench
P=$HOME/isal/prefix
CFLAGS="-O2 -g -Wall -Wextra -Wno-unused-parameter -I$P/include"
for prog in bench_inflate bench_deflate test_semantics; do
  [ -f $prog.c ] || continue
  gcc $CFLAGS -o $prog $prog.c $P/lib/libisal.a $P/lib/libz-ng.a /usr/lib/x86_64-linux-gnu/libdeflate.a
  echo "built $prog"
done
