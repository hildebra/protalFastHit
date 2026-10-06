#!/usr/bin/env bash
# cc8.sh <source.cpp> <output> [extra flags]: cc.sh against ~/perf8/work (or $H): compiles a bench with protal's flags and include
# directories; links what the header-only parts need (zstd from the system, zlib-ng from the build).
set -uo pipefail
H=${H:-$HOME/perf8/work}
src=$1; out=$2; shift 2
INC="-I$H/build/generated -I$H/src -I$H/lib -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/VarkitInterface -I$H/src/Alignment \
 -I$H/src/Utilities -I$H/src/Profiling -I$H/src/SequenceUtils -I$H/src/RandomForest -I$H/src/Taxonomy -I$H/src/SNPs -I$H/lib/robin \
 -I$H/lib/tsl -I$H/lib/cPMML/include -I$H/lib/wfa2-lib -I$H/lib/wfa2-lib/wavefront -I$H/lib/wfa2-lib/utils -I$H/lib/gzstream \
 -I$H/build/zlib-ng -I$H/lib/zlib-ng"
g++ -fopenmp -O3 -DNDEBUG -std=c++20 -pthread -march=x86-64 -mtune=generic -ffp-contract=off -DWITH_GZFILEOP -DZLIBNG_NATIVE_API \
  $INC "$@" "$src" -o "$out" -lzstd 2>&1 | grep -E "error" | head -20
ls -la "$out" 2>/dev/null || echo "no binary"
