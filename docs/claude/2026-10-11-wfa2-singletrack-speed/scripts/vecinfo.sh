#!/usr/bin/env bash
# vecinfo.sh: which loops of WFA2's hot files GCC vectorizes, with protal's flags at x86-64 (baseline), x86-64-v2 and
# x86-64-v3 (-fopt-info-vec-optimized / -missed for the loops in the hot functions).
set -uo pipefail
L=$HOME/wfasimd/src/lib/wfa2-lib
cd $L
for arch in x86-64 x86-64-v2 x86-64-v3; do
  echo "== $arch"
  for f in wavefront/wavefront_heuristic.c wavefront/wavefront_compute_affine.c wavefront/wavefront_compute.c wavefront/wavefront_extend_kernels.c wavefront/wavefront_termination.c; do
    gcc -D_FILE_OFFSET_BITS=64 -I. -Iwavefront -Iutils -fopenmp -O3 -DNDEBUG -march=$arch -fopt-info-vec-optimized -c $f -o /tmp/vec.o 2>&1 \
      | grep -E 'optimized: loop vectorized' | sed -E "s#^.*/##; s/: optimized: loop vectorized using / /; s/ vectors/v/" | sort -u | tr '\n' ' '
    echo "   <- $f"
  done
done
# Line numbers of the hot loops, to read the list above
grep -n "PRAGMA_LOOP_VECTORIZE" wavefront/wavefront_heuristic.c wavefront/wavefront_compute_affine.c | head -20
