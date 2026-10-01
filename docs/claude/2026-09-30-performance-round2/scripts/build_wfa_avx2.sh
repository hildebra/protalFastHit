#!/bin/bash
# Experiment build: WFA2-lib compiled for x86-64-v3 (lib/wfa2-lib.cmake builds one library with the
# baseline flags, for protal and protal_avx2 alike, so WFA2's own `#if __AVX2__` kernels are not
# in either binary). Applied to a copy of the tree; the repository is not changed.
set -e
source "$(dirname "$0")/env.sh"
rm -rf $PERF_DIR/src-wa; mkdir -p $PERF_DIR/src-wa
rsync -a --exclude '/data' --exclude '/build*/' --exclude '/.git' $PERF_DIR/src/ $PERF_DIR/src-wa/
cd $PERF_DIR/src-wa
printf '\ntarget_compile_options(wfa_lib PRIVATE -march=x86-64-v3)\n' >> lib/wfa2-lib.cmake
cd $PERF_DIR
cmake -S src-wa -B build-wa -G Ninja -DCMAKE_BUILD_TYPE=Release > build-wa.cmake.log 2>&1
cmake --build build-wa --target protal_avx2 -j "$(nproc)" > build-wa.log 2>&1
ls -la build-wa/protal_avx2
nm -C build-wa/protal_avx2 | grep -c "wavefront.*avx2" || true
