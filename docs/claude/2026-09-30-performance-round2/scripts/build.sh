#!/bin/bash
# Release binaries as shipped (build-rel) and the same flags plus debug info and frame pointers
# for callgrind (build-prof). The source is copied to $PERF_DIR/src first: on a case-insensitive
# file system (NTFS under WSL) cPMML's "options.h" would pick up protal's Options.h.
set -e
source "$(dirname "$0")/env.sh"
mkdir -p $PERF_DIR
rsync -a --delete --exclude '/data' --exclude '/build*/' --exclude '/.git' $PROTAL_SRC/ $PERF_DIR/src/
cd $PERF_DIR
cmake -S src -B build-rel -G Ninja -DCMAKE_BUILD_TYPE=Release > build-rel.cmake.log 2>&1
cmake --build build-rel --target protal protal_avx2 simulate_metagenomes -j "$(nproc)" > build-rel.log 2>&1
cmake -S src -B build-prof -G Ninja -DCMAKE_BUILD_TYPE=Release -DMY_FLAGS="-g -fno-omit-frame-pointer" > build-prof.cmake.log 2>&1
cmake --build build-prof --target protal_avx2 -j "$(nproc)" > build-prof.log 2>&1
ls -la build-rel/protal_avx2 build-rel/simulate_metagenomes build-prof/protal_avx2
