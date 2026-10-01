#!/bin/bash
# The exact content of a commit, to build and test it alone: HEAD (git archive) with the patches applied, in
# $PERF_DIR/pair/mine, with the tests. build_clean.sh FULL.patch [ZERO_CONTEXT.patch]
# (when other work is uncommitted in the same working tree, as with two sessions, a working-tree copy
# would test both).
set -e
source "$(dirname "$0")/env.sh"
repo=$PROTAL_SRC
rm -rf $PERF_DIR/pair/mine; mkdir -p $PERF_DIR/pair/mine/src
(cd $repo && git archive HEAD) | tar -x -C $PERF_DIR/pair/mine/src
cd $PERF_DIR/pair/mine/src
git apply "$1"
[ -n "$2" ] && git apply --unidiff-zero "$2"
cd $PERF_DIR/pair/mine
cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > cmake.log 2>&1
cmake --build build --target protal protal_avx2 simulate_metagenomes protal_tests -j "$(nproc)" > build.log 2>&1
echo "mine rc=$?"
