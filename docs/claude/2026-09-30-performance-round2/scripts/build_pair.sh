#!/bin/bash
# Two builds of the same tree for byte-identical comparisons: BASE, the committed HEAD (git archive), and NEW, the
# working tree, each in $PERF_DIR/pair/{base,new} with release binaries (protal, protal_avx2, simulate_metagenomes);
# NEW also with the unit tests. build_pair.sh [base|new|both]
set -e
source "$(dirname "$0")/env.sh"
what=${1:-both}
repo=$PROTAL_SRC
mkdir -p $PERF_DIR/pair
if [ $what != new ]; then
  rm -rf $PERF_DIR/pair/base; mkdir -p $PERF_DIR/pair/base/src
  (cd $repo && git archive HEAD) | tar -x -C $PERF_DIR/pair/base/src
  (cd $PERF_DIR/pair/base && cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 &&
   cmake --build build --target protal protal_avx2 simulate_metagenomes -j "$(nproc)" > build.log 2>&1; echo "base rc=$?") &
fi
if [ $what != base ]; then
  mkdir -p $PERF_DIR/pair/new/src
  rsync -a --delete --exclude '/data' --exclude '/build*/' --exclude '/.git' $repo/ $PERF_DIR/pair/new/src/
  (cd $PERF_DIR/pair/new && cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > cmake.log 2>&1 &&
   cmake --build build --target protal protal_avx2 simulate_metagenomes protal_tests -j "$(nproc)" > build.log 2>&1; echo "new rc=$?") &
fi
wait
