#!/bin/bash
# Experiment build: the same source as build-rel with Benchmark::Start/Stop empty (notimers.patch
# applied to a copy of the tree), to measure what the per-read stage timers cost.
set -e
source "$(dirname "$0")/env.sh"
rm -rf $PERF_DIR/src-nt; mkdir -p $PERF_DIR/src-nt
rsync -a --exclude '/data' --exclude '/build*/' --exclude '/.git' $PERF_DIR/src/ $PERF_DIR/src-nt/
cd $PERF_DIR/src-nt
perl -0pi -e 's/void Start\(bool new_sample=true\) \{\s*samplings \+= new_sample;\s*start_time = std::chrono::steady_clock::now\(\);\s*\}/void Start(bool new_sample=true) { samplings += new_sample; }/; s/void Stop\(\) \{\s*auto const stop_time = std::chrono::steady_clock::now\(\);\s*time_sum \+= duration_cast<std::chrono::nanoseconds>\(stop_time - start_time\)\.count\(\);\s*\}/void Stop() {}/' src/Utilities/Benchmark.h
grep -n "void Start\|void Stop" src/Utilities/Benchmark.h
cd $PERF_DIR
cmake -S src-nt -B build-nt -G Ninja -DCMAKE_BUILD_TYPE=Release > build-nt.cmake.log 2>&1
cmake --build build-nt --target protal_avx2 -j "$(nproc)" > build-nt.log 2>&1
ls -la build-nt/protal_avx2
