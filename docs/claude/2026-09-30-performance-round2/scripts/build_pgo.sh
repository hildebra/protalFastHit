#!/bin/bash
# Experiment: profile-guided optimisation of the current source (GCC -fprofile-generate / -fprofile-use), trained
# on 300k pairs of mix and 150k of w900, one thread. C (WFA2, zlib-ng) and C++ are both profiled.
set -e
source "$(dirname "$0")/env.sh"
cd $PERF_DIR
rm -rf src-pgo pgo-data build-pgo; mkdir -p pgo-data
rsync -a --exclude '/data' --exclude '/build*/' --exclude '/.git' src/ src-pgo/
GEN="-fprofile-generate=$PERF_DIR/pgo-data -fprofile-update=single"
cmake -S src-pgo -B build-pgo -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="$GEN" -DCMAKE_CXX_FLAGS="$GEN" > build-pgo.cmake.log 2>&1
cmake --build build-pgo --target protal_avx2 -j "$(nproc)" > build-pgo.log 2>&1
cp build-pgo/protal_avx2 pgo-gen-protal_avx2; echo "instrumented build done"
mkdir -p pgo-train; rm -rf pgo-train/*
head -n $((300000*4)) reads/mix_plain/mix_R1.fq > pgo-train/mix_R1.fq; head -n $((300000*4)) reads/mix_plain/mix_R2.fq > pgo-train/mix_R2.fq
head -n $((150000*4)) reads/w900_plain/w900_R1.fq > pgo-train/w900_R1.fq; head -n $((150000*4)) reads/w900_plain/w900_R2.fq > pgo-train/w900_R2.fq
for s in mix w900; do
  $PERF_DIR/pgo-gen-protal_avx2 --db db900n -1 pgo-train/${s}_R1.fq -2 pgo-train/${s}_R2.fq -o pgo-train/out_$s -t 1 --no_qcmsa > pgo-train/$s.log 2>&1
done
echo "training done: $(ls pgo-data | wc -l) profile files"
USE="-fprofile-use=$PERF_DIR/pgo-data -fprofile-correction -Wno-missing-profile -Wno-coverage-mismatch"
cmake -S src-pgo -B build-pgo -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="$USE" -DCMAKE_CXX_FLAGS="$USE" > build-pgo.cmake.log 2>&1
cmake --build build-pgo --target protal_avx2 -j "$(nproc)" > build-pgo.log 2>&1
ls -la build-pgo/protal_avx2
