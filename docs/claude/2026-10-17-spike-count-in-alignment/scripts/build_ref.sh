#!/usr/bin/env bash
# Build protal at HEAD (the tar git archive made, before the change) in ~/spikecount/ref, cores 0-3, niced, with ccache.
set -uo pipefail
TAR=$1
ROOT=$HOME/spikecount
rm -rf "$ROOT/ref"; mkdir -p "$ROOT/ref"
tar -xf "$TAR" -C "$ROOT/ref" || exit 1
cd "$ROOT/ref" || exit 1
start=$(date +%s)
taskset -c 0-3 nice -n 5 cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache > "$ROOT/ref_cmake.log" 2>&1 || { tail -30 "$ROOT/ref_cmake.log"; exit 1; }
taskset -c 0-3 nice -n 5 ninja -C build -j4 protal > "$ROOT/ref_ninja.log" 2>&1 || { tail -40 "$ROOT/ref_ninja.log"; exit 1; }
echo "ref build: $(( $(date +%s) - start )) s"; ./build/protal --version 2>&1 | head -2
