#!/usr/bin/env bash
# build.sh - copies the Windows checkout into ~/samcombine/src and builds protal and the unit tests there
# (docs/claude/2026-10-07-sam-combine). 4 cores, ccache. SRC_WIN: the checkout as WSL sees it.
set -eu
SRC_WIN=${SRC_WIN:-/mnt/c/Users/hildebra/Documents/locDev/protal}
DST=${DST:-$HOME/samcombine/src}
mkdir -p "$DST"
# --checksum --no-times: a changed file gets the copy's time, so ninja rebuilds what depends on it.
rsync -a --checksum --no-times --delete --exclude '/data' --exclude '/build*/' --exclude '/local' --exclude '.git' \
      --exclude '__pycache__' "$SRC_WIN/" "$DST/"
cd "$DST"
if ! grep -q "PROTAL_BUILD_TESTS:BOOL=ON" build/CMakeCache.txt 2>/dev/null; then
  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON \
        -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache > build_configure.log 2>&1
fi
taskset -c 0-3 nice -n 5 ninja -C build -j4 2>&1 | tail -25
ls -la build/protal
