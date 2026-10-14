#!/bin/bash
# protal at HEAD without this session's changes, in ~/keysfix/base, on cores 4,5 (the baseline of the timings).
set -e
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cf462e00-4c32-4801-b980-f96668bbdc3a/scratchpad/fix
T=~/keysfix/base
mkdir -p $T/src
if [ ! -f $T/src/CMakeLists.txt ]; then tar -xf $SP/head.tar -C $T/src; fi
cd $T
if [ ! -f build/build.ninja ]; then
  cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache > cmake.log 2>&1
fi
taskset -c 4,5 nice -n 5 ninja -C build -j2 protal > ninja.log 2>&1 || { tail -40 ninja.log; exit 1; }
ls -la build/protal
