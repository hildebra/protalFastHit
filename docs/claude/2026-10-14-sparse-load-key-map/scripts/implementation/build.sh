#!/bin/bash
# Build HEAD + this session's changes (overlay) with the unit tests in ~/keysfix/new, on cores 4,5.
set -e
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cf462e00-4c32-4801-b980-f96668bbdc3a/scratchpad/fix
T=~/keysfix/new
mkdir -p $T/src
if [ ! -f $T/src/CMakeLists.txt ]; then tar -xf $SP/head.tar -C $T/src; fi
rsync -a --checksum --no-times $SP/overlay/ $T/src/
cd $T
if [ ! -f build/build.ninja ]; then
  cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON \
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache > cmake.log 2>&1
fi
taskset -c 4,5 nice -n 5 ninja -C build -j2 protal simulate_metagenomes protal_tests > ninja.log 2>&1 || { grep -B3 -A15 " error\|error:" ninja.log | head -120; exit 1; }
tail -2 ninja.log
ls -la build/protal build/tests/protal_tests 2>/dev/null || find build -name protal_tests -type f
