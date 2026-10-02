#!/usr/bin/env bash
# CI's sanitizer job on the one-binary tree: Debug, ASan + UBSan, protal_tests, ctest. Niced, -j 2.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/mt-work/oneb-asan; rm -rf $W; mkdir -p $W/src
git -C $REPO archive d11381f | tar -x -C $W/src
(cd $W/src && patch -p1 -s < ${PATCH:-$S/oneb_full.patch}) || { echo "FAIL patch"; exit 1; }
cd $W/src
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
nice cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS="$FLAGS" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS" > $W/configure.log 2>&1 || { echo "FAIL configure"; exit 1; }
nice cmake --build build-asan --target protal_tests -j 2 > $W/build.log 2>&1 || { echo "FAIL build"; grep -E 'error' -A3 $W/build.log | head -40; exit 1; }
echo "clones in the ASan tests binary: $(nm build-asan/tests/protal_tests | grep -c '\.arch_x86_64_v3$')"
(cd build-asan && ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 nice ctest --output-on-failure -j 2 > $W/ctest.log 2>&1); ct=$?
echo "ctest: $(grep -E 'tests passed|tests failed' $W/ctest.log)"
[ $ct -eq 0 ] && echo "OK asan" || { grep -E 'Failed|FAILED|ERROR: AddressSanitizer|runtime error' $W/ctest.log | head -20; echo "FAIL asan"; }
