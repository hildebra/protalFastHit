#!/usr/bin/env bash
# asan.sh: the worktree's files (new.tar, the current test file) under AddressSanitizer and UBSan as CI builds them
# (Debug, -DPROTAL_NO_CLONES: WFA2's cloned functions compiled once too); every unit test, and the WFA2 replay of the
# pe and ont calls through the sanitized library. On cores 0-3.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfaspeed
WT=/mnt/c/Users/hildebra/Documents/locDev/protal-wfaspeed
W=$HOME/wfaspeed; A=$W/asan; C=$W/calls
pin() { taskset -c 0-3 nice -n 5 "$@"; }
rm -rf $A; mkdir -p $A; tar -xf $S/new.tar -C $A
for f in tests/test_WFA2Wrapper.cpp src/Alignment/WFA2Wrapper2.h; do cp $WT/$f $A/$f; done  # changed after new.tar
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
pin cmake -S $A -B $A/build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS="$FLAGS -DPROTAL_NO_CLONES" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS" > $W/asan.configure.log 2>&1 || { echo FAIL configure; exit 1; }
grep -q "PROTAL_NO_CLONES" $A/build/build.ninja && echo "PROTAL_NO_CLONES in build.ninja: $(grep -c 'DPROTAL_NO_CLONES' $A/build/build.ninja) lines"
pin cmake --build $A/build -j 4 --target protal_tests > $W/asan.build.log 2>&1 && echo "OK asan build" || { echo FAIL asan build; grep -E ' error' -A3 $W/asan.build.log | head; exit 1; }
# WFA2's objects compiled without the clones (no ifunc resolver symbols)
nm $(find $A/build -name libwfa_lib.a) 2>/dev/null | grep -c "\.resolver" | sed 's/^/resolver symbols in the sanitized libwfa_lib.a: /'
(cd $A && ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 pin ctest --test-dir build -j4 > $W/asan.ctest.log 2>&1); echo "asan ctest exit $?"
tail -4 $W/asan.ctest.log
H=$A; LIB=$(find $H/build -name libwfa_lib.a | head -1)
g++ -O1 -g $FLAGS -std=c++20 -I$H/src -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/SequenceUtils -I$H/src/Alignment \
  -I$H/src/Utilities -I$H/src/Taxonomy -I$H/lib -I$H/lib/robin -I$H/lib/tsl -I$H/lib/wfa2-lib -I$H/lib/gzstream -I$H/build/zlib-ng \
  -I$H/lib/zlib-ng -I$H/build/generated $S/wfareplay.cpp $LIB -o $W/replay/wfareplay_asan 2>&1 | grep -E ' error' | head
for s in pe pb ont; do
  ASAN_OPTIONS=detect_leaks=1 pin $W/replay/wfareplay_asan $C/$s.tsv 1 $W/replay/$s.asan classic singletrack > $W/replay/$s.asan.log 2>&1
  echo "asan replay $s: exit $?, $(grep -c ERROR $W/replay/$s.asan.log) sanitizer errors"
  for t in classic singletrack; do cmp -s $W/replay/$s.base.$t $W/replay/$s.asan.$t && echo "   $t identical to base" || echo "   $t DIFFERS from base"; done
done
echo "ASANDONE $(date +%T)"
