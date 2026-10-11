#!/usr/bin/env bash
# queue4.sh: after queue3: the slab allocator inlined again (static FORCE_INLINE); the checks of queue3 again
# (unit tests, replay identity, timing and callgrind against base, new against new1, whole runs, e2e, ASan).
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfaspeed
WT=/mnt/c/Users/hildebra/Documents/locDev/protal-wfaspeed
W=$HOME/wfaspeed; N=$W/new; A=$W/asan; R=$W/replay; C=$W/calls
pin() { taskset -c 0-3 nice -n 5 "$@"; }
until grep -q QUEUE3DONE $S/queue3.log; do sleep 10; done
for s in pe pb ont; do for t in classic singletrack; do cp $R/cg.$s.$t.new.out $R/cg.$s.$t.new2.out; done; done
cp $R/cg.tsv $R/cg2.tsv
for f in lib/wfa2-lib/wavefront/wavefront_slab.c lib/wfa2-lib/wavefront/wavefront_slab.h lib/wfa2-lib.cmake src/Alignment/WFA2Wrapper2.h tests/test_WFA2Wrapper.cpp; do
  cp $WT/$f $N/$f; cp $WT/$f $A/$f
done
pin cmake --build $N/build -j 4 --target protal simulate_metagenomes protal_tests > $W/new.build4.log 2>&1 && echo "OK new rebuilt" || { echo FAIL rebuild; grep -E ' error' -A3 $W/new.build4.log | head; exit 1; }
pin $N/build/tests/protal_tests --gtest_filter='WFA2Wrapper.*:WfaTraceback.*' > $W/unit.wfa4.log 2>&1; echo "wfa tests exit $?"
grep -E "FAILED|^\[  PASSED|Failure" $W/unit.wfa4.log | head
(cd $N && pin ctest --test-dir build -j4 > $W/unit.all4.log 2>&1); echo "ctest exit $?"; grep -E "tests passed" $W/unit.all4.log
echo "== replay ($(date +%T))"
bash $S/replay.sh 9 5
echo "== new against new1 (without free_intact), callgrind ($(date +%T))"
tot() { callgrind_annotate --inclusive=yes --threshold=100 $1 2>/dev/null | grep -E ':wavefront_align ' | head -1 | sed -E 's/^ *([0-9,]+) .*/\1/' | tr -d ,; }
for s in pe pb ont; do for t in classic singletrack; do
  a=$(tot $R/cg.$s.$t.new1.out); b=$(tot $R/cg.$s.$t.new.out)
  awk -v s=$s -v t=$t -v a=$a -v b=$b 'BEGIN {printf "   %-4s %-12s new1 %12d, new %12d, %+.2f%%\n", s, t, a, b, 100*(b-a)/a}'
done; done
echo "== whole runs ($(date +%T))"
bash $S/identity.sh
echo "== asan again ($(date +%T))"
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
pin cmake --build $A/build -j 4 --target protal_tests > $W/asan.build4.log 2>&1 && echo "OK asan rebuilt" || { echo FAIL asan rebuild; exit 1; }
(cd $A && ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 pin ctest --test-dir build -j4 > $W/asan.ctest4.log 2>&1); echo "asan ctest exit $?"
grep -E "tests passed" $W/asan.ctest4.log
H=$A; LIB=$(find $H/build -name libwfa_lib.a | head -1)
g++ -O1 -g $FLAGS -std=c++20 -I$H/src -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/SequenceUtils -I$H/src/Alignment \
  -I$H/src/Utilities -I$H/src/Taxonomy -I$H/lib -I$H/lib/robin -I$H/lib/tsl -I$H/lib/wfa2-lib -I$H/lib/gzstream -I$H/build/zlib-ng \
  -I$H/lib/zlib-ng -I$H/build/generated $S/wfareplay.cpp $LIB -o $R/wfareplay_asan 2>&1 | grep -E ' error' | head
for s in pe pb ont; do
  ASAN_OPTIONS=detect_leaks=1 pin $R/wfareplay_asan $C/$s.tsv 1 $R/$s.asan classic singletrack > $R/$s.asan.log 2>&1
  echo "asan replay $s: exit $?, $(grep -c ERROR $R/$s.asan.log) sanitizer errors"
  for t in classic singletrack; do cmp -s $R/$s.base.$t $R/$s.asan.$t && echo "   $t identical to base" || echo "   $t DIFFERS from base"; done
done
echo "QUEUE4DONE $(date +%T)"
