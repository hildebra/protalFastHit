#!/usr/bin/env bash
# queue1.sh: the fixed test file into ~/wfaspeed/new, rebuilt; the WFA tests; then replay.sh, lto.sh, identity.sh in turn.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfaspeed
WT=/mnt/c/Users/hildebra/Documents/locDev/protal-wfaspeed
W=$HOME/wfaspeed; N=$W/new
pin() { taskset -c 0-3 nice -n 5 "$@"; }
cp $WT/tests/test_WFA2Wrapper.cpp $N/tests/test_WFA2Wrapper.cpp
for f in src/Alignment/WFA2Wrapper2.h lib/wfa2-lib.cmake; do cmp -s $WT/$f $N/$f || { echo "WARNING: $f differs from the worktree"; }; done
pin cmake --build $N/build -j 4 --target protal_tests > $W/new.build2.log 2>&1 && echo "OK tests rebuilt" || { echo FAIL rebuild; grep -E ' error' -A3 $W/new.build2.log | head; exit 1; }
pin $N/build/tests/protal_tests --gtest_filter='WFA2Wrapper.*:WfaTraceback.*' > $W/unit.wfa2.log 2>&1; echo "wfa tests exit $?"
grep -E "^\[       OK|FAILED|^\[  PASSED|Failure" $W/unit.wfa2.log | head -30
echo "== replay ($(date +%T))"
bash $S/replay.sh 9 5
echo "== lto ($(date +%T))"
bash $S/lto.sh 9 5
echo "== identity ($(date +%T))"
bash $S/identity.sh
echo "QUEUE1DONE $(date +%T)"
