#!/usr/bin/env bash
# stk_build.sh: ~/stk/ref = git archive of 7d5df4b (WFA2-lib v2.3.6); ~/stk/work = the same with the Singletrack patch
# (st.patch: the fork's WFA2-lib, --wfa_traceback, build_settings.tsv). Release builds with the unit tests, on cores 0-3.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
PATCH=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/st.patch
W=$HOME/stk
mkdir -p $W
for N in ref work; do
  E=$W/$N
  if [ ! -d $E/build ]; then
    rm -rf $E; mkdir -p $E
    git -C $REPO archive 7d5df4b | tar -x -C $E
    if [ $N = work ]; then
      (cd $E && git apply --whitespace=nowarn $PATCH) || { echo "FAIL patch"; exit 1; }
    fi
    (cd $E && taskset -c 0-3 nice -n 5 cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/$N.configure.log 2>&1) \
      || { echo "FAIL configure $N"; tail -20 $W/$N.configure.log; exit 1; }
  fi
  (cd $E && taskset -c 0-3 nice -n 5 cmake --build build -j 4 > $W/$N.build.log 2>&1) && echo "OK build $N" \
    || { echo "FAIL build $N"; grep -E 'error|Error' -A4 $W/$N.build.log | head -60; }
done
echo "BUILDDONE $(date +%T)"
