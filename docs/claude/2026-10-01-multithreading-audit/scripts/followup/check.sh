#!/usr/bin/env bash
# check.sh NAME [PATCH] [BASE]: a clean tree (git archive BASE, default HEAD, plus PATCH) in ~/mt-work/NAME,
# Release build with tests, ctest, the mini database built by that binary, the e2e tests. Niced.
# Steps can be skipped: SKIP_E2E=1. Logs in ~/mt-work/NAME/*.log; the last line says what passed.
set -uo pipefail
NAME=$1; PATCH=${2:-}; BASE=${3:-HEAD}
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/mt-work/$NAME
rm -rf $W; mkdir -p $W/src
git -C $REPO archive $BASE | tar -x -C $W/src
if [ -n "$PATCH" ]; then (cd $W/src && patch -p1 -s < $PATCH) || { echo "FAIL patch"; exit 1; }; fi
cd $W/src
nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/configure.log 2>&1 || { echo "FAIL configure"; tail -20 $W/configure.log; exit 1; }
nice cmake --build build --target protal simulate_metagenomes protal_tests -j 4 > $W/build.log 2>&1 || { echo "FAIL build"; grep -E 'error|Error' $W/build.log | head -30; exit 1; }
(cd build && nice ctest --output-on-failure -j 4 > $W/ctest.log 2>&1); ct=$?
echo "ctest: $(grep -E 'tests passed|tests failed' $W/ctest.log)"
[ $ct -eq 0 ] || { grep -E 'Failed|FAILED' $W/ctest.log | head -20; echo "FAIL ctest"; exit 1; }
[ "${SKIP_E2E:-0}" = 1 ] && { echo "OK (unit only)"; exit 0; }
PROTAL=$W/src/build/protal nice bash scripts/mini_db/build_mini_db.sh $W/mini_db > $W/mini_db.log 2>&1 || { echo "FAIL mini_db"; tail -20 $W/mini_db.log; exit 1; }
PROTAL_TEST_DB=$W/mini_db/protal_db PROTAL=$W/src/build/protal SIMULATE=$W/src/build/simulate_metagenomes \
  nice python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e.log 2>&1; e2=$?
echo "e2e: $(tail -3 $W/e2e.log | tr '\n' ' ')"
[ $e2 -eq 0 ] || { grep -E '^(FAIL|ERROR):' $W/e2e.log | head -20; echo "FAIL e2e"; exit 1; }
echo "OK $NAME: unit and e2e passed"
