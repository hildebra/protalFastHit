#!/usr/bin/env bash
# The corrected StopsWhenTheSharesAreStable test: rebuild the tests in the fin tree and run ctest.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/mt-work/perf4; E=$W/fin
cp $REPO/tests/test_SampleContext.cpp $E/tests/
cd $E && nice cmake --build build --target protal_tests -j 5 > $W/fin3.build.log 2>&1 && echo "OK fin3 build" || { echo "FAIL fin3 build"; grep -E 'error' -A3 $W/fin3.build.log | head -30; exit 1; }
$E/build/tests/protal_tests --gtest_filter='SampleContext.*' 2>&1 | tail -2
ctest --test-dir build -j 4 > $W/fin3.ctest.log 2>&1; echo "ctest: $(grep -E 'tests passed|tests failed' $W/fin3.ctest.log | tail -1)"
grep -E '\*\*\*Failed' $W/fin3.ctest.log | head -3
