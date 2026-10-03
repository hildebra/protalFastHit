#!/usr/bin/env bash
# After em_build.sh: all unit tests (ctest) and the end-to-end tests on the em build, against the V073 database.
set -uo pipefail
W=$HOME/mt-work/perf4; E=$W/em
while [ ! -f $W/c.pb90M_t6.diff ]; do sleep 10; done
cd $E
nice cmake --build build --target simulate_metagenomes -j 5 > $W/em.build2.log 2>&1 || { echo "FAIL simulate build"; exit 1; }
ctest --test-dir build --output-on-failure -j 4 > $W/em.ctest.log 2>&1; echo "ctest: $(grep -E 'tests passed|tests failed' $W/em.ctest.log | tail -2 | tr '\n' ' ')"
PROTAL_TEST_DB=$HOME/bench071/V073/protal_db PROTAL=$E/build/protal SIMULATE=$E/build/simulate_metagenomes \
  python3 -m unittest -v tests/e2e/test_protal_e2e.py > $W/em.e2e.log 2>&1; echo "e2e: $(grep -E '^Ran |^OK|^FAILED' $W/em.e2e.log | tr '\n' ' ')"
grep -E '^(FAIL|ERROR):' $W/em.e2e.log | head
