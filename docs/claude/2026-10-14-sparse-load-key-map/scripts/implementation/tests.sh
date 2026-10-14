#!/bin/bash
# The unit tests (the touched suites first, timed, then all), the mini database and the end-to-end tests, on cores 4,5.
T=~/keysfix/new
cd $T
echo "== focus $(date +%T)"
start=$(date +%s.%N)
taskset -c 4,5 build/tests/protal_tests --gtest_filter='CompactKeys.*:IndexCodec.*:GridIndex.*:GridLevels.*:PackedIndex.*:Index.*:Database.*' > unit_focus.log 2>&1
echo "focus exit=$? $(echo "$(date +%s.%N) - $start" | bc) s"
grep -E "^\[  (PASSED|FAILED) |^\[       OK \] (CompactKeys|IndexCodec\.KeysAre)|FAILED" unit_focus.log | head -40
echo "== all unit tests $(date +%T)"
taskset -c 4,5 ctest --test-dir build --output-on-failure -j2 > unit_all.log 2>&1
echo "ctest exit=$?"
grep -E "tests passed|tests failed|Total Test time|\*\*\*Failed|Failed" unit_all.log | head -20
echo "== mini database $(date +%T)"
PROTAL=$T/build/protal taskset -c 4,5 bash src/scripts/mini_db/build_mini_db.sh $T/mini_db > mini_db.log 2>&1
echo "mini db exit=$?"
echo "== e2e $(date +%T)"
cd src
start=$(date +%s)
PROTAL_TEST_DB=$T/mini_db/protal_db PROTAL=$T/build/protal SIMULATE=$T/build/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
  taskset -c 4,5 python3 -m unittest tests/e2e/test_protal_e2e.py > ../e2e.log 2>&1
echo "e2e exit=$? $(( $(date +%s) - start )) s"
tail -4 ../e2e.log
grep -E "^(FAIL|ERROR):" ../e2e.log | head -20
echo "== done $(date +%T)"
