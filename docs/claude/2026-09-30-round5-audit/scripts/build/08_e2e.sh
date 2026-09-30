#!/bin/bash
# docs/testing.md's commands, literally where possible (a fake nproc caps the justfile's -j$(nproc)
# at 2; everything pinned to 2 cores): CI's Python step, just mini-db-test, just model-test,
# just e2e (baseline, fresh mini DB), and the e2e suite against protal_avx2 (not in CI or just e2e).
set -u
A=~/audit6/build
S=$A/src
export PATH=$A/fakebin:$PATH
cd $S
echo "== CI step: python3 -m unittest -v scripts/mini_db/test_mini_db.py"
{ time taskset -c 0,1 python3 -m unittest -v scripts/mini_db/test_mini_db.py ; } > $A/ci_minidb.log 2>&1; echo "rc=$?"; grep -E '^(Ran|OK|FAILED)' $A/ci_minidb.log
echo "== just mini-db-test"; taskset -c 0,1 just mini-db-test > $A/just_minidbtest.log 2>&1; echo "rc=$?"; grep -E '^(Ran|OK|FAILED)' $A/just_minidbtest.log
echo "== just model-test (scikit-learn installed? $(python3 -c 'import sklearn; print(sklearn.__version__)' 2>&1 | tail -1))"
taskset -c 0,1 just model-test > $A/just_modeltest.log 2>&1; echo "rc=$?"; grep -E '^(Ran|OK|FAILED)|Error|skipped' $A/just_modeltest.log | head -5
echo "== just e2e (builds data/mini_db with build/protal, then the suite)"
{ time taskset -c 0,1 just e2e ; } > $A/just_e2e.log 2>&1; echo "rc=$?"
grep -E '^(Ran|OK|FAILED)' $A/just_e2e.log | tail -3; grep -E '^(FAIL|ERROR):' $A/just_e2e.log | head; grep -E '^real' $A/just_e2e.log
echo "== e2e against protal_avx2 (fresh data/mini_db)"
{ time PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal_avx2 SIMULATE=build/simulate_metagenomes \
    taskset -c 0,1 python3 -m unittest tests/e2e/test_protal_e2e.py ; } > $A/e2e_avx2.log 2>&1; echo "rc=$?"
grep -E '^(Ran|OK|FAILED)' $A/e2e_avx2.log | tail -3; grep -E '^(FAIL|ERROR):' $A/e2e_avx2.log | head; grep -E '^real' $A/e2e_avx2.log
echo "== protal --help / --full_help / --map_help exit codes (docs/installation.md 'Checking the installation')"
for o in --help --full_help --map_help; do $A/prefix/bin/protal $o > /dev/null 2>&1; echo "$o rc=$?"; done
echo DONE
