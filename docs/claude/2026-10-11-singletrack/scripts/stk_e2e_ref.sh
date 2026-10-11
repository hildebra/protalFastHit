#!/usr/bin/env bash
# stk_e2e_ref.sh: the reference tree's end-to-end suite with its own build on the same mini database as stk_e2e.sh, and
# the tests that fail in one and not the other.
set -uo pipefail
R=$HOME/stk/ref/build; C=$HOME/stk/cmp
export PROTAL_TEST_DB=$HOME/bpar/src/data/mini_db/protal_db PROTAL=$R/protal SIMULATE=$R/simulate_metagenomes
cd $HOME/stk/ref
taskset -c 0-3 nice -n 5 timeout 3000 python3 -m unittest tests/e2e/test_protal_e2e.py > $C/e2e.ref.log 2>&1; echo "ref suite exit $?"
tail -3 $C/e2e.ref.log
grep -E '^(ERROR|FAIL):' $C/e2e.ref.log | sort > $C/e2e.ref.fails
grep -E '^(ERROR|FAIL):' $C/e2e.log | grep -v WfaTraceback | sort > $C/e2e.work.fails
echo "failing in work only:"; comm -13 $C/e2e.ref.fails $C/e2e.work.fails
echo "failing in ref only:"; comm -23 $C/e2e.ref.fails $C/e2e.work.fails
echo "E2EREFDONE $(date +%T)"
