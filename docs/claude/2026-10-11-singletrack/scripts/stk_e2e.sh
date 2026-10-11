#!/usr/bin/env bash
# stk_e2e.sh: the end-to-end tests (tests/e2e/test_protal_e2e.py) of the work tree on the mini database, on cores 0-3;
# first the new WfaTracebackTest alone, then the whole suite.
set -uo pipefail
K=$HOME/stk/work/build; C=$HOME/stk/cmp
export PROTAL_TEST_DB=$HOME/bpar/src/data/mini_db/protal_db PROTAL=$K/protal SIMULATE=$K/simulate_metagenomes
cd $HOME/stk/work
taskset -c 0-3 nice -n 5 timeout 1200 python3 -m unittest tests.e2e.test_protal_e2e.WfaTracebackTest -v > $C/e2e.wfa.log 2>&1; echo "WfaTracebackTest exit $?"
tail -8 $C/e2e.wfa.log
taskset -c 0-3 nice -n 5 timeout 3000 python3 -m unittest tests/e2e/test_protal_e2e.py > $C/e2e.log 2>&1; echo "suite exit $?"
tail -6 $C/e2e.log
echo "E2EDONE $(date +%T)"
