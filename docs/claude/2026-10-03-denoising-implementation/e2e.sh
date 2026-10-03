#!/usr/bin/env bash
# Mini database with the new binary, then the e2e tests and the mini-db unit tests (the GtdbBuildTest too).
set -uo pipefail
B=$HOME/protal-denoise
cd $B/src
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
PROTAL=$B/build/protal bash scripts/mini_db/build_mini_db.sh $B/mini > $B/mini.log 2>&1 || { echo "mini db failed"; tail -30 $B/mini.log; exit 1; }
echo "mini db done"
PROTAL_TEST_DB=$B/mini/protal_db PROTAL=$B/build/protal SIMULATE=$B/build/simulate_metagenomes \
  timeout 5400 python3 -m unittest tests/e2e/test_protal_e2e.py > $B/e2e.log 2>&1
echo "e2e exit $?"; tail -5 $B/e2e.log
PROTAL=$B/build/protal SIMULATE=$B/build/simulate_metagenomes PROTAL_TRAIN_PYTHON=$HOME/micromamba/envs/protal-db-build/bin/python \
  timeout 5400 python3 -m unittest scripts.mini_db.test_mini_db > $B/minitest.log 2>&1
echo "mini-db tests exit $?"; tail -5 $B/minitest.log
