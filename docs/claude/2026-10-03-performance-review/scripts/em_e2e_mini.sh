#!/usr/bin/env bash
# The end-to-end tests on the em build against a mini database (the one in ~/mt-work/head07c, built by 07c371b's tree).
set -uo pipefail
W=$HOME/mt-work/perf4; E=$W/em
cd $E
ls $HOME/mt-work/head07c/mini_db/protal_db | head -12
PROTAL_TEST_DB=$HOME/mt-work/head07c/mini_db/protal_db PROTAL=$E/build/protal SIMULATE=$E/build/simulate_metagenomes \
  python3 -m unittest -v tests/e2e/test_protal_e2e.py > $W/em.e2e_mini.log 2>&1; echo "e2e mini: $(grep -E '^Ran |^OK|^FAILED' $W/em.e2e_mini.log | tr '\n' ' ')"
grep -E '^(FAIL|ERROR):' $W/em.e2e_mini.log | head
