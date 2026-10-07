#!/usr/bin/env bash
# e2e.sh - the end-to-end tests of the build in ~/samcombine/src on the mini database (docs/claude/2026-10-07-sam-combine),
# on 4 cores. Missing prerequisites (numpy, art_illumina) skip their tests unless REQUIRED=1.
set -u
S=${S:-$HOME/samcombine/src}
cd "$S" || exit 1
PROTAL_TEST_DB=$HOME/samcombine/db/database.protal PROTAL=${PROTAL_BIN:-$S/build/protal} SIMULATE=$S/build/simulate_metagenomes \
  PROTAL_TESTS_REQUIRED=${REQUIRED:-0} taskset -c 0-3 nice -n 5 python3 -m unittest ${TESTS:-tests/e2e/test_protal_e2e.py} 2>&1 | tail -${LINES_SHOWN:-25}
