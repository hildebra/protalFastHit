#!/usr/bin/env bash
# unit.sh - the unit tests of the build in ~/samcombine/src (docs/claude/2026-10-07-sam-combine): the new ones first,
# then the whole suite as CI runs it, on 4 cores.
set -u
B=${B:-$HOME/samcombine/src/build}
cd "$B" || exit 1
taskset -c 0-3 nice -n 5 ./tests/protal_tests --gtest_filter='Options.*:LineSplitter.*' 2>&1 | tail -15
PROTAL_TESTS_REQUIRED=1 taskset -c 0-3 nice -n 5 ctest -j4 --output-on-failure 2>&1 | tail -8
