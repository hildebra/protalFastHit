#!/bin/bash
# Runs the split build's unit tests (protal_tests) on cores 0-3.
cd ~/split-build/build
taskset -c 0-3 nice -n 5 ./tests/protal_tests --gtest_brief=1 > ~/split-build/unit.log 2>&1
status=$?
echo "protal_tests exit $status"
grep -E '^\[  (PASSED|FAILED)|FAILED|tests ran' ~/split-build/unit.log | head -30
