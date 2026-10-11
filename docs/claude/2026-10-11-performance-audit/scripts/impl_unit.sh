#!/usr/bin/env bash
# Warnings from the changed files, then the unit tests (all, or FILTER), pinned to cores 0-3.
W=$HOME/impl/work
grep -E "warning" "$W/build.log" | grep -E "Database.h|BinaryTable.h|ColumnWeights.h|StrainAlleles.h|AncestrySites.h|RunProtal.h|Build.h|Options.h|Profiler.h" | sort -u | head -20
cd "$W/build/tests" 2> /dev/null || cd "$W/build"
bin=$(find "$W/build" -name protal_tests -type f | head -1)
taskset -c 0-3 nice -n 5 "$bin" --gtest_brief=1 ${FILTER:+--gtest_filter=$FILTER} 2>&1 | tail -25
