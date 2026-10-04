#!/usr/bin/env bash
# Builds the working tree (rsync of the checkout) in ~/perf-gtdb/work with tests; runs the tests given as arguments (gtest filter) or all.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/perf-gtdb; E=$W/work
mkdir -p $E
rsync -a --checksum --no-times --delete --exclude '/.git' --exclude '/build*/' --exclude '/data' --exclude '/docs' --exclude '/local' $REPO/ $E/
cd $E && { [ -d build ] || nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/work.configure.log 2>&1; } && \
  nice cmake --build build --target protal protal_tests -j 6 > $W/work.build.log 2>&1 && echo "OK build" || { echo FAIL; grep -E 'error|Error' -A3 $W/work.build.log | head -60; exit 1; }
if [ $# -gt 0 ]; then $E/build/tests/protal_tests --gtest_filter="$1" 2>&1 | tail -${2:-15}; else $E/build/tests/protal_tests 2>&1 | tail -12; fi
