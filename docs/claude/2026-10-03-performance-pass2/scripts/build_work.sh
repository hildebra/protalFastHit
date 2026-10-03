#!/usr/bin/env bash
# Builds the working tree (rsync of the checkout) in ~/perf-pass2/work with tests, runs the touched test suites.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/perf-pass2; E=$W/work
mkdir -p $E
rsync -a --checksum --no-times --delete --exclude '/.git' --exclude '/build*/' --exclude '/data' --exclude '/docs' $REPO/ $E/ 
cd $E && { [ -d build ] || nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/work.configure.log 2>&1; } && \
  nice cmake --build build --target protal protal_tests -j 6 > $W/work.build.log 2>&1 && echo "OK build" || { echo FAIL; grep -E 'error' -A3 $W/work.build.log | head -40; exit 1; }
$E/build/tests/protal_tests 2>&1 | tail -12
