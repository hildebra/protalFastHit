#!/usr/bin/env bash
# Builds HEAD (git archive, not the shared working tree) in ~/perf6/head: protal and the unit tests.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/perf6; E=$W/head
mkdir -p $W; rm -rf $E; mkdir -p $E
git -C $REPO rev-parse --short HEAD > $W/head.commit
git -C $REPO archive HEAD | tar -x -C $E
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $W/head.configure.log 2>&1 && \
  nice cmake --build build -j 6 > $W/head.build.log 2>&1 && echo "OK build $(cat $W/head.commit)" || { echo FAIL; grep -E 'error' -A3 $W/head.build.log | head -40; }
ls -la $E/build/protal $E/build/protal_tests 2>/dev/null
