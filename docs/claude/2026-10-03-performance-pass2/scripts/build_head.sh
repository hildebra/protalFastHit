#!/usr/bin/env bash
# Builds HEAD (git archive) in ~/perf-pass2/head with tests.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/perf-pass2; E=$W/head
mkdir -p $W; rm -rf $E; mkdir -p $E
git -C $REPO archive HEAD | tar -x -C $E
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF > $W/head.configure.log 2>&1 && \
  nice cmake --build build --target protal -j 6 > $W/head.build.log 2>&1 && echo "OK build" || { echo FAIL; grep -E 'error' -A3 $W/head.build.log | head -40; }
git -C $REPO rev-parse --short HEAD
