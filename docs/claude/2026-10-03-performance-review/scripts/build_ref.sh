#!/usr/bin/env bash
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/mt-work/perf4; mkdir -p $W
rm -rf $W/ref; mkdir -p $W/ref
git -C $REPO archive HEAD | tar -x -C $W/ref
git -C $REPO rev-parse --short HEAD > $W/ref/COMMIT
cd $W/ref && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF > $W/ref.configure.log 2>&1 && \
  nice cmake --build build --target protal -j 5 > $W/ref.build.log 2>&1 && echo "OK build $(cat $W/ref/COMMIT)" || echo "FAIL build"
