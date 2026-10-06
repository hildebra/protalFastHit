#!/usr/bin/env bash
# HEAD (same archive as ~/perf6/head) built with -g on top of the Release flags: the same code, with line tables, so
# callgrind can attribute inlined code (GetGenome, Gene::Window, the screen) to its source lines.
set -uo pipefail
W=$HOME/perf6; E=$W/headg
rm -rf $E; mkdir -p $E
cp -a $W/head/. $E/ && rm -rf $E/build
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-g -DCMAKE_C_FLAGS=-g > $W/headg.configure.log 2>&1 && \
  nice cmake --build build --target protal -j 4 > $W/headg.build.log 2>&1 && echo "OK build -g" || { echo FAIL; grep -E 'error' -A3 $W/headg.build.log | head -40; }
