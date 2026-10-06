#!/usr/bin/env bash
# build_proto.sh <name> <patch>: ~/perf6/<name> = HEAD's tree (~/perf6/head, from build_head.sh) with the patch applied,
# protal built (Release). The patches of this report: flex_ties_avx2.patch (name "seed"), screen_packed.patch ("screen").
set -uo pipefail
N=$1; PATCH=$(realpath "$2"); W=$HOME/perf6; E=$W/$N
rm -rf $E; mkdir -p $E; cp -a $W/head/. $E/; rm -rf $E/build
(cd $E && patch -p1 < "$PATCH") || { echo "patch failed"; exit 1; }
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $W/$N.configure.log 2>&1 && \
  nice cmake --build build --target protal -j 5 > $W/$N.build.log 2>&1 && echo "OK build $N" || { echo FAIL; grep -E 'error' -A3 $W/$N.build.log | head -30; }
