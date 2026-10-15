#!/usr/bin/env bash
# Build protal (Release) from ~/perf085/src (HEAD 4ec3f65) on cores 0-3; the binary is kept as bin/protal_base.
set -uo pipefail
W=~/perf085
cd $W/src
if [ ! -f $W/build/build.ninja ]; then
  cmake -S . -B $W/build -G Ninja -DCMAKE_BUILD_TYPE=Release > $W/configure.log 2>&1 || { tail -30 $W/configure.log; exit 1; }
fi
taskset -c 0-3 nice -n 5 ninja -C $W/build -j4 protal > $W/ninja.log 2>&1
echo "ninja exit $?"
tail -3 $W/ninja.log
cp $W/build/protal $W/bin/protal_base && ls -la $W/bin
echo BUILD_DONE
