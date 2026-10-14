#!/bin/bash
# The baseline: ab6bf25 (git archive) built in ~/split-base, protal and simulate_metagenomes, on cores 0-3.
set -o pipefail
S=$(cd "$(dirname "$0")" && pwd)
B=~/split-base
if [ ! -f $B/build/build.ninja ]; then
    rm -rf $B
    mkdir -p $B/src
    tar -xf $S/base.tar -C $B/src
    cmake -S $B/src -B $B/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON > $B/configure.log 2>&1 || { tail -30 $B/configure.log; exit 1; }
fi
cd $B/build
start=$(date +%s.%N)
taskset -c 0-3 nice -n 5 ninja -j4 ${@:-protal simulate_metagenomes} > $B/build.log 2>&1
status=$?
end=$(date +%s.%N)
echo "ninja exit $status, $(echo "$end - $start" | bc) s"
tail -3 $B/build.log
uptime
