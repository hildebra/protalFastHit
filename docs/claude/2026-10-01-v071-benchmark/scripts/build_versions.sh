#!/usr/bin/env bash
# The three versions, each from git (git archive into $B/src/<version>, so no working-tree change is in it),
# built in Release (Ninja, gcc of the machine): 0.6.0a (the tag), 0.7.0 (4b21427) and 0.7.1 (1c11a00), the
# binaries copied to $B/bin/protal-<version>; simulate_metagenomes of 0.7.1 simulates every sample.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
REPO=${REPO:-/mnt/c/Users/hildebra/Documents/locDev/protal}
mkdir -p $B/bin $B/src $B/logs
for pair in "0.6.0a 0.6.0a" "0.7.0 4b21427" "0.7.1 1c11a00"; do
  set -- $pair
  v=$1 rev=$2 S=$B/src/$1
  [ -x $B/bin/protal-$v ] && continue
  rm -rf $S; mkdir -p $S
  git -C $REPO archive $rev | tar -x -C $S || exit 1
  targets="protal"
  [ $v = 0.7.1 ] && targets="protal simulate_metagenomes"
  (cd $S && cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $B/logs/configure_$v.log 2>&1 &&
    cmake --build build --target $targets -j 4 > $B/logs/build_$v.log 2>&1) || { echo "build of $v failed"; tail -20 $B/logs/build_$v.log; exit 1; }
  cp $S/build/protal $B/bin/protal-$v
  [ $v = 0.7.1 ] && cp $S/build/simulate_metagenomes $B/bin/simulate_metagenomes
  echo "$v ($rev): $($B/bin/protal-$v --version 2>&1 | tail -1)"
done
