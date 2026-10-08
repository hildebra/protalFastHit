#!/usr/bin/env bash
# build.sh ref|work: ~/avx512/ref = git archive HEAD; ~/avx512/work = the same with this session's changed files
# copied over it (files.txt), so other sessions' uncommitted edits stay out. Release with the unit tests, 4 cores.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
S=$(cd "$(dirname "$0")" && pwd)
W=$HOME/avx512; N=$1; E=$W/$N
mkdir -p $W
if [ ! -d $E/build ]; then
  rm -rf $E; mkdir -p $E
  git -C $REPO archive HEAD | tar -x -C $E
  git -C $REPO rev-parse --short HEAD > $W/$N.commit
fi
if [ "$N" = work ]; then
  while read -r f; do
    [ -z "$f" ] && continue
    mkdir -p "$(dirname "$E/$f")"
    cmp -s "$REPO/$f" "$E/$f" || { cp "$REPO/$f" "$E/$f"; touch "$E/$f"; }
  done < $S/files.txt
fi
cd $E
[ -d build ] || taskset -c 0-3 nice -n 5 cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/$N.configure.log 2>&1
taskset -c 0-3 nice -n 5 cmake --build build -j 4 > $W/$N.build.log 2>&1 && echo "OK build $N ($(cat $W/$N.commit))" || { echo "FAIL build $N"; grep -E 'error|Error' -A6 $W/$N.build.log | head -80; }
