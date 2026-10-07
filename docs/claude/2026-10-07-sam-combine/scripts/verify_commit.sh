#!/usr/bin/env bash
# verify_commit.sh TAR NAME - builds a `git archive` of a commit in ~/samcombine/NAME and runs the unit and end-to-end
# tests on the mini database (docs/claude/2026-10-07-sam-combine), on 4 cores. The checkout holds other sessions'
# uncommitted work, so a commit is checked from its own files only.
set -u
TAR=$1
NAME=$2
D=$HOME/samcombine/$NAME
rm -rf "$D"; mkdir -p "$D"
tar -xf "$TAR" -C "$D"
cd "$D" || exit 1
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON \
      -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache > configure.log 2>&1 || { echo "configure failed"; exit 1; }
taskset -c 0-3 nice -n 5 ninja -C build -j4 > build.log 2>&1 || { echo "build failed"; tail -30 build.log; exit 1; }
echo "built"
(cd build && PROTAL_TESTS_REQUIRED=1 taskset -c 0-3 nice -n 5 ctest -j4 --output-on-failure 2>&1 | tail -6)
PROTAL_TEST_DB=$HOME/samcombine/db/database.protal PROTAL=$D/build/protal SIMULATE=$D/build/simulate_metagenomes \
  taskset -c 0-3 nice -n 5 python3 -m unittest tests/e2e/test_protal_e2e.py 2>&1 | tail -4
echo "done"
