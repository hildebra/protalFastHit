#!/usr/bin/env bash
# Build the working tree's protal in WSL (~/spikecount/work), run the unit tests, build the mini database and run the
# e2e tests; everything on cores 0-3, niced. Logs in ~/spikecount/*.log, a summary at the end.
# Usage: build_test.sh [build|unit|e2e|all] (default all)
set -uo pipefail
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal
ROOT=$HOME/spikecount
W=$ROOT/work
what=${1:-all}
mkdir -p "$W"
run() { taskset -c 0-3 nice -n 5 "$@"; }

if [ "$what" = all ] || [ "$what" = build ]; then
    rsync -a --checksum --no-times --delete --exclude='/data/' --exclude='/build*/' --exclude='/cmake-build-*/' \
        --exclude='/local/' --exclude='/misc/' --exclude='/.claude/' --exclude='/.git/' --exclude='/strain_test_out/' \
        "$SRC/" "$W/"
    cd "$W" || exit 1
    start=$(date +%s)
    run cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON \
        -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache > "$ROOT/cmake.log" 2>&1 \
        || { tail -30 "$ROOT/cmake.log"; exit 1; }
    run ninja -C build -j4 protal simulate_metagenomes protal_tests > "$ROOT/ninja.log" 2>&1 \
        || { grep -n -A5 "error" "$ROOT/ninja.log" | head -80; exit 1; }
    echo "build: $(( $(date +%s) - start )) s"
fi
cd "$W" || exit 1
if [ "$what" = all ] || [ "$what" = unit ]; then
    start=$(date +%s)
    run ctest --test-dir build -j4 --output-on-failure > "$ROOT/unit.log" 2>&1
    echo "unit: exit $?, $(( $(date +%s) - start )) s; $(grep -E 'tests passed|tests failed' "$ROOT/unit.log" | tail -1)"
    grep -E "SpikeGenes|AnObserver|\*\*\*Failed|Failed  " "$ROOT/unit.log" | head -20
fi
if [ "$what" = all ] || [ "$what" = e2e ]; then
    start=$(date +%s)
    if [ ! -e data/mini_db/protal_db/database.protal ]; then
        PROTAL=build/protal run bash scripts/mini_db/build_mini_db.sh data/mini_db > "$ROOT/minidb.log" 2>&1 \
            || { tail -30 "$ROOT/minidb.log"; exit 1; }
    fi
    PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes \
        run python3 -m unittest ${E2E_TESTS:-tests/e2e/test_protal_e2e.py} > "$ROOT/e2e.log" 2>&1
    echo "e2e: exit $?, $(( $(date +%s) - start )) s"
    tail -5 "$ROOT/e2e.log"
    grep -E "^(FAIL|ERROR):" "$ROOT/e2e.log" | head -20
fi
