#!/bin/bash
# Debug + ASan/UBSan build of protal_tests with -DPROTAL_NO_CLONES (as the CI sanitizer job), then ctest, on CPUs 0-3.
# Usage: asan.sh [ctest -R regex]
set -o pipefail
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal-testaudit/
DST=~/testaudit/src/
rsync -a --checksum --no-times --delete --exclude=/build*/ --exclude=/data --exclude=.git --exclude=__pycache__ "$SRC" "$DST"
cd "$DST" || exit 1
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
start=$(date +%s)
if [ ! -f build-asan/build.ninja ]; then
  taskset -c 0-3 cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
    -DCMAKE_CXX_FLAGS="$FLAGS -DPROTAL_NO_CLONES" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS" > ~/ti_asan_cmake.log 2>&1 || { tail -30 ~/ti_asan_cmake.log; exit 1; }
fi
taskset -c 0-3 nice -n 5 cmake --build build-asan --target protal_tests -- -j4 > ~/ti_asan_build.log 2>&1
rc=$?
echo "build exit $rc after $(( $(date +%s) - start )) s"
[ $rc -ne 0 ] && { grep -E 'error|Error' ~/ti_asan_build.log | head -40; exit $rc; }
start=$(date +%s)
args=(-j4 --output-on-failure --timeout 900)
[ -n "$1" ] && args+=(-R "$1")
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 taskset -c 0-3 nice -n 5 ctest --test-dir build-asan "${args[@]}" > ~/ti_asan_ctest.log 2>&1
rc=$?
echo "ctest exit $rc after $(( $(date +%s) - start )) s"
grep -E 'tests passed|tests failed|\(Failed\)|\(Timeout\)|\(SEGFAULT\)|\(Subprocess|ERROR: AddressSanitizer|runtime error' ~/ti_asan_ctest.log | head -40
grep -E 'Test +#[0-9]+: ' ~/ti_asan_ctest.log | sed -E 's/.*Test +#[0-9]+: ([^ ]+) .* ([0-9.]+) sec/\2 \1/' | sort -rn | head -12
exit $rc
