#!/usr/bin/env bash
# ThreadSanitizer build of the unit tests (clean tree + PATCH) in ~/mt-work/tsan; runs the profiling tests.
set -uo pipefail
PATCH=$1
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/mt-work/tsan
rm -rf $W; mkdir -p $W/src
git -C $REPO archive HEAD | tar -x -C $W/src
(cd $W/src && patch -p1 -s < $PATCH) || { echo "FAIL patch"; exit 1; }
cd $W/src
F="-fsanitize=thread -g -O1"
setarch "$(uname -m)" -R nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS="$F" -DCMAKE_C_FLAGS="$F" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread" > $W/configure.log 2>&1 || { echo "FAIL configure"; tail $W/configure.log; exit 1; }
setarch "$(uname -m)" -R nice cmake --build build --target protal_tests -j 4 > $W/build.log 2>&1 || { echo "FAIL build"; grep -E "error" -A3 $W/build.log | head -30; exit 1; }
cd build/tests
TSAN_OPTIONS="halt_on_error=0 second_deadlock_stack=1" setarch "$(uname -m)" -R ./protal_tests --gtest_filter='SamChunks.*:ProfileSam.*:MicrobialProfile.*:FromSam.*' > $W/tsan.log 2>&1
echo "exit $?; $(grep -c 'WARNING: ThreadSanitizer' $W/tsan.log) TSan warnings; $(grep -E 'PASSED|FAILED' $W/tsan.log | tr '\n' ' ')"
