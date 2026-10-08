#!/usr/bin/env bash
# The reader tests under ThreadSanitizer (Debug, the work tree, hot functions compiled once: -DPROTAL_NO_CLONES, as
# TSan cannot run the ifunc resolvers of target_clones before it starts), 4 cores.
G=$HOME/gzpar; E=$G/work; T="taskset -c 0-3 nice -n 5"
FL="-fsanitize=thread -fno-omit-frame-pointer -g -O1 -DPROTAL_NO_CLONES"
cd $E
rm -rf build-tsan
setarch "$(uname -m)" -R $T cmake -S . -B build-tsan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
   -DCMAKE_CXX_FLAGS="$FL" -DCMAKE_C_FLAGS="-fsanitize=thread -g -O1" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread" > $G/tsan.configure.log 2>&1
setarch "$(uname -m)" -R $T cmake --build build-tsan --target protal_tests -j 4 > $G/tsan.build.log 2>&1 || { echo "TSAN BUILD FAIL"; grep -E "error|Result" -A3 $G/tsan.build.log | head -20; }
[ -x build-tsan/tests/protal_tests ] || exit 1
cd build-tsan/tests
TSAN_OPTIONS="halt_on_error=0 second_deadlock_stack=1" setarch "$(uname -m)" -R $T ./protal_tests --gtest_filter='ThreadedGzStream.*' > $G/tsan.log 2>&1
echo "exit $?"
grep -E "^\[  (PASSED|FAILED)|WARNING: ThreadSanitizer|SUMMARY" $G/tsan.log | sort | uniq -c | head -20
echo TSAN_DONE
