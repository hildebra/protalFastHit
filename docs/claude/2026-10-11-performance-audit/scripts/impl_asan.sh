#!/usr/bin/env bash
# After impl_final.sh: the unit tests under ASan + UBSan as CI's sanitizer job (Debug, PROTAL_NO_CLONES), from the
# work tree as impl_final.sh left it; the system's libraries (no micromamba on PATH). 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
until grep -q "POST DONE" "$A/impl_post.out" 2> /dev/null; do sleep 15; done
W=$HOME/impl; S=$W/work/src
pin=(taskset -c 0-3 nice -n 5)
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
rm -rf "$W/asan"
start=$(date +%s)
echo "[$(date +%T)] ASan build"
if "${pin[@]}" cmake -S "$S" -B "$W/asan" -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
       -DCMAKE_CXX_FLAGS="$FLAGS -DPROTAL_NO_CLONES" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS" > "$W/asan_configure.log" 2>&1 &&
   "${pin[@]}" cmake --build "$W/asan" -j 4 --target protal_tests > "$W/asan_build.log" 2>&1; then
  echo "asan build in $(( $(date +%s) - start )) s"
  start=$(date +%s)
  ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 "${pin[@]}" ctest --test-dir "$W/asan" --output-on-failure -j4 > "$W/asan_tests.log" 2>&1
  echo "asan ctest exit $? in $(( $(date +%s) - start )) s"
  tail -6 "$W/asan_tests.log"
  grep -a -E "ERROR: (Address|Leak)Sanitizer|runtime error|\*\*\*Failed|FAILED" "$W/asan_tests.log" | head -20
else
  echo "asan build FAIL"; tail -12 "$W/asan_configure.log"; grep -a -E "error" -A4 "$W/asan_build.log" | head -40
fi
echo "[$(date +%T)] ASAN DONE"
