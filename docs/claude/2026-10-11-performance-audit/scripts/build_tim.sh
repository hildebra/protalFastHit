#!/usr/bin/env bash
# The audit's experiment build: 0.8.1 (b3ac271) with per-table timers and peak resets in the sequential start-up
# (scratchpad tim/RunProtal.h), in ~/perf9/tim. Pinned to the cores the other work uses (taskset -c 0-3).
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
E=$HOME/perf9/tim
mkdir -p "$E/src"; [ -f "$E/src/CMakeLists.txt" ] || tar -xf "$A/archives/b3ac271.tar" -C "$E/src"

cp "$A/tim/RunProtal.h" "$E/src/src/RunProtal.h"; cp "$A/tim/Profiler.h" "$E/src/src/Profiling/Profiler.h"
{ diff <(tar -xOf "$A/archives/b3ac271.tar" src/RunProtal.h) "$E/src/src/RunProtal.h"; diff <(tar -xOf "$A/archives/b3ac271.tar" src/Profiling/Profiler.h) "$E/src/src/Profiling/Profiler.h"; } > "$E/audit.diff"
launcher=()
command -v ccache > /dev/null && launcher=(-DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache)
start=$(date +%s)
if taskset -c 0-3 nice -n 5 cmake -S "$E/src" -B "$E/build" -G Ninja -DCMAKE_BUILD_TYPE=Release "${launcher[@]}" > "$E/configure.log" 2>&1 &&
   taskset -c 0-3 nice -n 5 cmake --build "$E/build" -j 4 --target protal > "$E/build.log" 2>&1; then
    echo "OK tim in $(( $(date +%s) - start )) s; diff lines $(wc -l < "$E/audit.diff")"
else
    echo "FAIL tim"; grep -E 'error' -A3 "$E/build.log" | head -30
fi
