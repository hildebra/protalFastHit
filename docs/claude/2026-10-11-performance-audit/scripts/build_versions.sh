#!/usr/bin/env bash
# Builds protal 0.8.1 (b3ac271), 0.8.0 (a82df8f) and 0.7.8 (5de3224) from git archives in ~/perf9/<version>,
# Release, on 4 pinned cores (taskset -c 0-3, nice 5), one after another. 0.8.1 also builds simulate_metagenomes.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad/archives
W=$HOME/perf9
mkdir -p "$W"
pin=(taskset -c 0-3 nice -n 5)
launcher=()
command -v ccache > /dev/null && launcher=(-DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache)
for spec in v081:b3ac271:"protal simulate_metagenomes" v080:a82df8f:protal v078:5de3224:protal; do
    IFS=: read -r name commit targets <<< "$spec"
    E=$W/$name
    rm -rf "$E"; mkdir -p "$E/src"
    tar -xf "$A/$commit.tar" -C "$E/src"
    echo "$commit" > "$E/commit"
    start=$(date +%s)
    if "${pin[@]}" cmake -S "$E/src" -B "$E/build" -G Ninja -DCMAKE_BUILD_TYPE=Release "${launcher[@]}" > "$E/configure.log" 2>&1 &&
       "${pin[@]}" cmake --build "$E/build" -j 4 --target $targets > "$E/build.log" 2>&1; then
        echo "OK $name $commit in $(( $(date +%s) - start )) s: $("$E/build/protal" --version 2>&1 | head -1)"
    else
        echo "FAIL $name $commit"; grep -E 'error' -A3 "$E/build.log" | head -30
    fi
done
echo "ALL DONE"
