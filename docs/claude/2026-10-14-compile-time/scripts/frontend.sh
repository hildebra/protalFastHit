#!/bin/bash
# Front-end cost (-fsyntax-only) of the main headers on their own: what a translation unit pays to include them.
W=~/ctime
B=$W/build
F=$W/frontend
rm -rf $F
mkdir -p $F
cd $F
CMD=$(grep -o '"command": "[^"]*protal.dir/src/main.cpp[^"]*"' $B/compile_commands.json | head -1 | sed 's/^"command": "//; s/"$//')
BASE=${CMD% -o *}
mk() { echo "#include \"$2\"" > $1.cpp; }
mk options Options.h
mk classify Classify.h
mk profiler Profiler.h
mk build Build.h
mk genomeloader GenomeLoader.h
mk strain Strain.h
mk seedmap Seedmap.h
mk zstd Zstd.h
mk runprotal RunProtal.h
cat > stdlibs.cpp <<'EOF'
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <random>
#include <ranges>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>
#include <immintrin.h>
#include <omp.h>
#include <cxxopts.hpp>
EOF
fe() {
    local u=$1
    local t0=$(date +%s.%N)
    taskset -c 0-3 nice -n 5 $(cat base_cmd.txt) -fsyntax-only $u.cpp 2> $u.err
    local t1=$(date +%s.%N)
    local lines=$($(cat base_cmd.txt) -E $u.cpp 2>/dev/null | wc -l)
    printf '%-14s %6.1f s  %7d preprocessed lines  %s\n' $u $(echo "$t1 - $t0" | bc) $lines "$(head -c 200 $u.err)"
}
echo "$BASE" > base_cmd.txt
export -f fe
uptime
printf '%s\n' stdlibs options classify profiler build genomeloader strain seedmap zstd runprotal | xargs -P4 -I{} bash -c 'fe {}'
uptime
