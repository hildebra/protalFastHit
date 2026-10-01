#!/bin/bash
# Experiment build: counts, in AnchoredAligner::Flank, how many WFA calls a flank would not need because its
# ungapped alignment has at most one mismatch (4 against at least 8 for any alignment with a gap), by
# mismatch count. Totals are printed to stderr at exit; run single-threaded. The repository is not changed.
set -e
source "$(dirname "$0")/env.sh"
rm -rf $PERF_DIR/src-fc; mkdir -p $PERF_DIR/src-fc
rsync -a --exclude '/data' --exclude '/build*/' --exclude '/.git' $PERF_DIR/src/ $PERF_DIR/src-fc/
cd $PERF_DIR/src-fc/src
cat > Alignment/FlankCounts.h <<'HDR'
#pragma once
#include <algorithm>
#include <cstdio>
namespace protal {
    struct FlankCounts {
        unsigned long calls = 0, fit = 0, m[6] = {}, lengths = 0;
        ~FlankCounts() {
            std::fprintf(stderr, "FLANKS calls=%lu fit=%lu mismatches(ungapped) 0:%lu 1:%lu 2:%lu 3:%lu 4:%lu 5+:%lu mean_len=%.1f\n",
                         calls, fit, m[0], m[1], m[2], m[3], m[4], m[5], calls ? double(lengths) / calls : 0.0);
        }
    };
    inline FlankCounts g_flank_counts;
}
HDR
perl -0pi -e 's/#include <cstddef>/#include <cstddef>\n#include "FlankCounts.h"/' Alignment/AnchoredAlignment.h
perl -0pi -e 's/(            aligner.Reset\(\);\n            aligner.Alignment\(m_query, m_ref, 0, std::min<int>\(std::max\(read_free, 0\), static_cast<int>\(r\)\),)/            { g_flank_counts.calls++; g_flank_counts.lengths += r; if (r <= g) { g_flank_counts.fit++; size_t mm = 0; for (size_t i = 0; i < r; i++) mm += m_query[i] != m_ref[i]; g_flank_counts.m[std::min<size_t>(mm, 5)]++; } }\n$1/' Alignment/AnchoredAlignment.h
grep -n "g_flank_counts" Alignment/AnchoredAlignment.h | cut -c1-80
cd $PERF_DIR
cmake -S src-fc -B build-fc -G Ninja -DCMAKE_BUILD_TYPE=Release > build-fc.cmake.log 2>&1
cmake --build build-fc --target protal_avx2 -j "$(nproc)" > build-fc.log 2>&1
ls -la build-fc/protal_avx2
