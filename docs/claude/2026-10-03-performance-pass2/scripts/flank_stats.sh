#!/usr/bin/env bash
# Experiment build: counts the flanks AnchoredAligner::Flank sends to WFA2 (by mismatches on the link's diagonal and
# length) and the resulting scores; runs 100k pairs at one thread.
set -uo pipefail
W=$HOME/perf-pass2; E=$W/exp; DB=$HOME/bench071/V073/protal_db
rm -rf $E; mkdir -p $E; cp -a $W/head/. $E/
F=$E/src/Alignment/AnchoredAlignment.h
perl -0pi -e 's|(            aligner\.Reset\(\);\n            aligner\.Alignment\(m_query, m_ref, 0, read_end_free, 0, ref_end_free, max_score - used\);\n            if \(!aligner\.Success\(\)\) return m_status = Status::Failed;)|            { static std::atomic<long> tot{0}, hist[64][4]; static struct P { ~P() { for (int i = 0; i < 64; i++) fprintf(stderr, "FLANKSTAT mm=%d ok=%ld ownscore<=%ld fail=%ld ungappedshape=%ld\\n", i, hist[i][0].load(), hist[i][1].load(), hist[i][2].load(), hist[i][3].load()); } } p; (void)p;\n              int mm = 0; for (size_t k = 0; k < m; k++) mm += m_query[k] != m_ref[k]; if (mm > 63) mm = 63; flank_mm = mm; flank_hist = &hist[0][0]; (void)tot; }\n$1|' $F
perl -0pi -e 's|(            used \+= -aligner\.GetAlignmentScore\(\);\n            aligner\.CigarInto\(out\);\n            return m_status = Status::Aligned;)|            { auto* h = reinterpret_cast<std::atomic<long>*>(flank_hist); int sc = -aligner.GetAlignmentScore(); h[flank_mm*4+0]++; if (sc <= 4*flank_mm) h[flank_mm*4+1]++; if (m_ungapped_flanks \&\& r - m <= (size_t)read_end_free \&\& g - m <= (size_t)ref_end_free) h[flank_mm*4+3]++; }\n$1|' $F
perl -0pi -e 's|(        std::string m_query, m_ref, m_left, m_right;)|        int flank_mm = 0; void* flank_hist = nullptr;\n$1|; s|#include <algorithm>|#include <algorithm>\n#include <atomic>\n#include <cstdio>|' $F
perl -0pi -e 's|(            if \(!aligner\.Success\(\)\) return m_status = Status::Failed;\n            used \+= -aligner)|            if (!aligner.Success()) { reinterpret_cast<std::atomic<long>*>(flank_hist)[flank_mm*4+2]++; return m_status = Status::Failed; }\n            used += -aligner|' $F
cd $E && rm -rf build && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF > $W/exp.configure.log 2>&1 && nice cmake --build build --target protal -j 6 > $W/exp.build.log 2>&1 || { grep error -A3 $W/exp.build.log | head -30; exit 1; }
R=$W/reads; rm -rf $W/oexp
$E/build/protal --db $DB -1 $R/pe100k_R1.fq.gz -2 $R/pe100k_R2.fq.gz --read_type pe --no_profile --prefix s -o $W/oexp -t 1 --no_qcmsa 2>&1 | grep FLANKSTAT | awk '{print}' | head -20
