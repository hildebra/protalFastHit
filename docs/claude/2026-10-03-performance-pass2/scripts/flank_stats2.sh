#!/usr/bin/env bash
# Re-run of the flank statistics with `score == 4 * mismatches` (the ungapped alignment was optimal) instead of <=.
set -uo pipefail
W=$HOME/perf-pass2; E=$W/exp; DB=$HOME/bench071/V073/protal_db
sed -i 's/if (sc <= 4\*flank_mm) h\[flank_mm\*4+1\]++;/if (sc == 4*flank_mm) h[flank_mm*4+1]++; else if (sc < 4*flank_mm) h[flank_mm*4+2]++;/' $E/src/Alignment/AnchoredAlignment.h
grep -c "sc == 4" $E/src/Alignment/AnchoredAlignment.h
cd $E && nice cmake --build build --target protal -j 6 > $W/exp.build2.log 2>&1 || { grep error -A3 $W/exp.build2.log | head -30; exit 1; }
R=$W/reads; rm -rf $W/oexp
$E/build/protal --db $DB -1 $R/pe100k_R1.fq.gz -2 $R/pe100k_R2.fq.gz --read_type pe --no_profile --prefix s -o $W/oexp -t 1 --no_qcmsa 2>&1 | grep FLANKSTAT | sed 's/ownscore<=/ungapped_optimal=/; s/fail=/cheaper_than_ungapped=/' | head -14
