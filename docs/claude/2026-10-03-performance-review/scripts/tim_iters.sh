#!/usr/bin/env bash
# After tim_build.sh: the EM's iterations and classes printed, rebuilt, run on the 500k and 5M SAMs.
set -uo pipefail
W=$HOME/mt-work/perf4; T=$W/tim
while ! grep -q "Run protal took" $W/o.tim_ont90M_t6.log 2>/dev/null; do sleep 5; done
cd $T
perl -0pi -e 's|(\n\s*)weight = std::move\(next\);\n(\s*)if \(change < kEmTolerance\) break;|${1}weight = std::move(next);\n${2}if (change < kEmTolerance \|\| iteration + 1 == kEmIterations) { std::cerr << "[tim] em_iterations " << iteration + 1 << " change " << change << std::endl; }\n${2}if (change < kEmTolerance) break;|' src/Profiling/SampleContext.h
grep -q '#include <iostream>' src/Profiling/SampleContext.h || perl -0pi -e 's|#include <algorithm>|#include <algorithm>\n#include <iostream>|' src/Profiling/SampleContext.h
grep -c 'em_iterations' src/Profiling/SampleContext.h
nice cmake --build build --target protal -j 5 > $W/tim2.build.log 2>&1 && echo "OK tim2 build" || { echo "FAIL tim2 build"; grep -E 'error' -A3 $W/tim2.build.log | head -30; exit 1; }
B=$T/build/protal; DB=$HOME/bench071/V073/protal_db
for spec in "500k 1 $W/o.pe500k_t1/s.sam.zst" "5M 6 $W/o.pe5M_t6/s.sam.zst"; do
  set -- $spec; n=$1; t=$2; sam=$3; rm -rf $W/o.tim2_${n}_t$t
  $B --db $DB --profile_only $sam --prefix s -o $W/o.tim2_${n}_t$t -t $t --no_qcmsa --verbose > $W/o.tim2_${n}_t$t.log 2> $W/o.tim2_${n}_t$t.err
  echo "== tim2 $n t=$t"; grep -E '^\[tim\]' $W/o.tim2_${n}_t$t.err | grep -E 'em_|congener'
done
