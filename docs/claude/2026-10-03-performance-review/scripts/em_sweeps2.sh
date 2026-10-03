#!/usr/bin/env bash
# Item 2: the EM's convergence sweep by sweep, printed by a copy of the p4 tree (sweeps_insert.cpp inserted before
# weight.swap(next) in AbundanceWeightedShares), on the 500k, 5M, ONT and PacBio SAMs.
set -uo pipefail
W=$HOME/mt-work/perf4; T=$W/sweeps3; DB=$HOME/bench071/V073/protal_db
SCR=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/85c3e3d5-6df2-4cf2-a5bc-63aafb61f818/scratchpad
rm -rf $T; mkdir -p $T; (cd $W/p4 && tar -c --exclude=./build .) | tar -x -C $T
cd $T
F=src/Profiling/SampleContext.h
awk -v f=$SCR/sweeps_insert.cpp 'BEGIN { while ((getline l < f) > 0) ins = ins l "\n" } /weight\.swap\(next\);/ && !d { printf "%s", ins; d = 1 } { print }' $F > $F.new && mv $F.new $F
grep -q '#include <iostream>' $F || sed -i 's|#include <algorithm>|#include <algorithm>\n#include <iostream>|' $F
echo "inserted: $(grep -c '\[sweep\]' $F)"
nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF > $W/sweeps3.configure.log 2>&1 && nice cmake --build build --target protal -j 5 > $W/sweeps3.build.log 2>&1 && echo "OK sweeps3 build" || { echo "FAIL sweeps3 build"; grep -E 'error' -A3 $W/sweeps3.build.log | head -30; exit 1; }
for spec in "500k $W/o.pe500k_t1/s.sam.zst" "5M $W/o.pe5M_t6/s.sam.zst" "ont90M $W/o.ont90M_t6/s.sam.zst" "pb90M $W/o.pb90M_t6/s.sam.zst"; do
  set -- $spec; n=$1; sam=$2; rm -rf $W/sw.$n
  $T/build/protal --db $DB --profile_only $sam --prefix s -o $W/sw.$n -t 6 --no_qcmsa > $W/sw.$n.log 2> $W/sw.$n.err
  echo "== $n: $(grep -c '^\[sweep\]' $W/sw.$n.err) sweeps"
  grep '^\[sweep\]' $W/sw.$n.err | awk '$2==1||$2==2||$2==3||$2==5||$2==10||$2==20||$2==30||$2==50||$2==75||$2==100||$2==150||$2==200'
  grep '^\[sweep\]' $W/sw.$n.err | tail -1
  cp $W/sw.$n.err $W/sweeps_$n.txt
done
