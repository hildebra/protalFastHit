#!/usr/bin/env bash
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/check.sh rtype $S/perf3/rtype.patch f343113 || exit 1
B=$HOME/mt-work/rtype/src/build
(cd $B && ./tests/protal_tests --gtest_filter='ReadTypeDetection.*:Options.TheReadTypeOf*' 2>&1 | grep -E "^\[ +(OK|FAILED|PASSED) " | head)
W=$HOME/mt-work/perf3; L=$HOME/bench071/samples_lr072/points
for f in $L/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz $L/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz $W/hifi/hifi_3M.fq.gz $HOME/bench071/samples/points/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; do
  rm -rf $W/rt; timeout 15 $B/protal --db $HOME/bench071/V072/protal_db -1 $f --prefix s -o $W/rt -t 2 --no_qcmsa --no_profile > $W/rt.log 2>&1
  echo "== $(basename $f): $(grep -a -E '^Note: Sample|^Error' $W/rt.log | head -2 | sed 's/ in \/home.*//')"
done
echo RTYPE2 DONE
