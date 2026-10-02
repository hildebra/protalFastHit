#!/usr/bin/env bash
# The clean long-read change (lr2.patch on f343113): check.sh (unit, mini DB, e2e), the new tests shown; short reads
# byte-identical to f343113 (100k pairs -t 1); long-read records against the reference and callgrind as before.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/check.sh lr2 $S/perf3/lr2.patch f343113 || { grep -E "Failure|FAILED" $HOME/mt-work/lr2/ctest.log | head; }
B=$HOME/mt-work/lr2/src/build
(cd $B && ./tests/protal_tests --gtest_filter='AnchoredAlignment.*' 2>&1 | grep -E "long reads:|^\[ +(OK|FAILED|PASSED) " | head -14)
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points; I=$HOME/mt-work/isacg
ID=$HOME/mt-work/perf3-id/src/build/protal; LR=$B/protal
rm -rf $W/l2.pe.id $W/l2.pe.lr
nice $ID --db $DB -1 $I/r1.fq -2 $I/r2.fq --no_profile --prefix s -o $W/l2.pe.id -t 1 --no_qcmsa > /dev/null 2>&1
nice $LR --db $DB -1 $I/r1.fq -2 $I/r2.fq --no_profile --prefix s -o $W/l2.pe.lr -t 1 --no_qcmsa > /dev/null 2>&1
cmp -s <(zstdcat $W/l2.pe.id/s.sam.zst) <(zstdcat $W/l2.pe.lr/s.sam.zst) && echo "short reads 100k -t 1: SAM text same as f343113" || echo "short reads: SAM DIFFER"
declare -A RD=( [ont_b3000000]="$P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz ont" [pb_b3000000]="$P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz pb" [hifi3M]="$W/hifi/hifi_3M.fq.gz pb" )
for s in ont_b3000000 pb_b3000000 hifi3M; do set -- ${RD[$s]}
  rm -rf $W/o.$s.lr2; nice $LR --db $DB -1 $1 --read_type $2 --no_profile --prefix s -o $W/o.$s.lr2 -t 1 --no_qcmsa > /dev/null 2>&1
  echo "== $s"; bash $S/perf3/lrcmp.sh $W/o.$s.ref/s.sam.zst $W/o.$s.lr2/s.sam.zst
  cmp -s <(zstdcat $W/o.$s.lr4/s.sam.zst) <(zstdcat $W/o.$s.lr2/s.sam.zst) && echo "   same SAM as the prototype" || echo "   SAM differs from the prototype"
done
echo LR2 DONE
