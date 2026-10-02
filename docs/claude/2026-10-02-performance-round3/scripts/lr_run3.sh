#!/usr/bin/env bash
# The prototype with the middle-first order: short reads must stay byte-identical (100k pairs, -t 1); long-read
# records against the reference; callgrind of both re-seeding variants (TAG=3).
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/perf3/build_exp.sh lr || exit 1
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points; I=$HOME/mt-work/isacg
LR=$HOME/mt-work/perf3-lr/src/build/protal
for b in ref lr; do bin=$W/ref/build/protal; [ $b = lr ] && bin=$LR
  rm -rf $W/pe.$b; nice $bin --db $DB -1 $I/r1.fq -2 $I/r2.fq --no_profile --prefix s -o $W/pe.$b -t 1 --no_qcmsa > $W/pe.$b.log 2>&1 &
done; wait
cmp -s <(zstdcat $W/pe.ref/s.sam.zst) <(zstdcat $W/pe.lr/s.sam.zst) && echo "short reads: SAM text same" || echo "short reads: SAM DIFFER"
for s in ont_b3000000 pb_b3000000; do t=ont; [ $s = pb_b3000000 ] && t=pb
  rm -rf $W/o.$s.lr3; nice $LR --db $DB -1 $P/$s/sim/reads/${s}_s_1.fq.gz --read_type $t --no_profile --prefix s -o $W/o.$s.lr3 -t 1 --no_qcmsa > /dev/null 2>&1
  echo "== $s reseed, middle first, vs ref"; bash $S/perf3/lrcmp.sh $W/o.$s.ref/s.sam.zst $W/o.$s.lr3/s.sam.zst
  cmp -s <(zstdcat $W/o.$s.reseed/s.sam.zst) <(zstdcat $W/o.$s.lr3/s.sam.zst) && echo "   same SAM as before the reorder" || echo "   SAM differs from before the reorder"
done
TAG=3 bash $S/perf3/lr_cg2.sh
