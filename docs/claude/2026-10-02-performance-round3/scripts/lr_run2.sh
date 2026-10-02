#!/usr/bin/env bash
# Reference, prototype without re-seeding (PROTAL_LR_NORESEED=1) and with it, on the Nanopore and PacBio 3 Mb
# samples, one thread, SAM only, alternated twice; records compared with the reference; callgrind of the re-seeded.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/perf3/build_exp.sh lr || exit 1
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
LR=$HOME/mt-work/perf3-lr/src/build/protal
run() { local s=$1 v=$2 t=ont bin=$LR env=""; [ $s = pb_b3000000 ] && t=pb; [ $v = ref ] && bin=$W/ref/build/protal; [ $v = noreseed ] && env="PROTAL_LR_NORESEED=1"
  rm -rf $W/o.$s.$v
  env $env nice $bin --db $DB -1 $P/$s/sim/reads/${s}_s_1.fq.gz --read_type $t --no_profile --prefix s -o $W/o.$s.$v -t 1 --no_qcmsa --verbose > $W/o.$s.$v.log 2>&1 || echo "FAIL $s $v"
  echo "$s $v: $(grep -E 'Alignment handler took' $W/o.$s.$v.log) | $(grep 'Total alignments' $W/o.$s.$v.log | tr -s '\t ' ' ')"; }
for rep in 1; do for s in ont_b3000000 pb_b3000000; do for v in ref noreseed reseed; do run $s $v; done; done; done
for s in ont_b3000000 pb_b3000000; do for v in noreseed reseed; do echo "== $s $v vs ref"; bash $S/perf3/lrcmp.sh $W/o.$s.ref/s.sam.zst $W/o.$s.$v/s.sam.zst; done; done
bash /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad/perf3/lr_cg2.sh
