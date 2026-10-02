#!/usr/bin/env bash
# callgrind of the prototype without and with re-seeding (PROTAL_LR_NORESEED), Nanopore and PacBio 3 Mb, one thread,
# after lr_run2.sh built it; alignment parts and WFA2's slab reset.
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
LR=$HOME/mt-work/perf3-lr/src/build/protal; TAG=${TAG:-2}
for s in ont_b3000000 pb_b3000000; do t=ont; [ $s = pb_b3000000 ] && t=pb
  for v in noreseed reseed; do env=""; [ $v = noreseed ] && env="PROTAL_LR_NORESEED=1"
    env $env nice valgrind --tool=callgrind --callgrind-out-file=$W/cg.$s.$v$TAG $LR --db $DB -1 $P/$s/sim/reads/${s}_s_1.fq.gz --read_type $t --no_profile --prefix s -o $W/cg.o.$s.$v$TAG -t 1 --no_qcmsa > /dev/null 2>&1 &
  done
done; wait
for s in ont_b3000000 pb_b3000000; do for v in noreseed reseed; do
  callgrind_annotate --inclusive=yes --threshold=99.5 $W/cg.$s.$v$TAG 2>/dev/null | c++filt > $W/cg.$s.$v$TAG.incl.txt
  echo "== $s $v$TAG: $(grep 'PROGRAM TOTALS' $W/cg.$s.$v$TAG.incl.txt | awk '{print $1}')"
  grep -E 'RunLongReads|AnchoredAligner::(Align|Flank|Gap|Reseed|Piece)|Between|alignEndsFree|alignEnd2End|slab_reap_repurpose' $W/cg.$s.$v$TAG.incl.txt | sed -E 's/\[\/home[^]]*\]//; s/\(([^()]|\([^()]*\))*\)//g' | awk '{ $1 = $1; print }' | cut -c1-120 | head -12
done; done
