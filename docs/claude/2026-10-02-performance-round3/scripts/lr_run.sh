#!/usr/bin/env bash
# The long-read prototype (lr.patch) built, then both builds on the Nanopore and PacBio 3 Mb samples, one thread,
# alternated twice, SAM only; stage times and the records compared.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/perf3/build_exp.sh lr || exit 1
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
declare -A BIN=( [ref]=$W/ref/build/protal [lr]=$HOME/mt-work/perf3-lr/src/build/protal )
for rep in 1 2; do
  for s in ont_b3000000 pb_b3000000; do
    t=ont; [ $s = pb_b3000000 ] && t=pb
    for b in ref lr; do
      rm -rf $W/o.$s.$b
      nice ${BIN[$b]} --db $DB -1 $P/$s/sim/reads/${s}_s_1.fq.gz --read_type $t --no_profile --prefix s -o $W/o.$s.$b -t 1 --no_qcmsa --verbose > $W/o.$s.$b.log 2>&1 || echo "FAIL $s $b"
      echo "$s $b rep $rep: $(grep -E 'Alignment handler took' $W/o.$s.$b.log) | $(grep -E 'Aligning reads took' $W/o.$s.$b.log) | $(grep 'Total alignments' $W/o.$s.$b.log | tr -s '\t ' ' ')"
    done
  done
done
for s in ont_b3000000 pb_b3000000; do echo "== $s"; bash $S/perf3/lrcmp.sh $W/o.$s.ref/s.sam.zst $W/o.$s.lr/s.sam.zst; done
