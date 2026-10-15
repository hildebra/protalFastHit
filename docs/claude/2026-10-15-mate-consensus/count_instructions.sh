#!/bin/bash
# Instructions of the paired-end alignment loop (the OpenMP body of RunPairedEnd, callgrind) on the first PAIRS pairs
# of a sample, with --mate_guidance off, guide and consensus (one binary): load-independent, unlike wall-clock times
# on the shared machine. One thread, lowest priority; no profile.
# usage: count_instructions.sh <protal> <database folder> <R1> <R2> <out folder> [PAIRS]     (WSL)
set -eu
BIN=$1; DB=$2; R1=$3; R2=$4; OUT=$5; PAIRS=${6:-100000}
mkdir -p $OUT
cat_any() { case $1 in *.zst) zstdcat $1 ;; *.gz) zcat $1 ;; *) cat $1 ;; esac; }
cat_any $R1 | head -n $((4 * PAIRS)) > $OUT/R1.fq || true
cat_any $R2 | head -n $((4 * PAIRS)) > $OUT/R2.fq || true
for mode in off guide consensus; do
  nice -n 19 taskset -c 0-3 valgrind --tool=callgrind --toggle-collect='*RunPairedEnd*_omp_fn*' \
      --callgrind-out-file=$OUT/cg.$mode $BIN --db $DB -1 $OUT/R1.fq -2 $OUT/R2.fq -o $OUT/out.$mode -t 1 --no_profile \
      --force --mate_guidance $mode > $OUT/log.$mode 2> $OUT/err.$mode
  total=$(callgrind_annotate $OUT/cg.$mode 2>/dev/null | grep -m1 "PROGRAM TOTALS" || true)
  echo "$mode: $total"
  callgrind_annotate --inclusive=yes $OUT/cg.$mode 2>/dev/null | grep -E "VoteMates|GuideMate|RescueMate|BestDiagonal" | head -6 || true
  grep -h "fragments had one mate\|Mate consensus" $OUT/log.$mode || true
done
rm -f $OUT/R1.fq $OUT/R2.fq
