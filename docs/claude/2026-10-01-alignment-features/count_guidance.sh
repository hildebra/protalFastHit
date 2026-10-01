#!/bin/bash
# Instructions of the paired-end alignment loop (the OpenMP body of RunPairedEnd, callgrind) of the first PAIRS pairs of a test sample, with
# and without mate guidance (the same binary, --no_mate_guidance): load-independent, unlike wall-clock times on a
# shared machine. One thread, pinned to CPU 0, lowest priority.
#     bash count_guidance.sh RUN SAMPLE_R1.fq.gz OUT [PAIRS]
set -eu
RUN=$1
R1=$2
R2=${R1/_R1/_R2}
OUT=$3
PAIRS=${4:-40000}
BIN=~/protal-cons/build/protal
mkdir -p $OUT
zcat $R1 | head -n $((4 * PAIRS)) > $OUT/R1.fq
zcat $R2 | head -n $((4 * PAIRS)) > $OUT/R2.fq
for variant in guided unguided; do
  extra=""; [ $variant = unguided ] && extra="--no_mate_guidance"
  nice -n 19 taskset -c 0 valgrind --tool=callgrind --toggle-collect='*RunPairedEnd*_omp_fn*' \
      --callgrind-out-file=$OUT/cg.$variant $BIN --db $RUN/training_db -1 $OUT/R1.fq -2 $OUT/R2.fq \
      -o $OUT/out.$variant -t 1 --no_profile $extra > $OUT/log.$variant 2> $OUT/err.$variant
  total=$(callgrind_annotate $OUT/cg.$variant 2>/dev/null | grep -m1 -E "PROGRAM TOTALS|Ir" | head -1)
  guide=$(callgrind_annotate --inclusive=yes $OUT/cg.$variant 2>/dev/null | grep -m1 "GuideMate" || true)
  echo "$variant: $total"
  echo "   GuideMate (inclusive): ${guide:-none}"
  grep -h "fragments had one mate" $OUT/log.$variant || true
done
