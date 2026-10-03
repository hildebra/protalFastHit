#!/usr/bin/env bash
# Every model on both test sets: the training tables of the two pipelines (base: species drawn uniformly; new: congener
# groups 0.25:2-5), each with the features before the relatives features and with them, each tested on both test sets.
set -euo pipefail
OUT=$HOME/denoise_bench
S=$HOME/protal-denoise/src/scripts
PY=$HOME/micromamba/envs/protal-db-build/bin/python
E=$OUT/eval2
mkdir -p $E
table() {  # pipeline kind(training|test) read-type-suffix
  local f
  f=$(find $OUT/$1 -path "*/$2/training_data$3.tsv" | head -1)
  [ -n "$f" ] || { echo "no $2 table of $1" >&2; exit 1; }
  echo $f
}
for rt in "" _se; do
  for train in base new; do
    TAX=$OUT/$train/training_db/internal_taxonomy.dmp
    [ -f $TAX ] || TAX=$(find $OUT/$train -name internal_taxonomy.dmp | head -1)
    for feats in normalized+adjacency normalized+adjacency+relatives; do
      for test in base new; do
        prefix=$E/${train}_${feats//+/-}_${test}$rt
        [ -s $prefix.metrics.json ] && continue
        $PY $HOME/denoise_bench/rf2/random_forest_cmdline.py --truth-file $(table $train training "$rt") --test-file $(table $test test "$rt") \
          --output-prefix $prefix --features $feats --ntree 64 --maxnodes 512 --seed 1 --threads 6 --taxonomy $TAX \
          --evaluation basic --depth-knobs --fdr-calls > $prefix.log 2>&1 || { echo "failed: $prefix"; tail -20 $prefix.log; exit 1; }
        echo "done $(basename $prefix)"
      done
    done
  done
done
