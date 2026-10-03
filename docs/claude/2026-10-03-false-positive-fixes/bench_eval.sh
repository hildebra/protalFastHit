#!/usr/bin/env bash
# Feature-set ablation on the benchmark-world pipeline's tables (~/fp_bench/new, built by the fp-anatomy protal): the
# default set against the sets without the depth and divergence groups, each with --depth-knobs (fitted only where
# the depth is not a feature), scored on the pipeline's test set. ~1 min per model.
set -euo pipefail
PY=$HOME/micromamba/envs/protal-db-build/bin/python
RF=$HOME/protal-fp/src/scripts/random_forest_cmdline.py
P=$HOME/fp_bench/new
OUT=$HOME/fp_bench/eval
mkdir -p $OUT
train() {  # read type suffix, feature set, tag
  local rt=$1 feats=$2 tag=$3
  [ -s $OUT/$tag$rt.log ] && grep -q "Independent test set" $OUT/$tag$rt.log && { echo "have $tag$rt"; return; }
  $PY $RF --truth-file $P/training/training_data$rt.tsv --output-prefix $OUT/$tag$rt --features $feats --ntree 64 \
    --maxnodes 512 --seed 1 --threads 6 --evaluation basic --depth-knobs --test-file $P/test/training_data$rt.tsv \
    > $OUT/$tag$rt.log 2>&1
  echo "done $tag$rt $(date +%T)"
}
for rt in "" _se; do
  train "$rt" normalized+adjacency+distance nad
  train "$rt" normalized+adjacency+distance+depth nad_depth
  train "$rt" normalized+adjacency+distance+divergence nad_div
  train "$rt" normalized+adjacency+distance+depth+divergence nad_depth_div
done
echo EVALDONE
