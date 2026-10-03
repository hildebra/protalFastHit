#!/usr/bin/env bash
# Feature-group ablation on the r226 v9 training tables (26b065c binary, sigma 1.3/2.0, 30% of the species held out),
# scored on v9's independent test set, with the trainer of the same commit (~/protal-fp/src = 26b065c + the
# SketchDistance clamp): the default set, without the priors, without the unfiltered group, v8's set
# (nad+depth+divergence) and the default without the depth feature. pe and se. ~1-2 min per model.
set -euo pipefail
PY=$HOME/micromamba/envs/protal-db-build/bin/python
RF=$HOME/protal-fp/src/scripts/random_forest_cmdline.py
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal/local/v9
T=$HOME/v9_tables
OUT=$HOME/v9_eval
mkdir -p $T/training $T/test $OUT
for rt in "" _se; do
  [ -s $T/training/training_data$rt.tsv ] || cp $SRC/training/training_data$rt.tsv $T/training/
  [ -s $T/test/training_data$rt.tsv ] || cp $SRC/test/training_data$rt.tsv $T/test/
done
train() {  # read type suffix, feature set, tag
  local rt=$1 feats=$2 tag=$3
  [ -s $OUT/$tag$rt.log ] && grep -q "Independent test set" $OUT/$tag$rt.log && { echo "have $tag$rt"; return; }
  $PY $RF --truth-file $T/training/training_data$rt.tsv --output-prefix $OUT/$tag$rt --features $feats --ntree 64 \
    --maxnodes 512 --seed 1 --threads 6 --evaluation basic --depth-knobs --test-file $T/test/training_data$rt.tsv \
    > $OUT/$tag$rt.log 2>&1
  echo "done $tag$rt $(date +%T)"
}
for rt in "" _se; do
  train "$rt" normalized+adjacency+distance+depth+divergence+unfiltered+priors default
  train "$rt" normalized+adjacency+distance+depth+divergence+unfiltered no_priors
  train "$rt" normalized+adjacency+distance+depth+divergence+priors no_unfiltered
  train "$rt" normalized+adjacency+distance+depth+divergence v8set
  train "$rt" normalized+adjacency+distance+divergence+unfiltered+priors no_depth
done
echo ABLATIONDONE
