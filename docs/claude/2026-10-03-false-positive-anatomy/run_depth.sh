#!/usr/bin/env bash
# Sample depth as a feature: trainer of 27423c6 + model_features_depth.patch, on the r226 v5 tables (local/v5).
set -euo pipefail
PY=/home/falk/micromamba/envs/protal-db-build/bin/python
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
S="/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/3611c5a5-4e2b-44f3-a684-74ae56dc8a7c/scratchpad"
D=$REPO/local/v5
W=$HOME/fpdepth
TAX=$D/internal_taxonomy.dmp
cd $W/rf && patch -p1 --forward < $S/model_features_depth.patch || true
for rt in "" _se _pb _ont; do
  for set in training test; do
    [[ -s $W/tables/${set}$rt.tsv ]] || $PY $S/add_sample_depth.py $D/$set/training_data$rt.tsv $W/tables/${set}$rt.tsv
  done
done
train() {  # read type, feature set, tag
  local rt=$1 feats=$2 tag=$3 leaves=512
  [[ $rt == _pb || $rt == _ont ]] && leaves=128
  $PY $W/rf/random_forest_cmdline.py --truth-file $W/tables/training$rt.tsv --output-prefix $W/out/$tag$rt \
    --features $feats --ntree 64 --maxnodes $leaves --seed 1 --threads 6 --taxonomy $TAX --evaluation basic \
    --depth-knobs --test-file $W/tables/test$rt.tsv > $W/out/$tag$rt.log 2>&1
  echo "done $tag$rt $(date +%T)"
}
for rt in "" _se; do
  train "$rt" normalized+adjacency na
  train "$rt" normalized+adjacency+depth na_depth
  train "$rt" normalized+adjacency+distance nad
  train "$rt" normalized+adjacency+distance+depth nad_depth
done
for rt in _pb _ont; do
  train "$rt" normalized+adjacency+distance nad
  train "$rt" normalized+adjacency+distance+depth nad_depth
done
echo ALLDONE
