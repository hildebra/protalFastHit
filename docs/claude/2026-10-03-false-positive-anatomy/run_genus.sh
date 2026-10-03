#!/usr/bin/env bash
# Genus size as a per-species feature, with and without the sample depth; runs after run_depth.sh (waits for ALLDONE).
set -euo pipefail
PY=/home/falk/micromamba/envs/protal-db-build/bin/python
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
D=$REPO/local/v5
W=$HOME/fpdepth
TAX=$D/internal_taxonomy.dmp
while ! grep -q ALLDONE $W/run.log; do sleep 30; done
cd $W/rf && patch -p1 --forward < $W/model_features_genus.patch || true
for rt in "" _se; do
  for set in training test; do
    [[ -s $W/tables/${set}$rt.g.tsv ]] || $PY $W/add_genus_size.py $TAX $W/tables/${set}$rt.tsv $W/tables/${set}$rt.g.tsv $W/rf
  done
done
train() {  # read type, feature set, tag
  local rt=$1 feats=$2 tag=$3 leaves=512
  $PY $W/rf/random_forest_cmdline.py --truth-file $W/tables/training$rt.g.tsv --output-prefix $W/out/$tag$rt \
    --features $feats --ntree 64 --maxnodes $leaves --seed 1 --threads 6 --taxonomy $TAX --evaluation basic \
    --depth-knobs --test-file $W/tables/test$rt.g.tsv > $W/out/$tag$rt.log 2>&1
  echo "done $tag$rt $(date +%T)"
}
for rt in "" _se; do
  train "$rt" normalized+adjacency+genus na_genus
  train "$rt" normalized+adjacency+depth+genus na_depth_genus
  train "$rt" normalized+adjacency+distance+depth+genus nad_depth_genus
done
echo GENUSDONE
