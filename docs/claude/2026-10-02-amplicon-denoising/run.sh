#!/usr/bin/env bash
# The experiment of this report, in WSL with the protal-db-build env (python 3.14, scikit-learn 1.9.1).
# base: normalized+adjacency; rel: + genus_skew, family_skew, genus_share, genus_spill; spill: + spill_divergent.
# The trainer is scripts/ of 6d2e241 with model_features.patch (two more feature sets). ~1 min per model.
set -euo pipefail
PY=${PY:-/home/falk/micromamba/envs/protal-db-build/bin/python}
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
R=$REPO/docs/claude/2026-10-02-amplicon-denoising
D=$REPO/local/protal0.7.3_r226_v3          # the r226 v3 build's training/ and test/ tables (git-ignored)
W=${W:-$HOME/denoise_r226}
TAX=$D/internal_taxonomy.dmp
mkdir -p $W/tables $W/out $W/rf
for f in random_forest_cmdline.py lineages.py model_features.py model_pmml.py; do
  git -C $REPO show 6d2e241:scripts/$f > $W/rf/$f
done
patch -d $W/rf -p1 < $R/model_features.patch

for rt in "" _se _pb _ont; do
  for set in training test; do
    $PY $R/add_relative_features.py $TAX $D/$set/training_data$rt.tsv $W/tables/${set}$rt.tsv "" $W/rf
    $PY $R/add_spill.py $W/tables/${set}$rt.tsv $W/tables/${set}$rt.spill.tsv
  done
done

train() {  # read type, feature set, tag, table suffix
  local rt=$1 feats=$2 tag=$3 sfx=$4 leaves=512
  [[ $rt == _pb || $rt == _ont ]] && leaves=128
  $PY $W/rf/random_forest_cmdline.py --truth-file $W/tables/training$rt$sfx.tsv --output-prefix $W/out/$tag$rt \
    --features $feats --ntree 64 --maxnodes $leaves --seed 1 --threads 6 --taxonomy $TAX --evaluation basic \
    --depth-knobs --test-file $W/tables/test$rt$sfx.tsv > $W/out/$tag$rt.log 2>&1
  echo "done $tag$rt"
}
for rt in "" _se _pb _ont; do
  train "$rt" normalized+adjacency base ""
  train "$rt" normalized+adjacency+relatives rel ""
  train "$rt" normalized+adjacency+spill spill .spill
done

$PY $R/describe.py $W/tables/training.tsv $W/out/base.predictions.tsv.gz p_species > $R/describe_species_held_out.txt
$PY $R/describe.py $W/tables/test.tsv $W/out/base.test_predictions.tsv.gz p > $R/describe_test.txt
$PY $R/ambiguity.py $W/tables > $R/ambiguity_output.txt
for rt in "" _se; do $PY $R/congeners.py $W/out $W/tables "$rt" base,rel,spill; done > $R/congeners_output.txt
$PY $R/compare.py $W/out $W/tables base,rel,spill > $R/compare_output.txt
