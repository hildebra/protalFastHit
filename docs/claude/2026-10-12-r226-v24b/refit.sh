#!/bin/bash
# Refit v24b's pe model on its own tables, whole or with one class of present strains thinned (thin_table.py), with the
# build's trainer (scripts at 4a42f3b, from git archive) and options; species-held-out folds only (ABL_SCHEMES, as in
# the v23 refits: the species folds are GroupKFold, the same either way). WSL, cores 0-3.
# usage: refit.sh <name> [real|insilico <share kept>]
set -u
W=~/v24babl
R=/mnt/c/Users/hildebra/Documents/locDev/protal
SRC=$R/local/v24b/protal0.8.1_r226_v24
mkdir -p $W/tables/training $W/tables/test $W/out
if [ ! -d $W/scripts ]; then
  git -C $R archive 4a42f3b scripts | tar x -C $W
  perl -pi -e 's/for scheme in \("rows", "samples", "species", \*\(CLADE_SCHEMES if opts.evaluation == "full" else \(\)\)\):/for scheme in [s for s in ("rows", "samples", "species") if s in os.environ.get("ABL_SCHEMES", "species")]:/' \
    $W/scripts/machine_learning_cmdline.py
  grep -q ABL_SCHEMES $W/scripts/machine_learning_cmdline.py || { echo "patch failed"; exit 1; }
fi
[ -f $W/internal_taxonomy.dmp ] || cp $SRC/work/internal_taxonomy.dmp $W/
for d in training test; do [ -f $W/tables/$d/training_data.tsv ] || cp $SRC/work/$d/training_data.tsv $W/tables/$d/; done
name=$1
table=$W/tables/training/training_data.tsv
if [ $# -ge 3 ]; then
  table=$W/tables/training/training_data_$name.tsv
  [ -f $table ] || python3 $R/docs/claude/2026-10-12-r226-v24b/thin_table.py $W/tables/training/training_data.tsv $table $2 $3 1
fi
[ -f $W/out/pe_$name.calls.tsv.gz ] && { echo "pe $name done before"; exit 0; }
PY=~/micromamba/envs/protal-db-build/bin/python
s=$(date +%s)
taskset -c 0-3 nice -n 5 $PY -I $W/scripts/machine_learning_cmdline.py --truth-file $table \
  --output-prefix $W/out/pe_$name \
  --features normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights \
  --model gbm --ntree 64 --maxnodes 63 --seed 1 --threads 4 --scenario-weight 0.25 --taxonomy $W/internal_taxonomy.dmp \
  --evaluation basic --depth-knobs --test-file $W/tables/test/training_data.tsv > $W/out/pe_$name.log 2>&1
echo "pe $name exit $? $(( $(date +%s) - s )) s"
