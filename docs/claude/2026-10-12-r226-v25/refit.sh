#!/bin/bash
# Refit v24b's pe or se model (the shift on, the hard fragment count) on its own tables, as they are or with the
# fragments feature replaced by the soft count's approximation (approx_soft.py write), with v24b's trainer (scripts at
# 4a42f3b, from git archive) and options; species-held-out folds only (ABL_SCHEMES, as in the v24b refits). WSL, cores 0-3.
# usage: refit.sh <pe|se> <hard|soft> [seed]
set -u
W=~/v25abl
R=/mnt/c/Users/hildebra/Documents/locDev/protal
SRC=$R/local/v24b/protal0.8.1_r226_v24
rt=$1; variant=$2; seed=${3:-1}
sfx=""; [ "$rt" = pe ] || sfx="_$rt"
mkdir -p $W/tables/training $W/tables/test $W/out
if [ ! -d $W/scripts ]; then
  git -C $R archive 4a42f3b scripts | tar x -C $W
  perl -pi -e 's/for scheme in \("rows", "samples", "species", \*\(CLADE_SCHEMES if opts.evaluation == "full" else \(\)\)\):/for scheme in [s for s in ("rows", "samples", "species") if s in os.environ.get("ABL_SCHEMES", "species")]:/' \
    $W/scripts/machine_learning_cmdline.py
  grep -q ABL_SCHEMES $W/scripts/machine_learning_cmdline.py || { echo "patch failed"; exit 1; }
fi
[ -f $W/internal_taxonomy.dmp ] || cp $SRC/work/internal_taxonomy.dmp $W/
for d in training test; do
  hard=$W/tables/$d/training_data${sfx}_hard.tsv
  [ -f $hard ] || cp $SRC/work/$d/training_data$sfx.tsv $hard
  soft=$W/tables/$d/training_data${sfx}_soft.tsv
  if [ "$variant" = soft ] && [ ! -f $soft ]; then
    taskset -c 0-3 nice -n 5 python3 $R/docs/claude/2026-10-12-r226-v25/approx_soft.py write $hard $soft || exit 1
  fi
done
name=${rt}_${variant}_s$seed
[ -f $W/out/$name.calls.tsv.gz ] && { echo "$name done before"; exit 0; }
PY=~/micromamba/envs/protal-db-build/bin/python
s=$(date +%s)
taskset -c 0-3 nice -n 5 $PY -I $W/scripts/machine_learning_cmdline.py --truth-file $W/tables/training/training_data${sfx}_$variant.tsv \
  --output-prefix $W/out/$name \
  --features normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights \
  --model gbm --ntree 64 --maxnodes 63 --seed $seed --threads 4 --scenario-weight 0.25 --taxonomy $W/internal_taxonomy.dmp \
  --evaluation basic --depth-knobs --test-file $W/tables/test/training_data${sfx}_$variant.tsv > $W/out/$name.log 2>&1
echo "$name exit $? $(( $(date +%s) - s )) s"
