#!/bin/bash
# fit v19's feature set on v19's training tables (fit only), then score v19's and v21's test tables with it and with
# v21's no_new model; waits for the ablation queue first
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/b5dc7ad1-8bd4-4d3a-a544-16771ad084bb/scratchpad
until grep -q "MORE DONE" $1; do sleep 20; done
W=~/v21abl; V19=/mnt/c/Users/hildebra/Documents/locDev/protal/local/v19/protal0.7.9_r226_v19
PY=~/micromamba/envs/protal-db-build/bin/python
B=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried
mkdir -p $W/v19tables $W/out19
for t in pb ont pe se; do
  suf=""; [ $t != pe ] && suf=_$t
  for d in training test; do mkdir -p $W/v19tables/$d; [ -f $W/v19tables/$d/training_data$suf.tsv ] || cp $V19/work/$d/training_data$suf.tsv $W/v19tables/$d/; done
  taskset -c 0-3 nice -n 5 $PY -I $W/scripts/machine_learning_cmdline.py --truth-file $W/v19tables/training/training_data$suf.tsv \
    --output-prefix $W/out19/${t}_v19 --features $B --model gbm --maxnodes 63 --seed 1 --threads 4 --scenario-weight 0.25 \
    --taxonomy $V19/work/internal_taxonomy.dmp --evaluation none > $W/out19/${t}_v19.log 2>&1
  echo "$t v19 fit exit $?"
  taskset -c 0-3 nice -n 5 $PY -I $S/cross_score.py $t v19_model=$W/out19/${t}_v19 v21_model=$W/out/${t}_no_new -- \
    v19_test=$W/v19tables/test/training_data$suf.tsv v21_test=$W/tables/test/training_data$suf.tsv 2>&1 | grep -v Warning
done
echo CROSS DONE
