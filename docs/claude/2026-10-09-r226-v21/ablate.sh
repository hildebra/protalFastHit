#!/bin/bash
# Retrain v21's models with feature groups left out; usage: ablate.sh <read type suffix: pe|se|pb|ont> <name> <feature set> [seed]
set -u
W=~/v21abl
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal/local/v21/protal0.7.9_r226_v21
SNAP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/b5dc7ad1-8bd4-4d3a-a544-16771ad084bb/scratchpad/v21snap/scripts
mkdir -p $W/tables/training $W/tables/test $W/out
[ -d $W/scripts ] || cp -r $SNAP $W/scripts
[ -f $W/internal_taxonomy.dmp ] || cp $SRC/work/internal_taxonomy.dmp $W/
t=$1; name=$2; set=$3; seed=${4:-1}
suf=""; [ $t != pe ] && suf=_$t
for d in training test; do [ -f $W/tables/$d/training_data$suf.tsv ] || cp $SRC/work/$d/training_data$suf.tsv $W/tables/$d/; done
PY=~/micromamba/envs/protal-db-build/bin/python
s=$(date +%s)
taskset -c 0-3 nice -n 5 $PY -I $W/scripts/machine_learning_cmdline.py --truth-file $W/tables/training/training_data$suf.tsv \
  --output-prefix $W/out/${t}_$name --features $set --model gbm --maxnodes 63 --seed $seed --threads 4 --scenario-weight 0.25 \
  --taxonomy $W/internal_taxonomy.dmp --evaluation basic --depth-knobs --test-file $W/tables/test/training_data$suf.tsv \
  > $W/out/${t}_$name.log 2>&1
echo "$t $name exit $? $(( $(date +%s) - s )) s"
