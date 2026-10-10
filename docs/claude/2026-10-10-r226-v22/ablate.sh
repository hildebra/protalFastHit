#!/bin/bash
# Retrain v22's models with feature groups left out; usage: ablate.sh <read type: pe|se|pb|ont> <name> <feature set> [seed]
# The trainer is the build's (scripts of 82bcc15, unpacked to ~/v22abl/scripts), with the build's options.
set -u
W=~/v22abl
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal/local/v22/protal0.7.9_r226_v22
TAR=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/ae3c85d7-6e03-4e7f-a312-d9bcaa792d21/scratchpad/scripts_82bcc15.tar
mkdir -p $W/tables/training $W/tables/test $W/out
[ -d $W/scripts ] || tar -xf $TAR -C $W
[ -f $W/internal_taxonomy.dmp ] || cp $SRC/work/internal_taxonomy.dmp $W/
t=$1; name=$2; set=$3; seed=${4:-1}
suf=""; [ $t != pe ] && suf=_$t
for d in training test; do [ -f $W/tables/$d/training_data$suf.tsv ] || cp $SRC/work/$d/training_data$suf.tsv $W/tables/$d/; done
PY=~/micromamba/envs/protal-db-build/bin/python
s=$(date +%s)
taskset -c 0-3 nice -n 5 $PY -I $W/scripts/machine_learning_cmdline.py --truth-file $W/tables/training/training_data$suf.tsv \
  --output-prefix $W/out/${t}_$name --features $set --model gbm --ntree 64 --maxnodes 63 --seed $seed --threads 4 --scenario-weight 0.25 \
  --taxonomy $W/internal_taxonomy.dmp --evaluation basic --depth-knobs --test-file $W/tables/test/training_data$suf.tsv \
  > $W/out/${t}_$name.log 2>&1
echo "$t $name exit $? $(( $(date +%s) - s )) s"
