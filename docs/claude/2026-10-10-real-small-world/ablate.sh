#!/bin/bash
# Refit the small world's model of one read type with groups of features left out (the build's trainer and options);
# usage: ablate.sh <build outdir> <pe|pb> <name> <feature set>. Outputs in <build outdir>_abl/.
set -u
B=$1; t=$2; name=$3; set=$4
W=${B}_abl
mkdir -p $W
suf=""; [ $t != pe ] && suf=_$t
PY=~/micromamba/envs/protal-db-build/bin/python
s=$(date +%s)
taskset -c 0-3 nice -n 5 $PY -I ~/bidx/src/scripts/machine_learning_cmdline.py --truth-file $B/work/training/training_data$suf.tsv \
  --output-prefix $W/${t}_$name --features $set --model gbm --ntree 64 --maxnodes 63 --seed 1 --threads 4 --scenario-weight 0.25 \
  --taxonomy $B/work/internal_taxonomy.dmp --evaluation basic --depth-knobs --test-file $B/work/test/training_data$suf.tsv \
  > $W/${t}_$name.log 2>&1
echo "$t $name exit $? $(( $(date +%s) - s )) s"
