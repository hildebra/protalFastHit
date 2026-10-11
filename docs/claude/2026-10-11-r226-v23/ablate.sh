#!/bin/bash
# Retrain v23's models with feature groups left out; usage: ablate.sh <read type: pe|se|pb|ont> <name> <feature set> [seed]
# The trainer is the build's (machine_learning_cmdline.py and model_features.py are the same at 82bcc15 and 8d8c89b,
# so ~/v22abl/scripts serves), with the build's options; 3 threads on cores 0-2 (the fourth core runs the profiler tests;
# CORES and THREADS change that). The copy in ~/v23abl/scripts was patched to fit only the species-held-out folds (the
# rows and samples cross-validations skipped: the species folds are GroupKFold, the same either way):
#   for scheme in [s for s in ("rows", "samples", "species", ...) if s in os.environ.get("ABL_SCHEMES", "species")]:
set -u
W=~/v23abl
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal/local/v23/protal0.7.9_r226_v23
mkdir -p $W/tables/training $W/tables/test $W/out
[ -d $W/scripts ] || cp -r ~/v22abl/scripts $W/scripts
[ -f $W/internal_taxonomy.dmp ] || cp $SRC/work/internal_taxonomy.dmp $W/
t=$1; name=$2; set=$3; seed=${4:-1}
suf=""; [ $t != pe ] && suf=_$t
for d in training test; do [ -f $W/tables/$d/training_data$suf.tsv ] || cp $SRC/work/$d/training_data$suf.tsv $W/tables/$d/; done
[ -f $W/out/${t}_$name.calls.tsv.gz ] && { echo "$t $name done before"; exit 0; }
PY=~/micromamba/envs/protal-db-build/bin/python
s=$(date +%s)
taskset -c ${CORES:-0-2} nice -n 5 $PY -I $W/scripts/machine_learning_cmdline.py --truth-file $W/tables/training/training_data$suf.tsv \
  --output-prefix $W/out/${t}_$name --features $set --model gbm --ntree 64 --maxnodes 63 --seed $seed --threads ${THREADS:-3} --scenario-weight 0.25 \
  --taxonomy $W/internal_taxonomy.dmp --evaluation basic --depth-knobs --test-file $W/tables/test/training_data$suf.tsv \
  > $W/out/${t}_$name.log 2>&1
echo "$t $name exit $? $(( $(date +%s) - s )) s"
