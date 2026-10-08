#!/bin/bash
# Profile the trainer (current defaults) on the r226 v14 paired-end table, pinned to 4 cores.
set -u
W=~/mlcpp
mkdir -p $W/data $W/out
rsync -a --checksum --no-times /mnt/c/Users/hildebra/Documents/locDev/protal/scripts/ $W/src/scripts/
V=/mnt/c/Users/hildebra/Documents/locDev/protal/local/v14
[ -f $W/data/train_pe.tsv ] || cp $V/training/training_data.tsv $W/data/train_pe.tsv
[ -f $W/data/test_pe.tsv ] || cp $V/test/training_data.tsv $W/data/test_pe.tsv
cp $V/internal_taxonomy.dmp $W/data/
P=~/micromamba/envs/protal-db-build/bin/python
cd $W/src
COMMON="--truth-file $W/data/train_pe.tsv --test-file $W/data/test_pe.tsv --taxonomy $W/data/internal_taxonomy.dmp --features normalized+adjacency+distance+depth+divergence+unfiltered+ref --model gbm --seed 1 --threads 4 --scenario-weight 0.25 --evaluation basic --depth-knobs"
echo "start seq $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" > $W/out/progress.txt
taskset -c 0-3 nice -n 5 /usr/bin/time -v $P -m cProfile -o $W/out/seq.prof scripts/machine_learning_cmdline.py $COMMON --fold-jobs 1 --output-prefix $W/out/seq > $W/out/seq.log 2>&1
echo "EXIT $?" >> $W/out/seq.log
echo "start par $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" >> $W/out/progress.txt
taskset -c 0-3 nice -n 5 /usr/bin/time -v $P scripts/machine_learning_cmdline.py $COMMON --output-prefix $W/out/par > $W/out/par.log 2>&1
echo "EXIT $?" >> $W/out/par.log
echo "done $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" >> $W/out/progress.txt
