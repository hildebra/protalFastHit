#!/bin/bash
# scikit-learn vs LightGBM vs XGBoost on the r226 v14 paired-end table, pinned to 4 cores, after the trainer profile.
set -u
W=~/mlcpp
while ! grep -q '^done' $W/out/progress.txt; do sleep 10; done
B=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-07-cpp-training/bench_libraries.py
cp $B $W/bench_libraries.py
export PROTAL_SCRIPTS=$W/src/scripts
P=$W/venv/bin/python
run() { echo "start $1 $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" >> $W/out/progress.txt; shift;
        taskset -c 0-3 nice -n 5 "$@"; }
run fit_times $P $W/bench_libraries.py --truth-file $W/data/train_pe.tsv --threads 1,2,4 --repeats 2 --out $W/out/fit_times.tsv > $W/out/fit_times.log 2>&1
run quality $P $W/bench_libraries.py --truth-file $W/data/train_pe.tsv --threads 4 --folds 5 --out $W/out/quality.tsv > $W/out/quality.log 2>&1
run small $P $W/bench_libraries.py --truth-file $W/data/train_pe.tsv --rows 4000 --threads 1,4 --repeats 2 --out $W/out/small.tsv > $W/out/small.log 2>&1
run forest $P $W/bench_libraries.py --truth-file $W/data/train_pe.tsv --threads 4 --folds 5 --forest --out $W/out/forest.tsv > $W/out/forest.log 2>&1
echo "bench done $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" >> $W/out/progress.txt
