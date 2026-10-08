#!/bin/bash
# Species folds re-timed (LightGBM's first run fell on another session's tests), twice, the libraries' order alternated.
set -u
W=~/mlcpp
export PROTAL_SCRIPTS=$W/src/scripts
P=$W/venv/bin/python
cp /mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-07-cpp-training/bench_libraries.py $W/bench_libraries.py
for i in 1 2; do
  case $i in 1) L=lightgbm,sklearn,xgboost ;; 2) L=xgboost,sklearn,lightgbm ;; esac
  echo "start folds$i $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" >> $W/out/progress.txt
  taskset -c 0-3 nice -n 5 $P $W/bench_libraries.py --truth-file $W/data/train_pe.tsv --libraries $L --threads 4 --repeats 0 \
    --folds 5 --out $W/out/folds$i.tsv > $W/out/folds$i.log 2>&1
  ps -eo pcpu,args --sort=-pcpu | head -4 | cut -c1-120 >> $W/out/progress.txt
done
echo "folds done $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg)" >> $W/out/progress.txt
