#!/usr/bin/env bash
# Feature groups (EM, distance, rank) on the congener-rich training tables, both test sets; and depth extrapolation:
# models trained on the training samples of up to 200,000 read pairs, tested on all (the 1M-pair test samples are
# beyond them), the knob curve against the calls at a target share of false calls.
set -euo pipefail
OUT=$HOME/denoise_bench
S=/mnt/c/Users/hildebra/Documents/locDev/protal-denoise/docs/claude/2026-10-03-denoising-implementation
PY=$HOME/micromamba/envs/protal-db-build/bin/python
E=$OUT/eval2
TAX=$OUT/new/internal_taxonomy.dmp
for rt in "" _se; do
  for feats in na+em na+distance na+rank; do
    for test in base new; do
      prefix=$E/new_${feats//+/-}_${test}$rt
      [ -s $prefix.metrics.json ] && continue
      $PY $HOME/denoise_bench/rf2/random_forest_cmdline.py --truth-file $OUT/new/training/training_data$rt.tsv --test-file $OUT/$test/test/training_data$rt.tsv \
        --output-prefix $prefix --features $feats --ntree 64 --maxnodes 512 --seed 1 --threads 6 --taxonomy $TAX \
        --evaluation basic --depth-knobs --fdr-calls > $prefix.log 2>&1 || { echo "failed $prefix"; tail -5 $prefix.log; exit 1; }
      echo "done $(basename $prefix)"
    done
  done
  # depth extrapolation
  shallow=$E/shallow_training$rt.tsv
  [ -s $shallow ] || $PY -c "
import pandas as pd, sys
t = pd.read_csv('$OUT/new/training/training_data$rt.tsv', sep='\t', low_memory=False)
t[t['meta_read_pairs'] <= 200000].to_csv('$shallow', sep='\t', index=False)"
  for feats in normalized+adjacency normalized+adjacency+relatives; do
    prefix=$E/shallow_${feats//+/-}_new$rt
    [ -s $prefix.metrics.json ] && continue
    $PY $HOME/denoise_bench/rf2/random_forest_cmdline.py --truth-file $shallow --test-file $OUT/new/test/training_data$rt.tsv \
      --output-prefix $prefix --features $feats --ntree 64 --maxnodes 512 --seed 1 --threads 6 --taxonomy $TAX \
      --evaluation basic --depth-knobs --fdr-calls > $prefix.log 2>&1 || { echo "failed $prefix"; tail -5 $prefix.log; exit 1; }
    echo "done $(basename $prefix)"
  done
done
