#!/usr/bin/env bash
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
F=$HOME/tiebench; for d in world samples samples_deep samples_lr073; do ln -sfn $HOME/bench071/$d $F/$d; done
bash $S/ties/alignment_accuracy.sh > $F/alignment_accuracy.txt; cat $F/alignment_accuracy.txt
cd $F && $HOME/micromamba/envs/protal-db-build/bin/python /mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-03-v075-benchmark/scripts/score.py $F short lr > $F/score.log 2>&1 || { echo "scoring failed"; tail -20 $F/score.log; }
for s in results_v075 results_lr075; do echo "== $s"; cat $F/$s/overall.md; cat $F/$s/paired.md; done
