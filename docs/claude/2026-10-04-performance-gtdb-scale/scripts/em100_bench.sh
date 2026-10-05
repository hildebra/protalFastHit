#!/usr/bin/env bash
# The read EM capped at 100 sweeps (after) against 200 (before, 836e4a3), on the 0.7.5 benchmark (tie_bench.sh's runs, in ~/embench).
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
F=$HOME/embench; mkdir -p $F/bin
[ -x $F/bin/protal-before ] || cp $HOME/perf-gtdb/work/build/protal $F/bin/protal-before
bash $S/build_work.sh "*" 8 || exit 1
cp $HOME/perf-gtdb/work/build/protal $F/bin/protal-after
TIES=$F bash $S/ties/tie_bench.sh
cd $F && $HOME/micromamba/envs/protal-db-build/bin/python /mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-03-v075-benchmark/scripts/score.py $F short lr > $F/score.log 2>&1 || { echo "scoring failed"; tail -20 $F/score.log; }
for s in results_v075 results_lr075; do echo "== $s"; cat $F/$s/overall.md; cat $F/$s/paired.md; done
grep -h "read EM:" $F/runs_v075/v07*.full.pe.rl150_p500000_s_1.log | cut -c1-200
