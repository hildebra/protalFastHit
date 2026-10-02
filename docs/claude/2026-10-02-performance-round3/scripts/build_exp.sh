#!/usr/bin/env bash
# build_exp.sh NAME: the perf3 worktree's changes (git diff HEAD, written to $S/perf3/NAME.patch by the caller)
# on HEAD's archive, built in ~/mt-work/perf3-NAME (protal).
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/build_wt.sh perf3-$1 $S/perf3/$1.patch $(cat $HOME/mt-work/perf3/ref/COMMIT) "${2:-protal}"
