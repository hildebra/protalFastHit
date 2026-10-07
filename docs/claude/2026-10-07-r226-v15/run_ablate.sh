#!/bin/bash
# The complexity ablation on v15's tables, 4 pinned cores (the machine is shared).
R=/mnt/c/Users/hildebra/Documents/locDev/protal
D=$R/docs/claude/2026-10-07-r226-v15
cd $R
taskset -c 0-3 nice -n 5 ~/soil13/venv/bin/python $D/ablate.py --build local/v15 --out ~/v15/ablate \
    --read-types pe,ont,se,pb -t 4
echo "ablate exit $?"
