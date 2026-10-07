#!/bin/bash
# The ablations and the test rows scored by species fold models, in one queue on 4 pinned cores (the machine is
# shared and was busy): per read type its missing variants, then its test rows by fold models; pe, ont, pb, se.
R=/mnt/c/Users/hildebra/Documents/locDev/protal
D=$R/docs/claude/2026-10-07-r226-v15
PY="taskset -c 0-3 nice -n 5 $HOME/soil13/venv/bin/python"
cd $R
mkdir -p ~/v15/tsho
for step in "pe no_ref" "ont no_ref" "pb v15,no_complexity,no_ref" "se no_complexity,no_ref"; do
    set -- $step
    $PY $D/ablate.py --build local/v15 --out ~/v15/ablate --read-types $1 --variants $2 -t 4
    $PY $D/test_species_held_out.py --build local/v15 --ablate ~/v15/ablate --variants v15,no_ref,no_complexity \
        --read-types $1 --out ~/v15/tsho -t 4
    echo "$1 done $?"
done
echo "queue done"
