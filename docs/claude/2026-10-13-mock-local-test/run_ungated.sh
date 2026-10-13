#!/bin/bash
set -eu
B=/mnt/c/Users/hildebra/Documents/locDev/protal/scripts/mock_benchmark
PY=$HOME/micromamba/envs/protal-db-build/bin/python
taskset -c 0-3 nice -n 5 $PY $B/score_mocks.py run --reads ~/Mocks --protal ~/mocktest/protal/build/protal \
    --db ~/mocktest/build/protal_db -t 4 --run-name 0.8.3_4dfc549_ungated --protal-args "--sparse_first 0" 2>&1 | grep -v "^  /"
$PY $B/score_mocks.py score --reads ~/Mocks > ~/Mocks/scores/by_run_print.tsv
