#!/bin/bash
# The local test: ten samples drawn at random (pick_samples.py, seed 20261013: one sample of each of ten sets with a
# truth table), their first 1M read pairs (long reads: 1 Gb) fetched into ~/Mocks, profiled by protal on the reduced
# database and scored.
# usage: run_test.sh <protal clone> <database folder> [download|run|score|all]     (WSL; 4 cores)
set -eu
SRC=$1; DB=$2; STEP=${3:-all}
MOCKS=$HOME/Mocks
PY=$HOME/micromamba/envs/protal-db-build/bin/python
B=$SRC/scripts/mock_benchmark
HERE=$(cd "$(dirname "$0")" && pwd)
SAMPLES=$($PY $HERE/pick_samples.py $B/samples.tsv | cut -f2 | paste -sd'|')
if [ $STEP = download ] || [ $STEP = all ]; then
  taskset -c 0-3 $PY $B/download_mocks.py -o $MOCKS --only "^($SAMPLES)\$" --max-pairs 1e6 --max-bases 1e9 --jobs 4
fi
if [ $STEP = run ] || [ $STEP = all ]; then
  taskset -c 0-3 nice -n 5 $PY $B/score_mocks.py run --reads $MOCKS --protal $SRC/build/protal --db $DB -t 4
fi
if [ $STEP = score ] || [ $STEP = all ]; then
  $PY $B/score_mocks.py score --reads $MOCKS
fi
