#!/bin/bash
# protal of one clone on the mock local test's twelve samples in ~/Mocks (the ten drawn and the two extra D6331
# paired-end ones, ../2026-10-13-mock-local-test/samples_extra.tsv) with one reduced database, into the run folder
# ~/Mocks/<run name>, then the scores of every run folder.
# usage: run_mocks.sh <protal clone> <database folder> <run name>     (WSL; 4 cores)
set -eu
SRC=$1; DB=$2; NAME=$3
HERE=$(cd "$(dirname "$0")" && pwd)
PY=$HOME/micromamba/envs/protal-db-build/bin/python
B=$SRC/scripts/mock_benchmark
SAMPLES=$HOME/Mocks/samples_plus.tsv
{ cat $B/samples.tsv; cat $HERE/../2026-10-13-mock-local-test/samples_extra.tsv; } > $SAMPLES
taskset -c 0-3 nice -n 5 $PY $B/score_mocks.py run --reads $HOME/Mocks --samples $SAMPLES --protal $SRC/build/protal \
    --db $DB -t 4 --run-name $NAME
$PY $B/score_mocks.py score --reads $HOME/Mocks --samples $SAMPLES > /dev/null
