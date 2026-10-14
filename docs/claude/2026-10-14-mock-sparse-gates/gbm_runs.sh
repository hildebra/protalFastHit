#!/bin/bash
# The same mock samples profiled with the database's gradient-boosted models (--model_kind gbm) instead of its default
# random forests, at the database's gates, ungated and on the grid's sparse index at C 4, by score_mocks.py run (a map
# per set), into run folders <label>_<level>-gbm, then scored with the forest runs of collect.py.
# usage: gbm_runs.sh <protal clone> <database folder> <label> <out folder>     (WSL; 4 cores)
set -u
SRC=$1; DB=$2; LABEL=$3; OUT=$4
HERE=$(cd "$(dirname "$0")" && pwd)
PY=$HOME/micromamba/envs/protal-db-build/bin/python
B=$SRC/scripts/mock_benchmark
SAMPLES=$HOME/Mocks/samples_plus.tsv
{ cat $B/samples.tsv; cat $HERE/../2026-10-13-mock-local-test/samples_extra.tsv; } > $SAMPLES
for spec in "dbgates:" "C0:--sparse_first 0" "C4g:--sparse_first 4 --sparse_first_pick grid"; do
  level=${spec%%:*}; args=${spec#*:}
  taskset -c 0-3 nice -n 5 $PY $B/score_mocks.py run --reads $HOME/Mocks --samples $SAMPLES --protal $SRC/build/protal \
      --db $DB -t 4 --run-name ${LABEL}_$level-gbm --protal-args "--model_kind gbm $args" > $OUT/gbm_$level.log 2>&1
  echo "$level: exit $?"
done
runs=$(printf '%s,' ${LABEL}_dbgates ${LABEL}_C0 ${LABEL}_C4s ${LABEL}_C4g ${LABEL}_dbgates-gbm ${LABEL}_C0-gbm ${LABEL}_C4g-gbm)
$PY $B/score_mocks.py score --reads $HOME/Mocks --samples $SAMPLES --runs ${runs%,} --no-sylph > /dev/null
rm -f $SAMPLES
cp $HOME/Mocks/scores/samples_by_run.tsv $OUT/samples_by_run_all.tsv
cp $HOME/Mocks/scores/summary_by_run.tsv $OUT/summary_by_run_all.tsv
