#!/bin/bash
# Every run of the evaluation, one binary, --mate_guidance guide and consensus each:
#  fin: the simulated test points on the finished reduced r226 database (its own models);
#  tdb: the same on the build's training database (the held-out species absent), with the finished database's pe model
#       (model_logs/trained_model.xml, the database's gate), as the build scored its independent test set;
#  mocks: the six paired-end mock samples on the finished database.
# usage: run_all.sh <protal> <out folder>     (WSL; 4 cores; the rl150 reads made by resim.sh first)
set -u
PROTAL=$1; OUT=$2
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p $OUT
for mode in guide consensus; do
  bash $HERE/run_sims.sh $PROTAL $mode $OUT/fin
  DB=$HOME/c4test/db/work/training_db bash $HERE/run_sims.sh $PROTAL $mode $OUT/tdb --model $HOME/c4test/db/model_logs/trained_model.xml
  bash $HERE/run_mocks.sh $PROTAL $mode $OUT/mocks
done
echo "all done $(date +%T)"
