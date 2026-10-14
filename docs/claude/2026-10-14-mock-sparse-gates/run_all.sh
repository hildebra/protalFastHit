#!/bin/bash
# The local mock test again with one protal commit, its reads gated three ways: protal and simulate_metagenomes of the
# commit (a clone), the reduced r226 database's inputs made again and the database built by it with models at gates 1
# and 4 for every read type, every mock sample in ~/Mocks timed at the database's gates and at --sparse_first 0, 4s and
# 4g (measure.sh), and each level's profiles scored against the truth (collect.py). Logs in ~/c4test/logs.
# usage: run_all.sh <commit>     (WSL; 4 cores)
set -u
COMMIT=$1
HERE=$(cd "$(dirname "$0")" && pwd)
M=$HERE/../2026-10-13-mock-local-test
W=$HOME/c4test
SRC=$W/protal
mkdir -p $W/logs
echo "start $(date +%T)"
# The clone first, so that the inputs (which need only its scripts) are made while it builds.
[ -d $SRC/.git ] || git clone --quiet /mnt/c/Users/hildebra/Documents/locDev/protal $SRC
(cd $SRC && git fetch --quiet origin && git checkout --quiet $COMMIT) || { echo "no commit $COMMIT"; exit 1; }
if [ ! -e $W/inputs/release/VERSION.txt ] || ! ls $W/inputs/genomes/*.f*a* > /dev/null 2>&1; then
  bash $M/make_inputs.sh $SRC $W/inputs > $W/logs/make_inputs.log 2>&1 &
  INPUTS=$!
fi
bash $M/setup_tools.sh $COMMIT $SRC > $W/logs/setup_tools.log 2>&1 || { echo "build failed"; tail -5 $W/logs/setup_tools.log; exit 1; }
echo "tools built $(date +%T): $(tail -2 $W/logs/setup_tools.log | head -1)"
if [ -n "${INPUTS:-}" ]; then wait $INPUTS || { echo "inputs failed"; tail -5 $W/logs/make_inputs.log; exit 1; }; fi
echo "inputs done $(date +%T)"
if [ ! -f $W/db/protal_db/database.protal ]; then
  bash $HERE/../2026-10-13-consensus-allele-identity/build_db.sh $SRC $W/inputs $W/db --model-gates pe:1/4,pb:1/4,ont:1/4 \
      > $W/logs/build.log 2>&1 || { echo "database build failed"; tail -5 $W/logs/build.log; exit 1; }
fi
echo "database built $(date +%T)"
bash $HERE/measure.sh $SRC $W/db/protal_db $W/measure > $W/logs/measure.log 2>&1
echo "measured $(date +%T)"
label=$($SRC/build/protal --version | tail -1 | sed -E 's/^\S+ v([0-9.]+) \(commit (.{7}).*/\1_\2/')
$HOME/micromamba/envs/protal-db-build/bin/python -I $HERE/collect.py $W/measure $HOME/Mocks $label $W/db/protal_db \
    $HERE/results > $W/logs/collect.log 2>&1 || { echo "collect failed"; tail -5 $W/logs/collect.log; exit 1; }
echo "all done $(date +%T)"
