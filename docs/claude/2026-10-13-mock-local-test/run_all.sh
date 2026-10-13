#!/bin/bash
# The whole local test: the tools first (the clone holds the repository's scripts), then the inputs and the samples'
# downloads at once, the reduced database, protal on the samples and the scores. Logs in ~/mocktest/logs.
# usage: run_all.sh <commit>     (WSL)
set -eu
COMMIT=$1
HERE=$(cd "$(dirname "$0")" && pwd)  # this folder: the test's own scripts
W=$HOME/mocktest
SRC=$W/protal
mkdir -p $W/logs
bash $HERE/setup_tools.sh $COMMIT $SRC > $W/logs/setup_tools.log 2>&1
echo "tools done $(date +%T)"
bash $HERE/make_inputs.sh $SRC $W/inputs > $W/logs/make_inputs.log 2>&1 &
INPUTS=$!
bash $HERE/run_test.sh $SRC - download > $W/logs/download.log 2>&1 &
DOWNLOAD=$!
wait $INPUTS
echo "inputs done $(date +%T)"
bash $HERE/run_build.sh $SRC $W/inputs $W/build > $W/logs/build.log 2>&1
echo "build done $(date +%T)"
wait $DOWNLOAD
echo "download done $(date +%T)"
bash $HERE/run_test.sh $SRC $W/build/protal_db run > $W/logs/run.log 2>&1
bash $HERE/run_test.sh $SRC $W/build/protal_db score > $W/logs/score.log 2>&1
echo "all done $(date +%T)"
