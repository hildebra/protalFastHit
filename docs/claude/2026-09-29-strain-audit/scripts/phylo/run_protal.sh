#!/bin/bash
# usage: run_protal.sh <run_name> <subdir> [extra protal args...]
# Re-runs protal on a run's simulated samples into runs/<run>/<subdir> with extra arguments.
set -euo pipefail
name=$1; sub=$2; shift 2
D=~/audit5/phylo/runs/$name
sed "s#^\#OUTPUT_DIR\t.*#\#OUTPUT_DIR\t$D/$sub#" $D/sim/protal.meta > $D/$sub.meta
rm -rf $D/$sub
PROTAL_QCMSA_SCRIPT=$HOME/audit5/src/scripts/qcmsa.py ~/audit5/bin/protal \
   --db ~/audit5/world/protal_db --map $D/$sub.meta -t 2 "$@" > $D/$sub.log 2>&1
echo "[run_protal] $name/$sub done ($*)"
