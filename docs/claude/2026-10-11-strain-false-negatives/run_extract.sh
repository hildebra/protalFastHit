#!/bin/bash
# Choose the strain-rich world (select_strain_world.py) and stream its genomes' marker genes out of the packed GTDB
# r226 archives (the real-ancestry report's extract_subset.py) into WSL ~/strainfn/world.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
PY=~/micromamba/envs/protal-db-build/bin/python
W=~/strainfn/world
mkdir -p $W
taskset -c 0-3 nice -n 5 $PY -I $HERE/select_strain_world.py ~/GTDB/r226 $W 1 | tee $HERE/selection.txt
cp $W/halves.tsv $HERE/halves.tsv
taskset -c 0-3 nice -n 5 $PY -I $HERE/../2026-10-10-real-ancestry/extract_subset.py ~/GTDB/r226 $W/taxa.tsv $W bac120 | tee -a $HERE/selection.txt
echo "== EXTRACT DONE $(date +%T)"
