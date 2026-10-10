#!/bin/bash
# Choose the taxa, stream their marker genes out of the packed GTDB r226 archives (~/GTDB/r226) and measure the
# ancestry sites on them; usage: run.sh <marker set: bac120|ar53> [seed]. Outputs next to this script (<set>_*),
# the extracted subset in ~/real_ancestry/<set>.
set -eu
SET=$1; SEED=${2:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=~/micromamba/envs/protal-db-build/bin/python
SUB=~/real_ancestry/$SET
mkdir -p $SUB
cd $HERE
taskset -c 0-3 nice -n 5 $PY -I select_taxa.py ~/GTDB/r226 ${SET}_taxa.tsv $SEED $SET | tee ${SET}_selection.txt
taskset -c 0-3 nice -n 5 $PY -I extract_subset.py ~/GTDB/r226 ${SET}_taxa.tsv $SUB $SET | tee -a ${SET}_selection.txt
taskset -c 0-3 nice -n 5 $PY real_ancestry.py $SUB ${SET}_taxa.tsv ${SET}_ancestry 4
echo "== RUN DONE $SET $(date +%T)"
