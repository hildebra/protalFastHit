#!/bin/bash
# usage: run_sim.sh <run_name> <make_manifest.py args...>
# Simulates the samples of a scenario and runs protal (defaults, qcmsa on) on them.
set -euo pipefail
name=$1; shift
P=~/audit5/phylo
D=$P/runs/$name
mkdir -p $D
python3 $P/scripts/make_manifest.py --out $D/manifest_in.tsv "$@"
if [ ! -f $D/sim/protal.meta ]; then
  ~/audit5/bin/simulate_metagenomes --from_manifest $D/manifest_in.tsv -o $D/sim --seed ${SIMSEED:-42} -t 2 \
     --protal_metafile $D/protal > $D/sim.log 2>&1
fi
rm -rf $D/protal
PROTAL_QCMSA_SCRIPT=$HOME/audit5/src/scripts/qcmsa.py /usr/bin/time -v ~/audit5/bin/protal \
   --db ~/audit5/world/protal_db --map $D/sim/protal.meta -t 2 \
   > $D/protal.log 2>&1
echo "[run_sim] $name done"
