#!/bin/bash
# The ceiling of more genomes per species: every genome of the world in the full reference, the test strains too, so
# that a test strain's own copies can be its species' alleles (up to 16 per gene, all indexed). Then the per-species
# table (species_table.py) over the runs.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
PROTAL=${1:-$HOME/b23/build/protal}
PY=~/micromamba/envs/protal-db-build/bin/python
S=~/strainfn
bash $HERE/make_variants.sh $PROTAL oracle:all:--strain_alleles,16,--index_alleles,0
bash $HERE/run_samples.sh $PROTAL oracle:oracle:
taskset -c 0-3 nice -n 5 $PY -I $HERE/evaluate.py $S $S/eval_oracle reps a4 a16 a16i0 oracle > $S/eval_oracle.txt 2>&1
taskset -c 0-3 nice -n 5 $PY -I $HERE/species_table.py $S reps a4 a16i0 oracle > $S/species_table.txt 2>&1
echo "== ORACLE DONE $(date +%T)"
