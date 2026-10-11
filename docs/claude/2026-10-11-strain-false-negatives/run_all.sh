#!/bin/bash
# The whole study after run_extract.sh: reads, database variants, runs and evaluation, the first three variants first
# so that their results come early. protal and its build at 8d8c89b (WSL ~/b23, Release).
# usage: run_all.sh [protal binary]
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
PROTAL=${1:-$HOME/b23/build/protal}
PY=~/micromamba/envs/protal-db-build/bin/python
S=~/strainfn
[ -s $S/reads/samples.tsv ] || taskset -c 0-3 nice -n 5 $PY -I $HERE/simulate_strain_reads.py $S/world $S/reads 1
echo "== reads done $(date +%T)"
bash $HERE/make_variants.sh $PROTAL reps:20:--strain_alleles,0 a4:20: a16i0:20:--strain_alleles,16,--index_alleles,0
bash $HERE/run_samples.sh $PROTAL reps:reps: a4:a4: a16i0:a16i0:
taskset -c 0-3 nice -n 5 $PY -I $HERE/evaluate.py $S $S/eval_first reps a4 a16i0 > $S/eval_first.txt 2>&1
echo "== first evaluation done $(date +%T)"
bash $HERE/make_variants.sh $PROTAL a4i01:20:--index_alleles,0.01 a4i0:20:--index_alleles,0 a16:20:--strain_alleles,16 \
    a16i0n5:5:--strain_alleles,16,--index_alleles,0
bash $HERE/run_samples.sh $PROTAL a4ns:a4:--no_allele_scores a4i01:a4i01: a4i0:a4i0: a16:a16: a16i0n5:a16i0n5: \
    a4c6:a4:-c,6 a4w16:a4:-w,16,-s,256 a16i0c6:a16i0:-c,6
taskset -c 0-3 nice -n 5 $PY -I $HERE/evaluate.py $S $S/eval_all reps a4 a4ns a4i01 a4i0 a16 a16i0 a16i0n5 a4c6 a4w16 a16i0c6 \
    > $S/eval_all.txt 2>&1
echo "== ALL DONE $(date +%T)"
