#!/bin/bash
# When do missing tools fail? E-j: --protal that does not exist; E-s: a simulator that does not exist;
# E-i: a Python without scikit-learn (the system python3 here); E-f: --inputs with --genome-table.
G=~/audit6/gtdb
S=$G/src/scripts
PY=~/protal-train/bin/python
COMMON="--inputs $G/dl/inputs --samples 2 --read-pairs 1000 --read-setups 100:HS20:300:40 --species-per-sample 8-12 --archaea 1 --holdout-max-share 0.1 --holdout-clades phylum:1,class:1,family:1,genus:2 -t 2 --evaluation none"
case "$1" in
j)
  rm -rf $G/ej; START=$(date +%s)
  $PY $S/build_gtdb_database.py $COMMON --outdir $G/ej --protal /nonexistent/protal --simulator ~/strain-build/bin/simulate_metagenomes 2>&1 | tail -4
  echo "exit ${PIPESTATUS[0]} after $(( $(date +%s) - START )) s; logs: $(ls $G/ej | tr '\n' ' ')"
  ;;
s)
  rm -rf $G/es; START=$(date +%s)
  $PY $S/build_gtdb_database.py $COMMON --outdir $G/es --protal ~/strain-build/bin/protal --simulator /nonexistent/simulate_metagenomes 2>&1 | tail -4
  echo "exit ${PIPESTATUS[0]} after $(( $(date +%s) - START )) s"; tail -3 $G/es/training_data.log
  pgrep -af -- "--db $G/es/" | cut -c1-100
  ;;
i)
  rm -rf $G/ei; START=$(date +%s)
  /usr/bin/python3 -c 'import numpy, pandas; print("system python3: numpy", numpy.__version__, "pandas", pandas.__version__)'
  /usr/bin/python3 $S/build_gtdb_database.py $COMMON --outdir $G/ei --protal ~/strain-build/bin/protal --simulator ~/strain-build/bin/simulate_metagenomes 2>&1 | tail -4
  echo "exit ${PIPESTATUS[0]} after $(( $(date +%s) - START )) s"; tail -3 $G/ei/classifier_training.log
  ;;
f)
  $PY $S/build_gtdb_database.py --inputs $G/dl/inputs --outdir $G/ef --genome-table $G/b2/genomes.tsv 2>&1 | tail -2
  echo "exit ${PIPESTATUS[0]}"
  ;;
esac
