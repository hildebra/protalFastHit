#!/bin/bash
# After E-s (the collector failed: no simulator; stop() terminated the background build), fix the
# simulator path and rerun into the same folder, as a user would.
exec > ~/audit6/gtdb/fix_rerun.out 2>&1
G=~/audit6/gtdb
S=$G/src/scripts
PY=~/protal-train/bin/python
echo "protal_db and training_db left by the failed run:"; ls $G/es/protal_db $G/es/training_db
$PY $S/build_gtdb_database.py --inputs $G/dl/inputs --samples 2 --read-pairs 1000 --read-setups 100:HS20:300:40 \
   --species-per-sample 8-12 --archaea 1 --holdout-max-share 0.1 --holdout-clades phylum:1,class:1,family:1,genus:2 \
   -t 2 --evaluation none --outdir $G/es --protal ~/strain-build/bin/protal --simulator ~/strain-build/bin/simulate_metagenomes 2>&1 | tail -3
echo "exit ${PIPESTATUS[0]}"
tail -1 $G/es/training_db_index.log
