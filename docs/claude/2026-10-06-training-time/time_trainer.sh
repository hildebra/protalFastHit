#!/bin/bash
# The trainer as build_gtdb_database.py runs it (--features auto, --evaluation full, --depth-knobs, the test table),
# on the r226 v13 tables, to see where its time goes. READ_TYPE: pe (default), se, pb or ont; FEATURES: auto (default).
RT=${READ_TYPE:-pe}
FEATURES=${FEATURES:-auto}
EVAL=${EVAL:-full}
TAG=${TAG:-$RT-$FEATURES-$EVAL}
case $RT in pe) T=training_data.tsv ;; *) T=training_data_$RT.tsv ;; esac
mkdir -p ~/timing
rsync -a --checksum --no-times /mnt/c/Users/hildebra/Documents/locDev/protal/scripts/ ~/protal-gbm/src/scripts/
V=/mnt/c/Users/hildebra/Documents/locDev/protal/local/v13
cd ~/protal-gbm/src
/usr/bin/time -v ~/micromamba/envs/protal-db-build/bin/python scripts/machine_learning_cmdline.py \
  --truth-file $V/training/$T --test-file $V/test/$T --taxonomy $V/internal_taxonomy.dmp \
  --output-prefix ~/timing/$TAG --features $FEATURES --evaluation $EVAL --depth-knobs --threads ${THREADS:-6} \
  $EXTRA > ~/timing/$TAG.log 2>&1
echo "EXIT $?" >> ~/timing/$TAG.log
