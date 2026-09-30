#!/usr/bin/env bash
# Build one world's database with protal 0.7 and 0.6.0a from the same converted release, timed
# (wall time, CPU, peak RSS: /usr/bin/time -v in logs/<world>.<step>.time).
#   build.sh WORLD RELEASE_DIR [MODEL_DIR]
# MODEL_DIR: trained_model*.xml of build_gtdb_database.py, added to the 0.7 database (pe, se, pb, ont).
# Writes $BENCH/WORLD/{conv,db07,db07raw,db06}. 0.6 gets the shipped model (random_forest.xml) as
# model.xml, the name it reads; 0.7 keeps it as model_pe.xml unless MODEL_DIR replaces it.
set -euo pipefail
WORLD=$1 REL=$2 MODELS=${3:-}
B=${BENCH:-$HOME/bench07}
SRC7=${SRC7:-$HOME/fix-build/src}
P7=${P7:-$HOME/fix-build/bin/protal}
P6=${P6:-$HOME/protal-0.6.0a/src/build/protal}
T=${T:-6}
W=$B/$WORLD
mkdir -p $W $B/logs
tm() { local step=$1; shift; echo "$(date +%T) $WORLD $step"; /usr/bin/time -v -o $B/logs/$WORLD.$step.time "$@" > $B/logs/$WORLD.$step.log 2>&1; }
build() { local p=$1 db=$2; shift 2; echo $p --build --no_profile -t $T --db $db --reference $db/reference.fna --full_reference $db/full_reference.fna "$@"; }

rm -rf $W/conv $W/db07 $W/db07raw $W/db06
tm convert python3 $SRC7/scripts/mini_db/gtdb_to_protal_db.py --gtdb $REL --outdir $W/conv -t $T
for d in db07 db07raw db06; do cp -r $W/conv $W/$d; done
mv $W/db06/model_pe.xml $W/db06/model.xml
tm build07 $(build $P7 $W/db07)                     # database.protal (zstd level 19), the default
tm build07raw $(build $P7 $W/db07raw --no_compress)  # raw files, as 0.6 writes them
tm build06 $(build $P6 $W/db06)
if [ -n "$MODELS" ]; then
  for t in pe se pb ont; do
    m=$MODELS/trained_model$([ $t = pe ] || echo _$t).xml
    [ -f $m ] && tm add_model_$t $P7 --add_model $m --read_type $t --db $W/db07 -t $T
  done
fi
du -sb $W/db07 $W/db07raw $W/db06 > $B/logs/$WORLD.sizes.txt
ls -l $W/db07 $W/db07raw $W/db06 >> $B/logs/$WORLD.sizes.txt
echo "$(date +%T) $WORLD done"
