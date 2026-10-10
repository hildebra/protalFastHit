#!/bin/bash
# Three builds of the small real world (seeds 1-3: other species and genera held out, other samples), each with its
# refits (ablate_all.sh), into ~/realworld_s<seed> and ~/realworld_s<seed>_abl. usage: run_all.sh
HERE=$(cd "$(dirname "$0")" && pwd)
for seed in 1 2 3; do
  OUT=~/realworld_s$seed
  [ -s $OUT/protal_db/database.protal ] || { rm -rf $OUT; bash $HERE/run_build.sh ~/realworld $OUT $seed > $OUT.log 2>&1 || { echo "build $seed failed"; exit 1; }; }
  echo "== build $seed done $(date +%T)"
  grep -A6 "^read type" $OUT.log
  bash $HERE/ablate_all.sh $OUT > ${OUT}_abl.log 2>&1
  echo "== refits $seed done $(date +%T)"
done
echo "== ALL SEEDS DONE"
