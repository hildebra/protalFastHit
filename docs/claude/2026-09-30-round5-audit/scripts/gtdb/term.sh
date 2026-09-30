#!/bin/bash
# E-b: SIGTERM to build_gtdb_database.py while the finished database builds in the background
# (as `kill PID`, or a scheduler that signals the main process only). Does the build survive?
OUT=~/audit6/gtdb/b3
rm -rf $OUT
PY=~/protal-train/bin/python
S=~/audit6/gtdb/src/scripts
$PY $S/build_gtdb_database.py --inputs ~/audit6/gtdb/dl/inputs --outdir $OUT \
    --protal ~/strain-build/bin/protal --simulator ~/strain-build/bin/simulate_metagenomes -t 2 \
    --samples 2 --read-pairs 1000 --read-setups 100:HS20:300:40 \
    --species-per-sample 8-12 --archaea 1 --holdout-max-share 0.1 \
    --holdout-clades phylum:1,class:1,family:1,genus:2 --evaluation none > $OUT.log 2>&1 &
MAIN=$!
until pgrep -f -- "--build --no_profile -t 2 --db $OUT/protal_db" > /dev/null; do sleep 0.5; done
sleep 3
echo "$(date +%T) builds running:"; pgrep -af -- "--build --no_profile" | grep "$OUT" | cut -c1-120
kill -TERM $MAIN
wait $MAIN; echo "$(date +%T) main exited $?"
sleep 2
echo "$(date +%T) after SIGTERM of the main script, still running:"; pgrep -af -- "--build --no_profile" | grep "$OUT" | cut -c1-120
echo "--- a rerun now, with the orphan still writing to protal_db"
$PY $S/build_gtdb_database.py --inputs ~/audit6/gtdb/dl/inputs --outdir $OUT \
    --protal ~/strain-build/bin/protal --simulator ~/strain-build/bin/simulate_metagenomes -t 2 \
    --samples 2 --read-pairs 1000 --read-setups 100:HS20:300:40 \
    --species-per-sample 8-12 --archaea 1 --holdout-max-share 0.1 \
    --holdout-clades phylum:1,class:1,family:1,genus:2 --evaluation none > $OUT.rerun.log 2>&1 &
RERUN=$!
sleep 8
echo "$(date +%T) builds running during the rerun:"; pgrep -af -- "--build --no_profile" | grep "$OUT" | cut -c1-120
wait $RERUN; echo "$(date +%T) rerun exited $?"
tail -5 $OUT.rerun.log
pkill -f -- "--db $OUT/" ; sleep 1
echo "left: $(pgrep -af -- "--db $OUT/" | wc -l)"
