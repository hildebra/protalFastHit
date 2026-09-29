#!/bin/bash
# A protal database from an extracted (synthetic or real) GTDB release:
#   prep_db.sh NAME GTDB_RELEASE_DIR
# e.g. prep_db.sh db900 world/  (world from simulate_gtdb_release.py --lineages, see README)
set -e
source "$(dirname "$0")/env.sh"
db=$1; gtdb=$2
cd $PERF_DIR
rm -rf $PERF_DIR/$db
/usr/bin/time -f "$db convert %e s %M KB" python3 $PROTAL_SRC/scripts/mini_db/gtdb_to_protal_db.py --gtdb $gtdb --outdir $PERF_DIR/$db > $db.convert.log 2>&1
/usr/bin/time -f "$db build %e s %M KB" $PERF_DIR/build-rel/protal --build --no_profile -t "$(nproc)" --db $PERF_DIR/$db \
  --reference $PERF_DIR/$db/reference.fna --full_reference $PERF_DIR/$db/full_reference.fna > $db.build.log 2>&1
tail -1 $db.convert.log; tail -1 $db.build.log
grep -E "Values size|^keys|^values" $db.build.log
ls -la $PERF_DIR/$db
