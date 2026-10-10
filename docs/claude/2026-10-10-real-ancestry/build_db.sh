#!/bin/bash
# A protal database of the extracted subset (~/real_ancestry/<set>), built by ~/bidx/build/protal (0d961e1: the pooled
# column weights and their log lines), with the strain alleles and the column weights; usage: build_db.sh <set> [aa|nt],
# the column weights' alignments (--column_weights_alignment; default aa).
# The build log goes next to this script (<set>_build.log), its column-weight and allele lines to <set>_build_lines.txt.
set -eu
SET=$1; ALIGN=${2:-aa}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=~/micromamba/envs/protal-db-build/bin/python
SUB=~/real_ancestry/$SET
DB=~/real_ancestry/${SET}_db_$ALIGN
REPO=$HERE/../../..
rm -rf $DB
taskset -c 0-3 nice -n 5 $PY -I $REPO/scripts/mini_db/gtdb_to_protal_db.py --gtdb $SUB --outdir $DB > $HERE/${SET}_convert.log 2>&1
FULL=$DB/full_reference.fna.zst; [ -s $FULL ] || FULL=$DB/full_reference.fna
taskset -c 0-3 nice -n 5 ~/bidx/build/protal --build --no_profile -t 4 --db $DB --reference $DB/reference.fna \
    --full_reference $FULL --column_weights_alignment $ALIGN > $HERE/${SET}_${ALIGN}_build.log 2>&1
grep -E "^(Column weights|Strain alleles|Index alleles|Congener gaps|Ancestry)" $HERE/${SET}_${ALIGN}_build.log | cut -c1-900 > $HERE/${SET}_${ALIGN}_build_lines.txt
cat $HERE/${SET}_${ALIGN}_build_lines.txt
echo "== BUILD DONE $SET $ALIGN"
