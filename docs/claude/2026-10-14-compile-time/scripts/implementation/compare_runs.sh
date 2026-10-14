#!/bin/bash
# The split binary against ab6bf25's: the mini database built by each from the same converted folder (files compared
# byte for byte), then the end-to-end tests run with each on its own database, their folders kept (TMPDIR per binary).
set -o pipefail
C=~/split-cmp
BASE=~/split-base/build
SPLIT=~/split-build/build
SRC=~/split-build/src
rm -rf $C
mkdir -p $C
T="taskset -c 0-3 nice -n 5"
python3 $SRC/scripts/mini_db/simulate_gtdb_release.py --outdir $C/gtdb_r226 --release 226 > $C/simulate.log 2>&1 || { tail $C/simulate.log; exit 1; }
python3 $SRC/scripts/mini_db/gtdb_to_protal_db.py --gtdb $C/gtdb_r226 --outdir $C/stage > $C/convert.log 2>&1 || { tail $C/convert.log; exit 1; }
python3 $SRC/scripts/mini_db/gene_neighbours.py --db $C/stage --genome_table $C/gtdb_r226/simulation/genomes.tsv > $C/gene_neighbours.log 2>&1
for v in base split; do
    bin=$BASE/protal; [ $v = split ] && bin=$SPLIT/protal
    cp -r $C/stage $C/db_$v
    /usr/bin/time -f "build $v: %e s" $T $bin --build --no_profile -t 4 --compress_level 3 --db $C/db_$v \
        --reference $C/db_$v/reference.fna --full_reference $C/db_$v/full_reference.fna > $C/build_$v.log 2>&1 || { echo "build $v failed"; tail -20 $C/build_$v.log; exit 1; }
done
echo "== database files that differ:"
(cd $C/db_base && find . -type f | sort) | while read -r f; do cmp -s "$C/db_base/$f" "$C/db_split/$f" || echo "  $f"; done
echo "== (end of list)"
diff <(sed -E 's/[0-9]+(\.[0-9]+)?(s| s|ms| GB| MB)\b/N\2/g' $C/build_base.log) <(sed -E 's/[0-9]+(\.[0-9]+)?(s| s|ms| GB| MB)\b/N\2/g' $C/build_split.log) | head -20
cd $SRC
for v in base split; do
    bin=$BASE; [ $v = split ] && bin=$SPLIT
    mkdir -p $C/tmp_$v
    start=$(date +%s)
    TMPDIR=$C/tmp_$v PROTAL_TEST_KEEP=1 PROTAL_TEST_DB=$C/db_$v PROTAL=$bin/protal SIMULATE=$bin/simulate_metagenomes PROTAL_TESTS_REQUIRED=1 \
        $T python3 -m unittest tests/e2e/test_protal_e2e.py > $C/e2e_$v.log 2>&1
    echo "e2e $v: exit $?, $(( $(date +%s) - start )) s; $(tail -1 $C/e2e_$v.log)"
    grep -E '^(FAIL|ERROR):' $C/e2e_$v.log | head -10
done
du -sh $C/tmp_base $C/tmp_split
