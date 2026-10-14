#!/bin/bash
# Every mock sample of ~/Mocks/reads/downloaded.tsv timed by the clone's scripts/measure_performance.sh on 4 threads: one
# run at the database's gates (pe 4, pb 16, ont 8, stride), then one at each --test_sparse level: 0 (ungated), 4s
# (--sparse_first 4 --sparse_seeding: the reads kept seeded from the hash pick's cores alone, as an index of only those
# would seed them; the whole index loaded) and 4g (--sparse_first 4 --sparse_first_pick grid: the grid's sparse index
# of level 2 loaded, the reads kept seeded from its entries). Without strain MSAs, as the mock runs.
# usage: measure.sh <protal clone> <database folder> <out folder>     (WSL; 4 cores)
set -u
SRC=$1; DB=$2; OUT=$3
mkdir -p $OUT
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
            { print $c["set"] "\t" $c["sample"] "\t" $c["read_type"] "\t" $c["first"] "\t" $c["second"] }' \
    $HOME/Mocks/reads/downloaded.tsv > $OUT/samples.tsv
while IFS=$'\t' read -r set sample type first second; do
  spec=$type:$first; [ "$second" != "-" ] && spec=$spec:$second
  SPARSE_SHORT="0 4s 4g" SPARSE_LONG="0 4s 4g" REPEATS=1 COHORT=0 PROTAL=$SRC/build/protal PROTAL_ARGS="--no_strains" \
      taskset -c 0-3 nice -n 5 bash $SRC/scripts/measure_performance.sh --threads 4 --test_sparse $OUT/$sample $DB "$spec" \
      > $OUT/$sample.log 2>&1 < /dev/null
  echo "$sample ($set, $type): exit $? $(date +%T)"
done < $OUT/samples.tsv
