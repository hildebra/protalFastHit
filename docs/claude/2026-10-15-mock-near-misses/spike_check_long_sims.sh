#!/bin/bash
# protal's sample check (--spike_genes auto, the default) on the reduced r226 database's simulated long-read test
# samples, whose SAMs are .sam.gz and which have no gold-standard files (the check needs no truth): each sample's log
# line where a taxon has spike genes (none printed: no taxon has any).
# usage: spike_check_long_sims.sh <protal> <database folder> <out folder>     (WSL; 4 cores)
set -u
PROTAL=$1; DB=$2; OUT=$3
mkdir -p $OUT
for sam in $HOME/c4test/db/work/test/points/{pb,ont}_*/protal/alignments/*.sam.gz; do
  point=$(basename $(dirname $(dirname $(dirname $sam))))
  sample=${point}_$(basename $sam .sam.gz)
  case $point in pb*) type=pb ;; ont*) type=ont ;; esac
  taskset -c 0-3 nice -n 5 $PROTAL --db $DB/database.protal --profile_only "$sam" --prefix $sample -o $OUT \
      --read_type $type -t 4 --no_strains --force > $OUT/$sample.run.log 2>&1 < /dev/null
  printf '%s\t%s\n' "$sample" "$(grep -h '^Sample .*spike mode' $OUT/$sample.run.log | sed 's/^Sample [^:]*: //')"
done
echo "runs: $(ls $OUT/*.run.log | wc -l), profiled: $(grep -l '^Profiling sample' $OUT/*.run.log | wc -l)"
