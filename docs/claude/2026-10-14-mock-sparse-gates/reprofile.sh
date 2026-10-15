#!/bin/bash
# Every mock sample's SAM of the database-gates run (score_mocks.py's run <label>_dbgates-gbm: the alignment does not
# depend on the model kind) profiled again with protal's default models (forests), so that its .profile.log lists every
# taxon with reads and its probability, as the forest run that collect.py scored.
# usage: reprofile.sh <protal> <database folder> <run folder> <out folder>     (WSL; 4 cores)
set -u
PROTAL=$1; DB=$2; RUN=$3; OUT=$4
mkdir -p $OUT
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } { print $c["sample"] "\t" $c["read_type"] }' \
    $HOME/Mocks/reads/downloaded.tsv | while IFS=$'\t' read -r sample type; do
  sam=$(find $RUN -name "$sample.sam.zst" | head -1)
  [ -n "$sam" ] || { echo "$sample: no SAM"; continue; }
  taskset -c 0-3 nice -n 5 $PROTAL --db $DB --profile_only "$sam" --prefix $sample -o $OUT --read_type $type -t 4 \
      --no_strains --force > $OUT/$sample.log 2>&1 < /dev/null
  echo "$sample: exit $? $(ls $OUT/$sample.profile.log 2> /dev/null)"
done
