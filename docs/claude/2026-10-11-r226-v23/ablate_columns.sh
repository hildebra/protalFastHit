#!/bin/bash
# The full feature set refitted on tables whose eight ancestry-site columns are blanked (blank_site_columns.py): what
# dropping the ancestry sites' construction would cost in F1 while the other columns of ancestry's neighbours stay.
# Core 3, one thread (the other three run ablate_all.sh); output ~/v23abl/out/<type>_no_site_columns.*
set -u
W=~/v23abl
D=$(dirname $0)
PY=~/micromamba/envs/protal-db-build/bin/python
FULL=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights
for t in ${TYPES:-pe pb}; do
  suf=""; [ $t != pe ] && suf=_$t
  [ -f $W/out/${t}_no_site_columns.calls.tsv.gz ] && continue
  for d in training test; do [ -f $W/tables/$d/training_data$suf.tsv ] || cp /mnt/c/Users/hildebra/Documents/locDev/protal/local/v23/protal0.7.9_r226_v23/work/$d/training_data$suf.tsv $W/tables/$d/; done
  ( cd /tmp && taskset -c 3 $PY -I $D/blank_site_columns.py $W/tables $W/tables_nosite $t )
  s=$(date +%s)
  taskset -c 3 nice -n 5 $PY -I $W/scripts/machine_learning_cmdline.py --truth-file $W/tables_nosite/training/training_data$suf.tsv \
    --output-prefix $W/out/${t}_no_site_columns --features $FULL --model gbm --ntree 64 --maxnodes 63 --seed 1 --threads 1 \
    --scenario-weight 0.25 --taxonomy $W/internal_taxonomy.dmp --evaluation basic --depth-knobs \
    --test-file $W/tables_nosite/test/training_data$suf.tsv > $W/out/${t}_no_site_columns.log 2>&1
  echo "$t no_site_columns exit $? $(( $(date +%s) - s )) s"
done
echo COLUMNS DONE
