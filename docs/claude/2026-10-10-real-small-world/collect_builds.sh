#!/bin/bash
# The three builds' summaries (hold-out, strains, samples, the models' tables, the column weights) into builds.txt.
HERE=$(cd "$(dirname "$0")" && pwd)
for s in 1 2 3; do
  echo "== seed $s (~/realworld_s$s)"
  grep -E 'training database \(training_db|in-silico strains of|simulated from:|training data \(|collected in|trained in' ~/realworld_s$s.log | cut -c21-330
  grep -A5 '^read type' ~/realworld_s$s.log
  grep -h '^Column weights codes\|^Column weights from' ~/realworld_s$s/logs/training_db_index.log | cut -c1-330
done > $HERE/builds.txt
wc -l $HERE/builds.txt
