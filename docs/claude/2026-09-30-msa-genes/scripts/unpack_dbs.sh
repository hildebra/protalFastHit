#!/bin/bash
# Unpack the unique k-mer tables (and taxonomy) of two realistic databases into ~/genes_study/<name>.
set -u
BIN=$HOME/audit5/bin/protal_P2
for pair in "tuneH2:$HOME/tune/H2/protal_db/database.protal" "db900:$HOME/protal-perf/db900/database.protal" \
            "world:$HOME/audit5/world/protal_db/database.protal"; do
  name=${pair%%:*}; db=${pair#*:}
  out=$HOME/genes_study/$name
  [ -s $out/unique_kmers.tsv ] || [ -s $out/unique_kmers.tsv.zst ] && { echo "$name: already unpacked"; continue; }
  mkdir -p $out
  $BIN --unpack_db --db $db --unpack_dir $out -t 4 > $out/unpack.log 2>&1 || { echo "$name: unpack failed"; tail -3 $out/unpack.log; continue; }
  echo "$name: $(ls $out | tr '\n' ' ')"
done
