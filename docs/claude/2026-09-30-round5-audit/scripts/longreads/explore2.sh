#!/bin/bash
set -u
cd ~/strain-build/mini_db/gtdb_r226
find . -maxdepth 3 | head -60
echo ---
for d in ~/strain-build/mini_db/gtdb_r226 ~/audit5/world/gtdb_r226; do
  echo "== $d"
  ls $d/simulation | head
  head -5 $d/simulation/marker_positions.tsv
  wc -l $d/simulation/marker_positions.tsv
  head -5 $d/simulation/genomes.tsv
  wc -l $d/simulation/genomes.tsv
  ls $d/simulation/genomes_nonreps | head -3
done
echo ---
cat ~/strain-build/mini_db/protal_db/genome2tiid.tsv | head -20
head -5 ~/strain-build/mini_db/protal_db/gene2geneid.tsv
grep -c '>' ~/strain-build/mini_db/protal_db/full_reference.fna
awk '/^>/{if(n)print n; n=0; next}{n+=length($0)}END{print n}' ~/strain-build/mini_db/protal_db/full_reference.fna | sort -n | tail -3
awk '/^>/{if(n)print n; n=0; next}{n+=length($0)}END{print n}' ~/audit5/world/protal_db/full_reference.fna | sort -n | tail -3
grep -c '>' ~/audit5/world/protal_db/full_reference.fna
head -2 ~/strain-build/mini_db/protal_db/full_reference.fna | cut -c1-100
cat ~/strain-build/mini_db.log
cat ~/strain-build/e2e.log
