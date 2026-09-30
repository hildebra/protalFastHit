#!/bin/bash
set -uo pipefail
A=$HOME/audit5; W=$A/world
grep -E "tests passed|tests failed" $A/src/build/Testing/Temporary/LastTest.log 2>/dev/null | tail -1
ctest --test-dir $A/src/build 2>&1 | grep -E "tests passed|tests failed"
head -3 $W/gtdb_r226/simulation/genomes.tsv | cut -c1-200; wc -l < $W/gtdb_r226/simulation/genomes.tsv; head -4 $W/gtdb_r226/simulation/divergence.tsv
rm -rf $A/smoke; mkdir -p $A/smoke
$A/bin/simulate_metagenomes --genome_table $W/gtdb_r226/simulation/genomes.tsv -o $A/smoke/sim -n 6 --total_read_pairs 200000 \
  --species_per_sample 5-8 --seed 3 -t 6 --protal_metafile $A/smoke/protal > $A/smoke/sim.log 2>&1; echo sim rc=$?
$A/bin/protal --db $W/protal_db --map $A/smoke/sim/protal.meta -t 6 > $A/smoke/protal.log 2>&1; echo protal rc=$?
ls $A/smoke/protal/strains | head -20; grep -a -i "error\|warning" $A/smoke/protal.log | head -5
for f in $A/smoke/protal/strains/*.raw.msa.fna; do echo "$(basename $f): $(grep -c '>' $f) rows, $(sed -n 2p $f | wc -c) columns"; done
