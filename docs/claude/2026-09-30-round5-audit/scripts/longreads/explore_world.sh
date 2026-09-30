#!/bin/bash
S=$(dirname "$0")
{
ls -la ~/audit5/world/protal_db
cat ~/audit5/world/gtdb_r226/simulation/divergence.tsv | head -60
cut -f2 ~/audit5/world/gtdb_r226/simulation/genomes.tsv | awk -F';' '{print $6, $7}' | sort | uniq -c | head -40
} > $S/out_world.txt 2>&1
