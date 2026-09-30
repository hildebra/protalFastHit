#!/bin/bash
set -euo pipefail
P=~/audit5/phylo
W=~/audit5/world
G=$W/gtdb_r226/simulation/genomes.tsv
run() { # acc tiid height seed short
  acc=$1; tiid=$2; h=$3; seed=$4; short=$5
  line=$(awk -F'\t' -v a=$acc '$1==a' $G)
  tax=$(echo "$line" | cut -f2); fa=$(echo "$line" | cut -f3)
  sp=$(echo "$tax" | sed 's/.*s__/s__/; s/ /_/g')
  python3 $P/scripts/evolve.py --genome $fa --accession $acc --taxonomy "$tax" \
     --markers $W/gtdb_r226/simulation/marker_positions.tsv --gene2id $W/protal_db/gene2geneid.tsv \
     --dbref $W/protal_db/full_reference.fna --tiid $tiid --ntips 12 --height $h --seed $seed \
     --outdir $P/strains/$short --refname ${sp}_reference
}
run GCF_999001001.1 5 0.01 11 Malpha
run GCF_999008001.1 1 0.005 12 Cferv
run GCF_999006001.1 4 0.002 13 Tone
