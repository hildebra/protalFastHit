#!/bin/bash
# E1: simulate from the mini DB world; count reads ART wrote per genome against the manifest.
set -u
W=~/audit6/simulator
SIM=~/strain-build/bin/simulate_metagenomes
GT=~/strain-build/mini_db/gtdb_r226/simulation/genomes.tsv
mkdir -p $W/e1
cd $W/e1
for f in $(tail -n +2 $GT | cut -f3); do echo "== $f"; zcat $f | grep '>' ; done > headers.txt
head -20 headers.txt
rm -rf sim
/usr/bin/time -v $SIM --genome_table $GT --output_dir sim --samples 2 --sample_prefix s \
  --total_read_pairs 20000 --species_per_sample 3 --strains_per_species "0.6,0.6" --seed 7 -t 2 \
  --keep_tmp --protal_metafile $W/e1/prot > sim.log 2> sim.err
echo "exit $?"
tail -25 sim.err | grep -v "ART command"
cat sim/manifest.tsv | cut -f1-8,12
ls sim sim/reads sim/s_1_tmp | head -40
cat sim/protal.meta
cat sim/protal_goldstd/s_1.profile_truth
cat sim/abundance_matrix.tsv
cat sim/run_params.tsv
