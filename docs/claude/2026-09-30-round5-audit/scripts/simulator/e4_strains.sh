#!/bin/bash
# E4: strain sharing spec: MIN_VCOV per strain? SAMPLE_FRACTION / N_STRAINS exclusive?
set -u
W=~/audit6/simulator
SIM=~/strain-build/bin/simulate_metagenomes
WORLD=~/audit5/world/gtdb_r226/simulation/genomes.tsv
mkdir -p $W/e4; cd $W/e4
printf "Mockella alpha\t0.5\t3\t1\t5\t1,1\n" > share.tsv
for sps in 8 4; do
rm -rf sim$sps
$SIM --genome_table $WORLD -o sim$sps -n 6 --species_per_sample $sps --total_read_pairs 20000 --strain_sharing_file share.tsv --seed 2 --test > /dev/null 2> sim$sps.err; echo "exit $?"
echo "--- species_per_sample $sps: Mockella alpha rows (sample genome read_pairs vcov)"
awk -F'\t' '$3=="Mockella alpha"{print $1"\t"$2"\t"$6"\t"$7}' sim$sps/manifest.tsv
echo "samples with Mockella alpha: $(awk -F'\t' '$3=="Mockella alpha"{print $1}' sim$sps/manifest.tsv | sort -u | wc -l) of 6 (SAMPLE_FRACTION 0.5 -> 3)"
echo "distinct alpha strains: $(awk -F'\t' '$3=="Mockella alpha"{print $2}' sim$sps/manifest.tsv | sort -u | wc -l) (N_STRAINS 3)"
echo "strain rows below MIN_VCOV 5: $(awk -F'\t' '$3=="Mockella alpha" && $7<5' sim$sps/manifest.tsv | wc -l)"
done
