#!/bin/bash
# E9: two genomes whose gzip FASTAs have the same file name (dirA/x.fna.gz, dirB/x.fna.gz) in one sample.
set -u
W=~/audit6/simulator
SIM=~/strain-build/bin/simulate_metagenomes
MINI=~/strain-build/mini_db/gtdb_r226
ALPHA=$MINI/genomic_files_reps/gtdb_genomes_reps_r226/database/GCF/999/001/001/GCF_999001001.1_genomic.fna.gz
GAMMA=$(grep GCF_999003001.1 $MINI/simulation/genomes.tsv | cut -f3)
TA="d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella alpha"
TG="d__Bacteria;p__Fictota;c__Fictia;o__Fictales;f__Fictaceae;g__Fakibacter;s__Fakibacter gamma"
rm -rf $W/e9; mkdir -p $W/e9/dirA $W/e9/dirB; cd $W/e9
cp $ALPHA dirA/x.fna.gz; cp $GAMMA dirB/x.fna.gz
printf "A\t$TA\t$W/e9/dirA/x.fna.gz\nG\t$TG\t$W/e9/dirB/x.fna.gz\n" > t.tsv
$SIM --genome_table t.tsv -o sim --species_per_sample 2 --total_read_pairs 3000 --seed 2 -t 2 --keep_tmp > /dev/null 2> sim.err; echo "exit $?"
cut -f2,6 sim/manifest.tsv
zcat sim/reads/sample_1_R1.fq.gz | awk 'NR%4==1' | sed 's/^@//; s/_contig.*//' | sort | uniq -c
ls sim/sample_1_tmp
