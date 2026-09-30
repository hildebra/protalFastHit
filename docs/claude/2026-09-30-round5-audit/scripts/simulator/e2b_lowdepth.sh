#!/bin/bash
# E2b: low depth x fragmentation. Requested read_pairs 10..1000 of the 3-contig genome and of the
# same genome cut into 5 kb and 20 kb contigs.
set -u
W=~/audit6/simulator
SIM=~/strain-build/bin/simulate_metagenomes
SRC=~/strain-build/mini_db/gtdb_r226/genomic_files_reps/gtdb_genomes_reps_r226/database/GCF/999/001/001/GCF_999001001.1_genomic.fna.gz
cd $W/e2
size=20000; name=frag$size
zcat $SRC | grep -v '>' | tr -d '\n' | awk -v s=$size -v n=$name '{L=length($0); for(i=1;i<=L;i+=s){c++; print ">"n"_contig"c; print substr($0,i,s)}}' > genomes/$name.fna
TAX="d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella alpha"
printf "sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\trelative_abundance\tfastq_r1\tfastq_r2\tfasta_path\tart_seed\n" > manifest_low.tsv
for rp in 10 30 50 100 300 1000; do
  for g in frag100000 frag20000 frag5000; do
    printf "rp$rp\t$g\tMockella alpha\t$TAX\t273960\t$rp\t\t\t\t\t$W/e2/genomes/$g.fna\t999\n" >> manifest_low.tsv
  done
done
rm -rf simlow
$SIM --from_manifest manifest_low.tsv --output_dir simlow -t 2 --seed 1 > simlow.log 2> simlow.err
echo "exit $?"
bash $(dirname "$0")/count_reads.sh simlow
