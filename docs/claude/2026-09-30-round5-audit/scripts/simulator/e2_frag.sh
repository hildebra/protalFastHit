#!/bin/bash
# E2: fragmented assemblies. One genome cut into contigs of a fixed size; ART asked for 2000 pairs
# of each via a crafted manifest (--from_manifest replays exactly the requested read_pairs).
set -u
W=~/audit6/simulator
SIM=~/strain-build/bin/simulate_metagenomes
SRC=~/strain-build/mini_db/gtdb_r226/genomic_files_reps/gtdb_genomes_reps_r226/database/GCF/999/001/001/GCF_999001001.1_genomic.fna.gz
mkdir -p $W/e2/genomes
cd $W/e2
TAX="d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella alpha"
printf "sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\trelative_abundance\tfastq_r1\tfastq_r2\tfasta_path\tart_seed\n" > manifest.tsv
for size in 100000 5000 2000 1000 500 400 300 200 140; do
  name=frag$size
  zcat $SRC | grep -v '>' | tr -d '\n' | awk -v s=$size -v n=$name '{L=length($0); for(i=1;i<=L;i+=s){c++; print ">"n"_contig"c; print substr($0,i,s)}}' > genomes/$name.fna
  len=$(grep -v '>' genomes/$name.fna | tr -d '\n' | wc -c)
  nc=$(grep -c '>' genomes/$name.fna)
  echo "$name contigs=$nc len=$len"
  printf "x\t$name\tMockella alpha\t$TAX\t$len\t2000\t\t\t\t\t$W/e2/genomes/$name.fna\t12345\n" >> manifest.tsv
done
rm -rf sim
$SIM --from_manifest manifest.tsv --output_dir sim -t 2 --seed 1 > sim.log 2> sim.err
echo "exit $?"
grep -v "ART command\|ART override" sim.err | sort | uniq -c | head -20
bash $W/../simulator/count_reads.sh sim 2>/dev/null || bash $(dirname "$0")/count_reads.sh sim
