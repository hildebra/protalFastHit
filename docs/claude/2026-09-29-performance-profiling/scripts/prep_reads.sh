#!/bin/bash
# Simulated paired reads (ART HS25, 2x150, fragment 350+-50), 1M pairs each:
#   prep_reads.sh NAME GENOME_TABLE SPECIES     one sample from simulate_metagenomes
#   prep_reads.sh --mix NAME FROM               every 20th pair of reads/FROM interleaved with 950k pairs
#                                               from 40 Mb of random sequence no database contains: ~5%
#                                               of pairs from database species, as in a real metagenome
#                                               whose marker genes are a few % of each genome
set -e
source "$(dirname "$0")/env.sh"
R=$PERF_DIR/reads; mkdir -p $R; cd $R
export LC_ALL=C
if [ "$1" != "--mix" ]; then
  name=$1; table=$2; species=$3
  rm -rf $R/$name
  /usr/bin/time -f "$name simulate %e s" $SIM --genome_table $table -o $R/$name -n 1 --sample_prefix $name \
     --total_read_pairs 1000000 --species_per_sample $species --strains_per_species 0.3,0.1 --seed 11 -t "$(nproc)" \
     --protal_metafile $R/$name/protal > $R/$name.log 2>&1
  tail -1 $R/$name.log
  ls -la $R/$name/reads
  exit 0
fi
name=$2; from=$3
MAP=$(printf 'ACGT%.0s' $(seq 64))
rm -f bg.fa
for i in $(seq 1 10); do
  echo ">bg_contig_$i" >> bg.fa
  head -c 4000000 /dev/urandom | tr '\000-\377' "$MAP" | fold -w 80 >> bg.fa
done
# 950k pairs * 300 bp / 40 Mb = 7.125x
art_illumina -ss HS25 -i bg.fa -p -l 150 -f 7.125 -m 350 -s 50 -rs 13 -na -q -o bg_ > bg.art.log 2>&1
zcat $R/$from/reads/*_R1.fq.gz | paste - - - - > w1.tsv
zcat $R/$from/reads/*_R2.fq.gz | paste - - - - > w2.tsv
paste w1.tsv w2.tsv > from.pairs.tsv
paste - - - - < bg_1.fq > b1.tsv
paste - - - - < bg_2.fq > b2.tsv
paste b1.tsv b2.tsv > bg.pairs.tsv
awk -F'\t' -v OFS='\t' 'NR==FNR { if (FNR % 20 == 0) on[++n] = $0; next }
  { print; if (FNR % 19 == 0 && k < n) print on[++k] }
  END { while (k < n) print on[++k] }' from.pairs.tsv bg.pairs.tsv > mix.pairs.tsv
mkdir -p $name
cut -f1-4 mix.pairs.tsv | tr '\t' '\n' | pigz -p "$(nproc)" > $name/${name}_R1.fq.gz
cut -f5-8 mix.pairs.tsv | tr '\t' '\n' | pigz -p "$(nproc)" > $name/${name}_R2.fq.gz
echo "$name pairs: $(wc -l < mix.pairs.tsv)"
rm -f w1.tsv w2.tsv b1.tsv b2.tsv from.pairs.tsv bg.pairs.tsv mix.pairs.tsv bg_1.fq bg_2.fq
ls -la $name
