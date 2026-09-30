#!/bin/bash
# E7: output-format consistency on real protal outputs, ART -l/-f overrides, protal on a flattened SE map.
set -u
W=~/audit6/simulator
SIM=~/strain-build/bin/simulate_metagenomes
MINI=~/strain-build/mini_db/gtdb_r226
cd $W
echo "=== column counts: header vs rows (E1 protal output, sample s_1)"
P=$W/e1/prot/profiles
for f in $P/s_1.profile $P/s_1.profile.log $P/s_1.profile.gene.log $P/s_1.profile.genes.log $P/s_1.profile.truth_annotated $W/e1/prot/misc/s__Mockella_alpha.statistics.tsv; do
  echo "$(basename $f): $(awk -F'\t' '{print NF}' $f | sort | uniq -c | tr '\n' ' ')"
done
echo "--- .profile.log: GeneCov columns and ones that are 0 in every row"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i ~ /^GeneCov[0-9]+$/){c[i]=$i; n++}; print n" GeneCov columns: "$16" .. "$(15+n); next}
            {for(i in c) if($i!=0) nz[i]=1} END{s=""; for(i in c) if(!(i in nz)) s=s" "c[i]; print "always 0 in s_1 and s_2:"s}' $P/s_1.profile.log $P/s_2.profile.log
echo "--- statistics.tsv ends with a newline?"; tail -c 1 $W/e1/prot/misc/s__Mockella_alpha.statistics.tsv | od -c | head -1; cat $W/e1/prot/misc/s__Mockella_alpha.statistics.tsv; echo "|EOF"
echo "--- genes.log header, first row"; head -2 $P/s_1.profile.genes.log | cut -f1-8
echo "--- gene.log header, first row"; head -2 $P/s_1.profile.gene.log | cut -f1-8

echo "=== ART overrides -l 100 and -f 1 in --extra_art_args: manifest vs reads"
GT=$MINI/simulation/genomes.tsv
for extra in "-l 100" "-f 1"; do
  rm -rf e7sim
  $SIM --genome_table $GT -o e7sim --species_per_sample 2 --total_read_pairs 4000 --seed 4 -t 2 --extra_art_args "$extra" > /dev/null 2> e7sim.err; echo "[$extra] exit $?"
  bash $(dirname "$0")/count_reads.sh e7sim
  echo "read length in R1: $(zcat e7sim/reads/sample_1_R1.fq.gz | awk 'NR%4==2{print length($0)}' | sort -u | tr '\n' ' ')"
done

echo "=== protal on a map flattened by protal_map_utils (SE row: SECOND '-')"
mkdir -p e7map; cd e7map
zcat $W/e1/sim/reads/s_1_R1.fq.gz | head -4000 > se.fq
zcat $W/e1/sim/reads/s_1_R1.fq.gz | head -4000 > pe_1.fq
zcat $W/e1/sim/reads/s_1_R2.fq.gz | head -4000 > pe_2.fq
printf "#INPUT_DIR\t$W/e7map\n#OUTPUT_DIR\t$W/e7map/out\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\nse1\tse.fq\t-\tse1\npe1\tpe_1.fq\tpe_2.fq\tpe1\n" > orig.map
python3 ~/strain-build/src/scripts/protal_map_utils flatten --map orig.map > flat.map; cat flat.map
timeout 300 ~/strain-build/bin/protal --db ~/strain-build/mini_db/protal_db/database.protal --map flat.map -t 2 --no_strains > flat.log 2> flat.err; echo "protal on flat.map exit $?"
grep -i "se1\|not exist\|cannot\|error\|single" flat.log flat.err | grep -v "^\s*$" | head -8
