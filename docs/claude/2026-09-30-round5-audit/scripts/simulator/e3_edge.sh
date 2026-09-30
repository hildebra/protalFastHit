#!/bin/bash
# E3: simulator edge cases.
set -u
W=~/audit6/simulator
SIM=~/strain-build/bin/simulate_metagenomes
MINI=~/strain-build/mini_db/gtdb_r226
WORLD=~/audit5/world/gtdb_r226/simulation/genomes.tsv
ALPHA=$MINI/genomic_files_reps/gtdb_genomes_reps_r226/database/GCF/999/001/001/GCF_999001001.1_genomic.fna.gz
GAMMA=$(grep GCF_999003001.1 $MINI/simulation/genomes.tsv | cut -f3)
TA="d__Bacteria;p__Simulatota;c__Simulatia;o__Simulales;f__Simulaceae;g__Mockella;s__Mockella alpha"
TG="d__Bacteria;p__Fictota;c__Fictia;o__Fictales;f__Fictaceae;g__Fakibacter;s__Fakibacter gamma"
rm -rf $W/e3; mkdir -p $W/e3/g; cd $W/e3

echo "=== (a) duplicate genome names: DUP = alpha (273960 bp) and DUP = first 50 kb of gamma"
zcat $GAMMA | awk '/^>/{print; next}{s=s $0}END{print substr(s,1,50000)}' | sed 's/_contig1/_short/' > g/gamma50k.fna
head -c 60 g/gamma50k.fna; echo
printf "DUP\t$TA\t$ALPHA\nDUP\t$TG\t$W/e3/g/gamma50k.fna\n" > dup.tsv
$SIM --genome_table dup.tsv -o simdup --species_per_sample 2 --total_read_pairs 4000 --seed 3 -t 2 > /dev/null 2> simdup.err; echo "exit $?"
cut -f1-7 simdup/manifest.tsv
zcat simdup/reads/sample_1_R1.fq.gz | awk 'NR%4==1' | sed 's/^@//; s/-[0-9]*\/1$//' | sort | uniq -c

echo "=== (b) repeated --taxon: g__Mockella:2 then g__Fakibacter:1, 3 species per sample, 6 samples (--test)"
$SIM --genome_table $WORLD -o simtaxon -n 6 --species_per_sample 3 --taxon g__Mockella:2 --taxon g__Fakibacter:1 --seed 5 --test > /dev/null 2> simtaxon.err; echo "exit $?"
tail -n +2 simtaxon/manifest.tsv | awk -F'\t' '{split($4,t,";"); print $1"\t"t[6]"\t"$3}' | sort -u | awk -F'\t' '{g[$1]=g[$1]" "$2":"$3} END{for(k in g) print k"\t"g[k]}' | sort
echo "--- same, as one comma-separated --taxon"
$SIM --genome_table $WORLD -o simtaxon2 -n 6 --species_per_sample 3 --taxon "g__Mockella:2,g__Fakibacter:1" --seed 5 --test > /dev/null 2> simtaxon2.err; echo "exit $?"
tail -n +2 simtaxon2/manifest.tsv | awk -F'\t' '{split($4,t,";"); print $1"\t"t[6]"\t"$3}' | sort -u | awk -F'\t' '{g[$1]=g[$1]" "$2":"$3} END{for(k in g) print k"\t"g[k]}' | sort

echo "=== (c) -rs in --extra_art_args: 2 samples of the same single genome"
printf "A1\t$TA\t$ALPHA\n" > one.tsv
$SIM --genome_table one.tsv -o simrs -n 2 --species_per_sample 1 --total_read_pairs 2000 --seed 1 --extra_art_args "-rs 42" -t 2 > /dev/null 2> simrs.err; echo "exit $?"
cut -f1,2,6,12 simrs/manifest.tsv
for s in 1 2; do zcat simrs/reads/sample_${s}_R1.fq.gz | md5sum; done
$SIM --genome_table one.tsv -o simnors -n 2 --species_per_sample 1 --total_read_pairs 2000 --seed 1 -t 2 > /dev/null 2> simnors.err; echo "exit $?"
cut -f1,2,6,12 simnors/manifest.tsv
for s in 1 2; do zcat simnors/reads/sample_${s}_R1.fq.gz | md5sum; done

echo "=== (d) genome names as paths: 'sub/g1' and '../../escape'"
printf "sub/g1\t$TA\t$ALPHA\n../../escape\t$TG\t$GAMMA\n" > paths.tsv
mkdir -p out_d
$SIM --genome_table paths.tsv -o out_d/sim --species_per_sample 2 --total_read_pairs 2000 --seed 1 -t 2 > /dev/null 2> simpaths.err; echo "exit $?"
cut -f1,2,6 out_d/sim/manifest.tsv
find $W/e3/out_d $W/e3 -maxdepth 2 -name "escape*" -o -maxdepth 2 -name "sub" | sort
ls -la $W/e3/out_d

echo "=== (e) --plot_png with a sample prefix containing \$(...) (--test)"
mkdir -p plot/scripts; cp ~/strain-build/src/scripts/plot_abundances.R plot/scripts/
cd plot
$SIM --genome_table ../one.tsv -o simplot --species_per_sample 1 --total_read_pairs 100 --seed 1 --test --plot_png --sample_prefix 'p$(touch INJECTED)' > plot.out 2> plot.err; echo "exit $?"
ls; ls simplot simplot/manifests simplot/plots 2>/dev/null
tail -3 plot.err
cd ..

echo "=== (f) headerless table whose path contains 'genome', 'tax' and 'file'"
mkdir -p taxa_genomes/files
zcat $ALPHA > taxa_genomes/files/a.fna; zcat $GAMMA > taxa_genomes/files/c.fna
printf "GCF_1\t$TA\t$W/e3/taxa_genomes/files/a.fna\nGCF_2\t$TG\t$W/e3/taxa_genomes/files/c.fna\n" > headerless.tsv
$SIM --genome_table headerless.tsv -o simhl --species_per_sample 5 --total_read_pairs 100 --seed 1 --test > /dev/null 2> simhl.err; echo "exit $?"
cut -f2-4 simhl/manifest.tsv; tail -2 simhl.err

echo "=== (g) --test vs real run: same design?"
$SIM --genome_table $MINI/simulation/genomes.tsv -o simT -n 3 --species_per_sample 2 --strains_per_species 0.5 --total_read_pairs 3000 --seed 11 --test > /dev/null 2>&1
$SIM --genome_table $MINI/simulation/genomes.tsv -o simR -n 3 --species_per_sample 2 --strains_per_species 0.5 --total_read_pairs 3000 --seed 11 -t 1 > /dev/null 2>&1
$SIM --genome_table $MINI/simulation/genomes.tsv -o simR2 -n 3 --species_per_sample 2 --strains_per_species 0.5 --total_read_pairs 3000 --seed 11 -t 2 > /dev/null 2>&1
diff <(cut -f1-8,12 simT/manifest.tsv) <(cut -f1-8,12 simR/manifest.tsv) && echo "test == real design"
for s in 1 2 3; do echo "sample_$s t1 $(md5sum < simR/reads/sample_${s}_R1.fq.gz | cut -c1-12) t2 $(md5sum < simR2/reads/sample_${s}_R1.fq.gz | cut -c1-12)"; done
