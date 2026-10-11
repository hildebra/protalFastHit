#!/bin/bash
# The strain alleles' cost in the alignment (ScoreAlleles / ShiftOf on every candidate on a copy with alleles): one new
# paired-end sample of the small real world (200,000 pairs, seed 7, the build's design otherwise), aligned and profiled
# by protal 8d8c89b (~/b23) on one thread (core 3) under callgrind, with and without --no_allele_scores.
set -u
B=~/b23/build
W=~/fcost/align
DB=~/realworld_sp_s1/work/training_db
mkdir -p $W
if [ ! -f $W/sim/reads/al_s_1_R1.fq.zst ]; then
  taskset -c 3 $B/simulate_metagenomes --genome_table ~/realworld_sp_s1/work/genomes_simulated.tsv -o $W/sim -n 1 \
    --sample_prefix al_s --total_read_pairs 200000 --species_per_sample 15-60 --read_length 150 --sequencer HSXt \
    --fragment_mean 350 --fragment_stdev 50 --seed 7 -t 1 --distribution lognormal --pln_sigma 1.3,2.0 \
    --reads_compression zstd --strains_per_species 0.5,0.2 --pick_random_demand_if_fail --congener_groups 0.5:2-5 \
    --genome_store ~/realworld_sp_s1/work/genome_store > $W/simulate.log 2>&1 || { echo "simulation failed"; exit 1; }
fi
R1=$(ls $W/sim/reads/*_R1.fq.zst | head -1); R2=${R1/_R1/_R2}
OUT=$W/align_costs.tsv
[ -f $OUT ] || echo -e "variant\ttotal_Ir\tprofile_Ir\tstrain_alleles_line" > $OUT
for v in default no_allele_scores; do
  grep -q "^$v	" $OUT && continue
  extra=""; [ $v = no_allele_scores ] && extra=--no_allele_scores
  d=$W/$v; rm -rf $d; mkdir -p $d
  ( cd $d && taskset -c 3 valgrind --tool=callgrind --callgrind-out-file=cg.out $B/protal -1 $R1 -2 $R2 --db $DB -o $d/out -t 1 \
      --no_strains $extra > run.log 2> vg.log )
  tot=$(grep -oP "Collected : \K[0-9]+" $d/vg.log)
  prof=$(callgrind_annotate --inclusive=yes $d/cg.out 2>/dev/null | grep -m1 "protal::ProfileSample(" | awk '{gsub(",", "", $1); print $1}')
  callgrind_annotate --inclusive=yes $d/cg.out 2>/dev/null | head -60 > $d/top.txt
  line=$(grep -m1 "strain alleles:" $d/run.log | cut -c1-300)
  echo -e "$v\t$tot\t$prof\t$line" >> $OUT
  rm -f $d/cg.out
  echo "$v done"
done
echo ALIGN COST RUNS DONE
