#!/usr/bin/env bash
# stk_sim_big.sh: larger samples for the timing runs, of the same 40 genomes as stk_sim.sh's (its pe manifest):
# 3M read pairs, 300 Mb of HiFi reads, 150 Mb of Nanopore reads. On cores 0-3.
set -uo pipefail
W=$HOME/stk; D=$W/data; SIM=$W/ref/build/simulate_metagenomes
TABLE=$HOME/realworld_s1/work/genomes_simulated.tsv
ONT_MODEL=$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model
run() { taskset -c 0-3 nice -n 5 "$@"; }
echo "== paired-end 3M ($(date +%T))"
run $SIM --genome_table $TABLE -o $D/pe3m -n 1 --sample_prefix pe3m_s --total_read_pairs 3000000 --species_per_sample 40 \
  --read_length 150 --seed 7 -t 4 --protal_metafile $D/pe3m/protal.meta > $D/pe3m.log 2>&1; echo "exit $?"
ls -la $D/pe3m/reads
mkdir -p $D/big
printf "sample\tout\tbases\tseed\npb\t$D/big/pb.fq.gz\t300000000\t303\n" > $D/big/samples_pb.tsv
printf "sample\tout\tbases\tseed\nont\t$D/big/ont.fq.gz\t150000000\t404\n" > $D/big/samples_ont.tsv
for name in pb ont; do
  echo "== $name ($(date +%T))"
  if [ $name = pb ]; then setup="hifi:15000:3000:3"; model=""; else setup="qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36"; model="--long_model $ONT_MODEL"; fi
  run $SIM --long_samples $D/big/samples_$name.tsv --long_genomes $D/long/genomes_$name.tsv --long_setup "$setup" $model \
    --long_stats $D/big/stats_$name.tsv -t 4 > $D/big/$name.log 2>&1; echo "exit $?"
  ls -la $D/big/$name.fq.gz
done
echo "SIMDONE $(date +%T)"
