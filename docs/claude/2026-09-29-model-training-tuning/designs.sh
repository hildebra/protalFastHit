#!/bin/bash
# The build-and-train workflow on the tuning release in five designs (two lanes of 4 threads), and an
# independent test set profiled against the full database.
T=$HOME/tune; R=/mnt/c/Users/hildebra/Documents/locDev/protal; B=$HOME/audit4/bin
export PATH=$HOME/protal-train/bin:$PATH
STRAINS="--extra-genomes $T/world/simulation/genomes_nonreps"
design() {  # name, options
  local name=$1; shift
  [ -f $T/$name/model_logs/build_metadata.tsv ] && { echo "$name already done"; return; }
  rm -rf $T/$name
  /usr/bin/time -f "$name %e s" -o $T/$name.time python3 $R/scripts/build_gtdb_database.py --gtdb $T/release_p --outdir $T/$name \
      --protal $B/protal --simulator $B/simulate_metagenomes -t 4 --samples 8 \
      --read-pairs 1000,3000,10000,30000,100000 --species-per-sample 10-40 --archaea 2 --seed 1 "$@" > $T/$name.log 2>&1 \
      && echo "$name done $(cat $T/$name.time)" || echo "$name FAILED"
}
testset() {  # name, seed, options
  local name=$1 seed=$2; shift 2
  until [ -f $T/A/model_logs/build_metadata.tsv ]; do sleep 30; done
  python3 $R/scripts/collect_training_data.py --db $T/A/protal_db --genome_table $T/world/simulation/genomes.tsv \
      -o $T/$name --protal $B/protal --simulator $B/simulate_metagenomes --samples 8 \
      --read_pairs 1000,3000,10000,30000,100000 --species_per_sample 10-40 --archaea 2 -t 4 --seed $seed \
      --novel_species $T/unknown.txt --taxonomy $T/A/internal_taxonomy.dmp "$@" > $T/$name.log 2>&1 \
      && echo "$name done" || echo "$name FAILED"
}
( design A --holdout 0; design C --holdout 0.1 $STRAINS; design E --holdout 0.1; design F --holdout 0.1 --congeners 3 $STRAINS ) &
( design B --holdout 0 $STRAINS; testset test 909; testset test_congeners 919 --congeners 3; design D --holdout 0.25 $STRAINS ) &
wait
echo "all designs done"
