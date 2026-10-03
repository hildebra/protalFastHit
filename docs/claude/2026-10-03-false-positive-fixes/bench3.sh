#!/usr/bin/env bash
# The benchmark-world pipeline again with the refined suspect-copies rule (new3), then the feature-group ablation on
# its tables and a summary of the models.
set -uo pipefail
B=$HOME/protal-fp
W=$HOME/bench071/world
OUT=$HOME/fp_bench
PY=$HOME/micromamba/envs/protal-db-build/bin/python
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
mkdir -p $OUT/logs
common=(--gtdb $W/release_p --protal $B/build/protal --simulator $B/build/simulate_metagenomes -t 6 --seed 1
        --extra-genomes $W/full/simulation/genomes_nonreps --no-binary-check --read-types pe,se --samples 8
        --read-pairs 1000,20000,200000,1000000:4 --test-samples 4 --test-read-pairs 500,10000,100000,1000000:2
        --evaluation basic)
if [ ! -s $OUT/new3/model_logs/summary.txt ]; then
  echo "$(date +%T) pipeline new3"
  /usr/bin/time -v -o $OUT/logs/new3.time python3 $B/src/scripts/build_gtdb_database.py "${common[@]}" --outdir $OUT/new3 \
    > $OUT/logs/new3.log 2>&1 || { echo "pipeline new3 failed"; tail -30 $OUT/logs/new3.log; exit 1; }
  grep -E "Elapsed|Maximum resident" $OUT/logs/new3.time
fi
cat $OUT/new3/model_logs/summary.txt
grep -h "Suspect copies:" $OUT/new3/index_and_package.log $OUT/new3/training_db_index.log | cut -c1-240
grep -h "suspect_copies\|model_pe_depth_knobs\|classifier_features" $OUT/new3/protal_db/build_metadata.tsv
# The ablation on new3's tables.
P=$OUT/new3
E=$OUT/eval3
mkdir -p $E/tables
for rt in "" _se; do
  ln -sf $P/training/training_data$rt.tsv $E/tables/training$rt.tsv
  ln -sf $P/test/training_data$rt.tsv $E/tables/test$rt.tsv
done
RF=$B/src/scripts/random_forest_cmdline.py
train() {
  local rt=$1 feats=$2 tag=$3
  [ -s $E/$tag$rt.log ] && grep -q "Independent test set" $E/$tag$rt.log && return
  $PY $RF --truth-file $E/tables/training$rt.tsv --output-prefix $E/$tag$rt --features $feats --ntree 64 --maxnodes 512 \
    --seed 1 --threads 6 --evaluation basic --depth-knobs --test-file $E/tables/test$rt.tsv > $E/$tag$rt.log 2>&1
  echo "done $tag$rt $(date +%T)"
}
for rt in "" _se; do
  train "$rt" normalized+adjacency+distance nad
  train "$rt" normalized+adjacency+distance+divergence nad_div
  train "$rt" normalized+adjacency+distance+divergence+unfiltered nad_div_unf
  train "$rt" normalized+adjacency+distance+divergence+unfiltered+priors nad_div_unf_pri
  train "$rt" normalized+adjacency+distance+depth+divergence+unfiltered+priors default
done
$PY /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/3611c5a5-4e2b-44f3-a684-74ae56dc8a7c/scratchpad/ablation_summary.py $E $E/tables 2>&1 | grep -v -E "Warning|^\s+df\[" > $E/ablation_output.txt
cat $E/ablation_output.txt
echo BENCH2DONE $(date +%T)
