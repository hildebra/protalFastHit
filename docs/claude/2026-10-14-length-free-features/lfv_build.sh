#!/usr/bin/env bash
# Validation of the length-free features (2026-10-14): the small real world (~/realworld) built twice with the working
# tree's protal, simulator and scripts (snapshot in ~/lfv): paired-end reads of 75, 100, 150 and 250 bp (15 samples a point)
# and of 150 bp alone (60 a point, as the strain study's model was trained), both with the new default features and a
# random forest; then each training table trained again with the previous default set. 4 cores, niced.
set -euo pipefail
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
V=~/lfv
mkdir -p $V/bin
rsync -a --delete --exclude '/data' --exclude '/local/' ~/grid-build/src/ $V/src/
cp ~/grid-build/build/protal ~/grid-build/build/simulate_metagenomes $V/bin/
PREVIOUS="normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights+strain_identity"
build() {
  local name=$1 samples=$2 setups=$3 out=$V/$1
  if [ ! -s $out/protal_db/database.protal ]; then
    rm -rf $out
    taskset -c 0-3 nice -n 5 python3 $V/src/scripts/build_gtdb_database.py --inputs ~/realworld --outdir $out \
        --protal $V/bin/protal --simulator $V/bin/simulate_metagenomes -t 4 --seed 1 \
        --samples $samples --read-pairs 2000,10000,50000,200000,1000000 --read-setups $setups \
        --species-per-sample 15-60 --archaea 2 --congeners 0.5:2-5 --read-types pe \
        --test-samples 15 --test-read-pairs 10000,100000,1000000 --test-species-per-sample 15-60 \
        --scenarios none --holdout-clades genus:3 --holdout-max-share 0.1 \
        --error-reads none --evaluation basic --no-binary-check --release 226 --model forest > $out.log 2>&1
  fi
  echo "== $name built $(date +%T)"
  # The build's own trainer command, with the previous default set.
  local command
  command=$(grep -h -m1 -o "machine_learning_cmdline.py --truth-file .*" $out/logs/classifier_training.log $out/console.log 2>/dev/null | head -1 || true)
  echo "$command" > $out.trainer_command.txt
}
build mix 15 75:HS25:250:40,100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50
build b150 60 150:HSXt:350:50
echo "== builds done $(date +%T)"
