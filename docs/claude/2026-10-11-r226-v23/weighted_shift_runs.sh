#!/bin/bash
# The weighted site shift retested on the small real world (docs/claude/2026-10-10-real-small-world): the build of
# ~/realworld_sp_s1..3 again with this session's protal and scripts (~/bnew: the copy-constructor fix and the prebuilt
# ancestry sites), with and without --weighted-site-shift, three seeds each, on cores 0-3. The in-silico strains as
# those builds had them (--insilico-strains 1 --insilico-multi 0.5, the defaults until 2026-10-11), so that the arm
# without the shift can be compared with them file by file. Run: seed 1 without the shift was byte-identical to
# ~/realworld_sp_s1, so ~/rw_ws_off_s2 and _s3 were made links to ~/realworld_sp_s2 and _s3 before the loop reached them
# (a folder with protal_db/build_metadata.tsv counts as done).
set -u
B=~/bnew
for seed in 1 2 3; do
  for arm in off on; do
    out=~/rw_ws_${arm}_s$seed
    [ -f $out/protal_db/build_metadata.tsv ] && { echo "$out done before"; continue; }
    extra=""; [ $arm = on ] && extra=--weighted-site-shift
    s=$(date +%s)
    taskset -c 0-3 nice -n 5 ~/micromamba/envs/protal-db-build/bin/python $B/src/scripts/build_gtdb_database.py \
      --inputs ~/realworld --outdir $out --protal $B/build/protal --simulator $B/build/simulate_metagenomes -t 4 \
      --seed $seed --samples 60 --read-pairs 2000,10000,50000,200000,1000000 --read-setups 150:HSXt:350:50 \
      --species-per-sample 15-60 --archaea 2 --congeners 0.5:2-5 --read-types pe,pb \
      --long-read-bases 1000000,10000000,100000000 --long-read-samples 24 --test-samples 15 \
      --test-read-pairs 10000,100000,1000000 --test-species-per-sample 15-60 --test-long-read-bases 3000000,30000000 \
      --test-long-read-samples 12 --scenarios none --holdout-clades genus:3 --holdout-max-share 0.1 --error-reads none \
      --evaluation basic --no-binary-check --release 226 --insilico-strains 1 --insilico-multi 0.5 $extra \
      > ~/rw_ws_${arm}_s$seed.log 2>&1
    echo "$out exit $? $(( $(date +%s) - s )) s"
  done
done
echo WS RUNS DONE
