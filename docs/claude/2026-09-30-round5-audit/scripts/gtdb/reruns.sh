#!/bin/bash
# E-a: rerun of a completed build (same options); E-d: rerun with another --seed (another holdout)
# on a copy; E-k: the collector against another database reuses the dumps.
T=$(dirname "$0")
G=~/audit6/gtdb
PY=~/protal-train/bin/python
S=$G/src/scripts
case "$1" in
a)
  rm -rf $G/b6; cp -a $G/b2 $G/b6
  bash $T/snap.sh b6 b6_before > /dev/null
  bash $T/build1.sh b6 > $G/b6.out 2>&1
  grep -v '^\s' $G/b6.out | tail -25
  bash $T/snap.sh b6 b6_after
  ;;
d)
  rm -rf $G/b7; cp -a $G/b2 $G/b7
  bash $T/build1.sh b7 --seed 2 > $G/b7.out 2>&1
  grep -v '^\s' $G/b7.out | tail -25
  diff <(cut -f1 $G/b2/heldout_species.txt) <(cut -f1 $G/b7/heldout_species.txt) | head -5
  echo "points simulated: $(grep -c 'simulating' $G/b7/training_data.log) lines 'simulating'; profiling: $(grep -c 'profiling' $G/b7/training_data.log)"
  cat $G/b7/training_data.log
  cat $G/b7/training/parity/parity.txt 2>/dev/null
  cmp $G/b2/training/training_data.tsv $G/b7/training/training_data.tsv && echo "training_data.tsv identical to the seed-1 run's"
  ;;
k)
  rm -rf $G/k; mkdir -p $G/k; cp -a $G/b2/training $G/k/training
  rm $G/k/training/training_data.tsv
  $PY $S/collect_training_data.py --db $G/b2/protal_db --genome_table $G/b2/genomes.tsv -o $G/k/training \
     --protal ~/strain-build/bin/protal --simulator ~/strain-build/bin/simulate_metagenomes --samples 2 \
     --read_pairs 1000,20000 --read_setups 100:HS20:300:40,150:HS25:350:50 -t 2 > $G/k/collect.log 2>&1
  echo "exit $?"; cat $G/k/collect.log
  cmp $G/b2/training/training_data.tsv $G/k/training/training_data.tsv && echo "same table as against training_db (dumps reused, no warning)"
  ;;
esac
