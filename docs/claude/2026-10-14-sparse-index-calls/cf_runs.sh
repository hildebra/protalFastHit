#!/usr/bin/env bash
# Counterfactual profiles of the strain world's grid C 4 run (sparse index, g4) and the ungated run (g0): their SAMs with
# fields swapped (swap_tags.pl), re-profiled with --profile_only and the same model as run_samples.sh.
# usage: cf_runs.sh name:base:donor:spec ...   (base/donor = run names g0/g4; spec see swap_tags.pl)
set -eu
S=$(cd "$(dirname "$0")" && pwd)
W=~/strainfn
CF=$W/cf
PROTAL=$CF/protal
MODEL=$W/model_s1/model_pe.xml
mkdir -p $CF/truth
samples=$(cut -f1 $W/reads/samples.tsv)
for s in $samples; do
  [ -s $CF/truth/$s.txt ] || awk -F'\t' -v s=$s '$1 == s {print $4}' $W/reads/truth.tsv > $CF/truth/$s.txt
done
for item in "$@"; do
  IFS=: read -r name base donor spec <<< "$item"
  R=$W/run_cf_$name
  [ -s $R/done ] && { echo "== $name exists"; continue; }
  rm -rf $R; mkdir -p $R/alignments $R/profiles
  for s in $samples; do
    perl $S/swap_tags.pl $W/run_$base/alignments/$s.sam.zst $W/run_$donor/alignments/$s.sam.zst $spec $R/alignments/$s.sam.zst 2>> $R/swap.log &
    while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 0.2; done
  done
  wait
  sams=$(for s in $samples; do printf '%s,' $R/alignments/$s.sam.zst; done); sams=${sams%,}
  truths=$(for s in $samples; do printf '%s,' $CF/truth/$s.txt; done); truths=${truths%,}
  ( cd $R && taskset -c 0-3 nice -n 5 $PROTAL --db $W/db_grid --profile_only "$sams" --profile_truth "$truths" -o $R/profiles -t 4 \
      --no_strains --no_qcmsa --model $MODEL > $R/protal.log 2>&1 )
  echo "== $name ($base, $donor's $spec): $(grep -c 'Write truth to' $R/protal.log) profiles"
  date > $R/done
done
