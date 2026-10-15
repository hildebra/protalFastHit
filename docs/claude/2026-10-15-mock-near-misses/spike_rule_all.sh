#!/bin/bash
# The spike-gene rule's prototype (spike_rule.py) on every mock SAM of the database-gates run and on the reduced r226
# database's simulated test samples (their database-gate SAMs), with each sample's truth.
# usage: spike_rule_all.sh <reference.map of the database> <genome2tiid.tsv> <out folder> [spike_rule.py options]
set -u
MAP=$1; G2T=$2; OUT=$3; shift 3
HERE=$(cd "$(dirname "$0")" && pwd)
PY="python3 -I"
mkdir -p $OUT
T=$HERE/../../../scripts/mock_benchmark/truth
RUN=$HOME/Mocks/0.8.4_5369c99_dbgates-gbm
# mocks: sample, truth table
for spec in d6331_dna_SRR38158090:d6331 d6331_hku_SRR23926922:d6331 d6331_auckland_SRR33710519:d6331 \
            treichel_ERR14788576:treichel_even70 hcom2_SRR21819229:hcom2_closed atcc1003_hifi_SRR11606871:atcc_msa1003 \
            d6331_hifi_SRR13128012:d6331 d6331_ont_ERR15459768:d6331 strain_mixes_ERR13964287:ecoli_kpneumoniae; do
  IFS=: read -r sample truth <<< "$spec"
  sam=$(find $RUN -name "$sample.sam.zst")
  $PY $HERE/../2026-10-13-ecoli-hifi/truth_taxids.py $T/$truth.tsv $G2T > $OUT/$sample.truth 2> /dev/null
  echo "== $sample"
  zstdcat $sam | $PY $HERE/spike_rule.py $MAP --truth $OUT/$sample.truth "$@"
done
# simulated test samples: lineages in sim/protal_goldstd/<sample>.profile_truth
for point in $HOME/c4test/db/work/test/points/*; do
  for sam in $point/protal/alignments/*.sam.zst; do
    [ -e "$sam" ] || continue
    sample=$(basename $sam .sam.zst)
    awk -F'\t' 'NR == FNR { if (NF >= 4) t[$4] = $2; next } ($1 in t) { print t[$1] }' $G2T \
        $point/sim/protal_goldstd/$sample.profile_truth > $OUT/$sample.truth
    echo "== sim $sample"
    zstdcat $sam | $PY $HERE/spike_rule.py $MAP --truth $OUT/$sample.truth --quiet "$@"
  done
done
