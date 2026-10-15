#!/bin/bash
# The spike genes' effect on calls: every mock SAM of the database-gates run and every simulated test sample's SAM of
# the reduced r226 database, profiled again with a protal that has spike mode, once per --spike_genes mode (auto, the
# default; on; off), each with --profile_truth (features of every taxon), into <out>/<mode>. The calls are compared by
# spike_calls.py. (The first version of this script ran a protal whose rule had no sample check, as "rule" and "keep"
# for --keep_spike_genes.)
# usage: spike_compare.sh <protal> <database folder> <out folder> [modes, default "auto on off"]     (WSL; 4 cores)
set -u
PROTAL=$1; DB=$2; OUT=$3; MODES=${4:-auto on off}
HERE=$(cd "$(dirname "$0")" && pwd)
T=$HERE/../../../scripts/mock_benchmark/truth
RUN=$HOME/Mocks/0.8.4_5369c99_dbgates-gbm
mkdir -p $OUT/truth
profile() {  # sample sam type truth_file
  for mode in $MODES; do
    mkdir -p $OUT/$mode
    taskset -c 0-3 nice -n 5 $PROTAL --db $DB/database.protal --profile_only "$2" --prefix $1 -o $OUT/$mode --read_type $3 \
        --profile_truth $4 -t 4 --no_strains --force --spike_genes $mode > $OUT/$mode/$1.log 2>&1 < /dev/null
  done
}
for spec in d6331_dna_SRR38158090:d6331:pe d6331_hku_SRR23926922:d6331:pe d6331_auckland_SRR33710519:d6331:pe \
            treichel_ERR14788576:treichel_even70:pe hcom2_SRR21819229:hcom2_closed:pe hcom1a_SRR24799833:hcom1a_full:pe \
            atcc1003_hifi_SRR11606871:atcc_msa1003:pb d6331_hifi_SRR13128012:d6331:pb d6331_ont_ERR15459768:d6331:ont \
            strain_mixes_ERR13964287:ecoli_kpneumoniae:ont; do
  IFS=: read -r sample truth type <<< "$spec"
  sam=$(find $RUN -name "$sample.sam.zst")
  python3 -I $HERE/../2026-10-13-ecoli-hifi/truth_taxids.py $T/$truth.tsv $DB/genome2tiid.tsv > $OUT/truth/$sample.truth 2> /dev/null
  profile $sample $sam $type $OUT/truth/$sample.truth
  echo "mock $sample done"
done
for point in $HOME/c4test/db/work/test/points/*; do
  case $(basename $point) in rl*) type=pe ;; pb*) type=pb ;; ont*) type=ont ;; esac
  for sam in $point/protal/alignments/*.sam.zst; do
    [ -e "$sam" ] || continue
    sample=$(basename $sam .sam.zst)
    awk -F'\t' 'NR == FNR { if (NF >= 4) t[$4] = $2; next } ($1 in t) { print t[$1] }' $DB/genome2tiid.tsv \
        $point/sim/protal_goldstd/$sample.profile_truth > $OUT/truth/$sample.truth
    profile $sample $sam $type $OUT/truth/$sample.truth
  done
  echo "sim $(basename $point) done"
done
