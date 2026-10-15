#!/bin/bash
# Profile the paired-end mock samples' SAMs of the database-gates run (docs/claude/2026-10-14-mock-sparse-gates, run
# 0.8.4_5369c99_dbgates-gbm: the alignment does not depend on the model kind) again with --profile_truth, so that
# protal writes every taxon's features and probability (<profile>.truth_annotated), with the database's default
# (forest) models, as the scored run.
# usage: dump_features.sh <protal> <database folder> <out folder> [extra protal options ...]     (WSL; 4 cores)
set -eu
PROTAL=$1; DB=$2; OUT=$3; shift 3
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=$HOME/Mocks/0.8.4_5369c99_dbgates-gbm
T=$HERE/../../../scripts/mock_benchmark/truth
mkdir -p $OUT
for spec in d6331_dna_SRR38158090:d6331 d6331_hku_SRR23926922:d6331 d6331_auckland_SRR33710519:d6331 \
            treichel_ERR14788576:treichel_even70 hcom2_SRR21819229:hcom2_closed hcom1a_SRR24799833:hcom1a_full; do
  IFS=: read -r sample truth <<< "$spec"
  sam=$(find $RUN -name "$sample.sam.zst")
  python3 -I $HERE/../2026-10-13-ecoli-hifi/truth_taxids.py $T/$truth.tsv $DB/genome2tiid.tsv > $OUT/$sample.truth 2> $OUT/$sample.truth.log
  taskset -c 0-3 nice -n 5 $PROTAL --db $DB/database.protal --profile_only "$sam" --prefix $sample -o $OUT --read_type pe \
      --profile_truth $OUT/$sample.truth -t 4 --no_strains --force "$@" > $OUT/$sample.log 2>&1 < /dev/null
  echo "$sample: $(ls $OUT/$sample.profile.truth_annotated 2>/dev/null || echo missing)"
done
