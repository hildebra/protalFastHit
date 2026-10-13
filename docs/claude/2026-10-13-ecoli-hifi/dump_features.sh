#!/bin/bash
# Re-profile the local mock test's SAMs (docs/claude/2026-10-13-mock-local-test, run 0.8.3_4dfc549) with
# --profile_truth, so that protal writes every taxon's features (<profile>.truth_annotated), on the reduced database.
# usage: dump_features.sh <protal> <out folder>     (WSL; 4 cores)
set -eu
PROTAL=$1; OUT=$2
HERE=$(cd "$(dirname "$0")" && pwd)
DB=$HOME/Mocks/reduced_db_r226
RUN=$HOME/Mocks/0.8.3_4dfc549
T=$HERE/../../../scripts/mock_benchmark/truth
mkdir -p $OUT
for spec in atcc1003_hifi_SRR11606871:atcc_msa1003:pb d6331_hifi_SRR13128012:d6331:pb \
            d6331_ont_ERR15459768:d6331:ont d6331_hku_SRR23926922:d6331:pe; do
  IFS=: read -r sample truth type <<< "$spec"
  sam=$(find $RUN -name "$sample.sam.zst")
  python3 -I $HERE/truth_taxids.py $T/$truth.tsv $DB/genome2tiid.tsv > $OUT/$sample.truth
  taskset -c 0-3 $PROTAL --db $DB/database.protal --profile_only "$sam" --prefix $sample -o $OUT --read_type $type \
      --profile_truth $OUT/$sample.truth -t 4 --no_strains > $OUT/$sample.log 2>&1
  echo "$sample: $(ls $OUT/$sample.profile.truth_annotated 2>/dev/null || echo missing)"
done
