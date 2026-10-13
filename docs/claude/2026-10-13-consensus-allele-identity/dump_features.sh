#!/bin/bash
# Re-profile one run's long-read SAMs with --profile_truth, so that protal writes every taxon's features
# (<profile>.truth_annotated), on the run's database (../2026-10-13-ecoli-hifi/dump_features.sh for any run).
# usage: dump_features.sh <protal> <database folder> <run folder> <out folder>     (WSL; 4 cores)
set -eu
PROTAL=$1; DB=$2; RUN=$3; OUT=$4
HERE=$(cd "$(dirname "$0")" && pwd)
T=$HERE/../../../scripts/mock_benchmark/truth
mkdir -p $OUT
for spec in atcc1003_hifi_SRR11606871:atcc_msa1003:pb d6331_hifi_SRR13128012:d6331:pb d6331_ont_ERR15459768:d6331:ont; do
  IFS=: read -r sample truth type <<< "$spec"
  sam=$(find $RUN -name "$sample.sam.zst")
  python3 -I $HERE/../2026-10-13-ecoli-hifi/truth_taxids.py $T/$truth.tsv $DB/genome2tiid.tsv > $OUT/$sample.truth
  taskset -c 0-3 $PROTAL --db $DB/database.protal --profile_only "$sam" --prefix $sample -o $OUT --read_type $type \
      --profile_truth $OUT/$sample.truth -t 4 --no_strains > $OUT/$sample.log 2>&1
  echo "$sample: $(ls $OUT/$sample.profile.truth_annotated 2>/dev/null || echo missing)"
done
