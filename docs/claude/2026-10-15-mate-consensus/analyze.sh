#!/bin/bash
# The tables of the report from run_all.sh's runs (guide against each consensus run named): the calls
# (score_sims.py, and score_mocks.py for the mocks), the records that changed (compare_records.py), the fragments that
# moved (changed_fragments.py, simulated samples) and the features and calls that moved (feature_shift.py, over the
# model's features), into RESULTS.
# usage: analyze.sh <runs folder> <results folder> [run label ...]   (default label: consensus)     (WSL)
set -u
RUNS=$1; RES=$2; shift 2
LABELS=${*:-consensus}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$HERE/../../..
PY=$HOME/micromamba/envs/protal-db-build/bin/python
FIN=$HOME/c4test/db/protal_db
TDB=$HOME/c4test/db/work/training_db
WORK=$HOME/c4test/db/work
REPORT=$HOME/c4test/db/model_logs/trained_model.report.txt  # the pe model's features
mkdir -p $RES
for db in fin tdb; do
  dbdir=$FIN; [ $db = tdb ] && dbdir=$TDB
  $PY -I $HERE/score_sims.py $RUNS/$db $dbdir $RES/sims_$db guide $LABELS > $RES/sims_$db.txt
  for label in $LABELS; do
    specs=(); tables=()
    for dir in $RUNS/$db/guide/*/; do
      point=$(basename $dir)
      for sam in $dir/alignments/*.sam.zst; do
        s=$(basename $sam .sam.zst)
        specs+=("$s:$sam:$RUNS/$db/$label/$point/alignments/$s.sam.zst")
        tables+=("$s:$dir/profiles/$s.profile.truth_annotated:$RUNS/$db/$label/$point/profiles/$s.profile.truth_annotated")
      done
    done
    $PY -I $HERE/compare_records.py $RES/records_${db}_$label.tsv guide $label $dbdir --work $WORK "${specs[@]}" 2> /dev/null
    $PY -I $HERE/changed_fragments.py $RES/changed_${db}_$label.tsv $dbdir $WORK "${specs[@]}" > $RES/changed_${db}_$label.txt
    $PY -I $HERE/feature_shift.py $RES/shift_${db}_$label guide $label --features-from $REPORT "${tables[@]}" > $RES/shift_${db}_$label.txt
  done
done
# The mocks: calls by score_mocks.py, records and features as above.
SAMPLES=$RES/samples_plus.tsv
{ cat $REPO/scripts/mock_benchmark/samples.tsv; cat $REPO/docs/claude/2026-10-13-mock-local-test/samples_extra.tsv; } > $SAMPLES
runs=mc_guide; for label in $LABELS; do runs=$runs,mc_$label; done
$PY -I $REPO/scripts/mock_benchmark/score_mocks.py score --reads $HOME/Mocks --out $RUNS/mocks/score --samples $SAMPLES \
    --runs $runs --no-sylph > $RES/mocks_score.txt 2>&1
cp $RUNS/mocks/score/scores/samples_by_run.tsv $RES/mocks_samples_by_run.tsv
cp $RUNS/mocks/score/scores/summary_by_run.tsv $RES/mocks_summary_by_run.tsv
rm -f $SAMPLES
for label in $LABELS; do
  specs=(); tables=()
  for sam in $RUNS/mocks/guide/alignments/*.sam.zst; do
    s=$(basename $sam .sam.zst)
    specs+=("$s:$sam:$RUNS/mocks/$label/alignments/$s.sam.zst")
    tables+=("$s:$RUNS/mocks/guide/profiles/$s.profile.truth_annotated:$RUNS/mocks/$label/profiles/$s.profile.truth_annotated")
  done
  $PY -I $HERE/compare_records.py $RES/records_mocks_$label.tsv guide $label $FIN "${specs[@]}" 2> /dev/null
  $PY -I $HERE/feature_shift.py $RES/shift_mocks_$label guide $label --features-from $REPORT "${tables[@]}" > $RES/shift_mocks_$label.txt
done
# The aligner's own counts of the vote and of guidance, per run.
logs=()
for f in $RUNS/fin/*/*/protal.log $RUNS/tdb/*/*/protal.log $RUNS/mocks/*/protal.log; do
  [ -f $f ] || continue
  rel=${f#$RUNS/}
  logs+=("${rel%/protal.log}:$f")
done
$PY -I $HERE/aligner_counts.py "${logs[@]}" > $RES/aligner_counts.tsv
echo "done"
