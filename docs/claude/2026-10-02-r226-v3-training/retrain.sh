#!/bin/bash
# The four r226 v3 models trained again on v3's own training and test tables, with the settings of suggestions 1 and 4
# (leaves per read type: 512 pe/se, 128 pb/ont; knob points of at least 6 samples), as build_gtdb_database.py
# would run the trainer otherwise (64 trees, seed 1, normalized+adjacency, full evaluation, depth knobs, the taxonomy).
# Usage: retrain.sh TRAINER_DIR V3_DIR OUT_DIR  (TRAINER_DIR: the scripts/ folder with the trainer to run)
set -e
TRAINER=$1/random_forest_cmdline.py
V3=$2
OUT=$3
PY=${PY:-python3}
mkdir -p "$OUT"
for t in pe se pb ont; do
    case $t in pe) suffix="" leaves=512 ;; se) suffix=_se leaves=512 ;; pb) suffix=_pb leaves=128 ;; ont) suffix=_ont leaves=128 ;; esac
    start=$(date +%s)
    $PY "$TRAINER" --truth-file "$V3/training/training_data$suffix.tsv" --output-prefix "$OUT/trained_model$suffix" \
        --features normalized+adjacency --ntree 64 --maxnodes $leaves --seed 1 --threads "${THREADS:-6}" \
        --taxonomy "$V3/internal_taxonomy.dmp" --evaluation full --depth-knobs \
        --test-file "$V3/test/training_data$suffix.tsv" > "$OUT/classifier_training$suffix.log" 2>&1
    echo "$t: trained in $(( $(date +%s) - start )) s"
done
