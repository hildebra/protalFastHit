#!/bin/bash
# Retrain v31's pe and se GBM models on the build's own tables as built (the default feature set) and with the sample
# tail features (the default set + tail), whose two columns tail_from_tables.py computes from the tables as protal
# computes them (context::SampleComplexityOf: the rows' fragments and identity, the sample's identity). The trainer's
# options are the build's (trained_model*.report.txt line 3) except --threads 4, the feature set, and without the C 1 /
# C 8 extra test tables, which the share archive lacks.
#
# Usage (WSL): retrain_with_tail.sh REPO V31_FOLDER WORK_DIR [MARGIN]
#   REPO: the protal checkout (scripts/ with the tail group, model_features.TAIL_FEATURES)
#   V31_FOLDER: the unpacked share archive (protal0.8.5_r226_v31)
#   WORK_DIR: a new folder for the tables and the models (on the Linux file system); a model whose calls are in
#   WORK_DIR/out already is kept (the as-built models do not read the tail columns: they are the same at any margin)
#   MARGIN: the close thin taxa's identity margin (tail_from_tables.py --margin; default 0.01, protal's)
set -euo pipefail
REPO=$1
A=$2
W=$3
MARGIN=${4:-0.01}
DEFAULT=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights+strain_identity+length_free+read_length
HERE=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$W"/tables/training "$W"/tables/test "$W"/out
for s in training test; do
    for t in "" _se; do
        python "$HERE/tail_from_tables.py" --margin "$MARGIN" add "$A/work/$s/training_data$t.tsv" "$W/tables/$s/training_data$t.tsv"
    done
done
python "$HERE/tail_from_tables.py" --margin "$MARGIN" summary "$W/tables/test/training_data.tsv" "$W/tables/test/training_data_se.tsv" \
    > "$W/out/tail_test_samples.txt"

cd "$W/out"
for t in "" _se; do
    for v in as_built with_tail; do
        name=${v}${t:-_pe}
        features=$DEFAULT
        [ "$v" = with_tail ] && features=$DEFAULT+tail
        if [ -f "$W/out/$name.calls.tsv.gz" ]; then echo "=== $name kept"; continue; fi
        echo "=== $name $(date +%T)"
        /usr/bin/time -v python "$REPO/scripts/machine_learning_cmdline.py" \
            --truth-file "$W/tables/training/training_data$t.tsv" --output-prefix "$W/out/$name" \
            --features "$features" --model gbm --ntree 64 --maxnodes 63 --seed 1 --threads 4 --scenario-weight 0.25 \
            --taxonomy "$A/work/internal_taxonomy.dmp" --evaluation basic \
            --feature-versions "$A/model_logs/feature_versions.tsv" --depth-knobs \
            --test-file "$W/tables/test/training_data$t.tsv" > "$W/out/$name.log" 2> "$W/out/$name.err"
        grep -E "Maximum resident|Elapsed" "$W/out/$name.err" || true
    done
done
echo "=== done $(date +%T)"
