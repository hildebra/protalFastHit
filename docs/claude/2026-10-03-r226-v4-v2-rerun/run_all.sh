#!/bin/bash
# The comparisons of this report, from the repository root, on the runs extracted into local/ (git-ignored):
#   local/protal0.7.3_r226_v3 (v3, with v3/retrain from docs/claude/2026-10-02-r226-v3-training/retrain.sh),
#   local/r226_v4 (local/protal_r226_v4_training.tgz), local/r226_v2_new (local/protal_r226_v2_training.tgz).
# Trains three pairs of long-read models into local/r226_cross (about 5 minutes on 6 threads) and writes the outputs
# next to this script. Usage: PY=python3 bash docs/claude/2026-10-03-r226-v4-v2-rerun/run_all.sh
set -e
PY=${PY:-python3}
HERE=$(cd "$(dirname "$0")" && pwd)
L=local
$PY "$HERE/compare_runs.py" v3=$L/protal0.7.3_r226_v3 v3-retrain=$L/protal0.7.3_r226_v3/retrain v4=$L/r226_v4 v2-new=$L/r226_v2_new \
    > "$HERE/compare_runs_output.txt"
{
    echo "## pe, v3 against v4"; $PY "$HERE/table_diff.py" $L/protal0.7.3_r226_v3/training/training_data.tsv $L/r226_v4/training/training_data.tsv
    echo "## pe, v4 against v2-new"; $PY "$HERE/table_diff.py" $L/r226_v4/training/training_data.tsv $L/r226_v2_new/training/training_data.tsv
    echo "## ont, v4 against v2-new"; $PY "$HERE/table_diff.py" $L/r226_v4/training/training_data_ont.tsv $L/r226_v2_new/training/training_data_ont.tsv
    echo "## ont test, v4 against v2-new"; $PY "$HERE/table_diff.py" $L/r226_v4/test/training_data_ont.tsv $L/r226_v2_new/test/training_data_ont.tsv
} > "$HERE/table_diff_output.txt" 2>&1
$PY "$HERE/common_rows.py" $L/r226_v4 $L/r226_v2_new > "$HERE/common_rows_output.txt"
# Long-read models on one test set (v2-new's): A v4's tables at 256 leaves (v4's own settings), B v4's tables at 128,
# C v2-new's tables at 128 (v2-new's own model; its test predictions equal the HPC run's).
for run in "A r226_v4 256" "B r226_v4 128" "C r226_v2_new 128"; do
    set -- $run
    mkdir -p $L/r226_cross/$1
    for t in pb ont; do
        $PY scripts/random_forest_cmdline.py --truth-file $L/$2/training/training_data_$t.tsv \
            --output-prefix $L/r226_cross/$1/trained_model_$t --features normalized+adjacency --ntree 64 --maxnodes $3 \
            --seed 1 --threads ${THREADS:-6} --taxonomy $L/$2/internal_taxonomy.dmp --evaluation full --depth-knobs \
            --test-file $L/r226_v2_new/test/training_data_$t.tsv > $L/r226_cross/$1/classifier_training_$t.log 2>&1
    done
done
$PY "$HERE/compare_runs.py" A=$L/r226_cross/A B=$L/r226_cross/B C=$L/r226_cross/C > "$HERE/cross_test_output.txt"
{
    for r in r226_v2_new r226_v4 protal0.7.3_r226_v3; do echo "######## $r"; $PY "$HERE/curve_shrink.py" scripts $L/$r; done
    for r in A B C; do echo "######## cross $r (on v2-new's test set)"; $PY "$HERE/curve_shrink.py" scripts $L/r226_cross/$r $L/r226_v2_new/test; done
} > "$HERE/curve_shrink_output.txt"
