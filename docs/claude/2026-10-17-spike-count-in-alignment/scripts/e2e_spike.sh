#!/usr/bin/env bash
# The e2e SpikeCheckTest alone, verbose, its work folder kept; then the spike lines, spike genes and spike timers of
# the aligned run and the --profile_only run.
cd ~/spikecount/work || exit 1
export PROTAL_TEST_DB=data/mini_db/protal_db PROTAL=build/protal SIMULATE=build/simulate_metagenomes PROTAL_TEST_KEEP=1
export TMPDIR=$HOME/spikecount/tmp; rm -rf "$TMPDIR"; mkdir -p "$TMPDIR"
taskset -c 0-3 python3 tests/e2e/test_protal_e2e.py -v SpikeCheckTest 2>&1 | tail -8
w=$(ls -d "$TMPDIR"/protal_e2e_SpikeCheckTest_* | head -1)
for d in aligned only ahead; do
    echo "== $d"
    cat "$w/$d/spiky.profile.spike_genes.tsv" 2>/dev/null
    grep "spike genes" "$w/$d/misc/spiky_runtime.tsv"
done
grep -h "^Sample spiky: " "$w"/aligned/*.log 2>/dev/null | head -3
"$PROTAL" --version 2>&1 | head -1
