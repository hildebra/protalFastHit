#!/bin/bash
# The whole build_gtdb_database.py run on the tuning world (765 species) with real-sized genomes to simulate
# from (pad_genomes.sh): the default design (read setups, depths, 20-200 species, strains, long-read depths,
# the test set's design) but SAMPLES training and TEST_SAMPLES test samples per design point; ART and pbsim3
# timed per call (timed_tool.sh).
# usage: build_real_size.sh OUT [THREADS]
# needs: ~/bprof/world/genomes.tsv, the tuning world's release ~/tune/release_p, protal and simulate_metagenomes
# in $BIN and the scripts of the same commit in $SCRIPTS (default ~/build-prof/{build,src/scripts})
set -euo pipefail
out=$1 threads=${2:-4}
here=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-$HOME/build-prof/build}
SCRIPTS=${SCRIPTS:-$HOME/build-prof/src/scripts}
env=$HOME/micromamba/envs/protal-db-build
mkdir -p "$out"
bash "$here/make_wrappers.sh" "$out/bin" "$out/tool_times.tsv"
: > "$out/tool_times.tsv"
export PATH="$out/bin:$env/bin:$PATH"  # the timed art_illumina first; zstd, pbsim from the environment
/usr/bin/time -v -o "$out/build.time" "$env/bin/python3" "$SCRIPTS/build_gtdb_database.py" \
  --gtdb ~/tune/release_p --outdir "$out/db" --genome-table ~/bprof/world/genomes.tsv \
  --protal "$BIN/protal" --simulator "$BIN/simulate_metagenomes" -t "$threads" --seed 1 \
  --samples "${SAMPLES:-1}" --test-samples "${TEST_SAMPLES:-1}" \
  --pbsim "$out/bin/pbsim" --pbsim-models "$env/data" > "$out/build.log" 2>&1
echo BUILD_DONE
