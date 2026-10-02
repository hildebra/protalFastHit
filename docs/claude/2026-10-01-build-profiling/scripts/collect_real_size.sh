#!/bin/bash
# The training data collection of build_gtdb_database.py (its defaults: 20-200 species, strains 0.3,0.1, two
# archaea, one species of a held-out clade per rank, pe, se, pb, ont) on the real-sized genomes, with ART and
# pbsim3 timed per call; fewer samples and depths than the defaults (SAMPLES, READ_PAIRS, LONG_BASES).
#
# usage: collect_real_size.sh OUT [THREADS]
# needs: ~/bprof/world/genomes.tsv (pad_genomes.sh), the tuning world's training database ~/tune/V3/training_db
# and its heldout_species.txt, internal_taxonomy.dmp; protal and simulate_metagenomes in $BIN (default
# ~/build-prof/build)
set -euo pipefail
out=$1 threads=${2:-4}
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
BIN=${BIN:-$HOME/build-prof/build}
SAMPLES=${SAMPLES:-2}
READ_PAIRS=${READ_PAIRS:-1000,500000}
READ_SETUPS=${READ_SETUPS:-150:HSXt:350:50}
LONG_BASES=${LONG_BASES:-300000,150000000}
SPECIES=${SPECIES:-20-200}
mkdir -p "$out"
bash "$here/make_wrappers.sh" "$out/bin" "$out/tool_times.tsv"
: > "$out/tool_times.tsv"
models=$(find "$HOME/micromamba/envs/protal-db-build" -name "QSHMM-ONT-HQ.model" -printf '%h\n' | head -1)
export PATH="$out/bin:$PATH"  # simulate_metagenomes runs art_illumina from PATH
/usr/bin/time -v -o "$out/collector.time" python3 "$repo/scripts/collect_training_data.py" \
  --db ~/tune/V3/training_db --genome_table ~/bprof/world/genomes.tsv -o "$out/collect" \
  --protal "$BIN/protal" --simulator "$BIN/simulate_metagenomes" --samples "$SAMPLES" \
  --read_pairs "$READ_PAIRS" --read_setups "$READ_SETUPS" --archaea 2 --species_per_sample "$SPECIES" \
  --seed 1 -t "$threads" --taxonomy ~/tune/V3/internal_taxonomy.dmp --congeners 0 --read_types pe,se,pb,ont \
  --long_read_bases "$LONG_BASES" --pbsim "$out/bin/pbsim" --pbsim_models "$models" \
  --strains_per_species 0.3,0.1 --novel_species ~/tune/V3/heldout_species.txt --novel_clades 1 \
  > "$out/collector.log" 2>&1
echo COLLECT_DONE
