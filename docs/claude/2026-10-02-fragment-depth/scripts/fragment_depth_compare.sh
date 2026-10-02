#!/usr/bin/env bash
# The depth with a fragment's bases once (overlapping mates' overlap once) against both mates' bases, on the operon
# world's paired-end test samples (build ~/opw/b_gn2: 2x100 at 300 +- 40, 2x150 at 350 +- 50, 2x250 at 550 +- 50 bp
# fragments; 1,000 to 100,000 pairs; 6 samples each), profiled from their SAMs (--profile_only) against the build's
# training database with its paired-end model, by the binary before (OLD) and after (NEW) the change.
# usage: fragment_depth_compare.sh OUT OLD NEW [THREADS]
set -euo pipefail
out=$1 old=$2 new=$3 threads=${4:-6}
build=${BUILD:-$HOME/opw/b_gn2}
py=$HOME/micromamba/envs/protal-db-build/bin/python3
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
for arm in old new; do
  dir=$out/$arm
  [ -f "$dir/done" ] && continue
  rm -rf "$dir" && mkdir -p "$dir"
  { echo -e "#OUTPUT_DIR\t$dir"; grep '^#SAMPLEID' "$build/test/profile_all/samples.map"
    grep -v '^#' "$build/test/profile_all/samples.map" | grep -P '^rl\d+_p\d+_s_\d+\t' | \
      awk -F'\t' -v d="$dir" 'BEGIN { OFS = "\t" } { $5 = d "/" $1; $6 = d "/" $1 ".profile"; print }'
  } > "$dir/samples.map"
  bin=$([ $arm = old ] && echo "$old" || echo "$new")
  nice "$bin" --db "$build/training_db" --map "$dir/samples.map" --profile_only -t "$threads" --no_strains --no_qcmsa \
    --model "$build/trained_model.xml" > "$dir/protal.log" 2>&1
  touch "$dir/done"
done
"$py" "$here/fragment_depth_summary.py" "$out"
