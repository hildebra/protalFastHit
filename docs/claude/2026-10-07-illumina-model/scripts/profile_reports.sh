#!/usr/bin/env bash
# simulate_metagenomes --illumina_report for every instrument at its usual read setup: the summary lines (mean quality,
# share >= Q30, substitutions, insertions, deletions, Ns, the written qualities' shares) and the mean quality of a few
# cycles against the curve the profile was given ("pub").
# Usage: profile_reports.sh SIMULATOR OUT_DIR [PAIRS]
set -uo pipefail
sim=$1 out=$2 pairs=${3:-20000}
mkdir -p "$out"
for spec in HS20:100:300 HS25:125:300 HSXt:150:350 NovaSeq:150:350 MSv3:250:550 MSv3:300:550; do
  IFS=: read -r prof len frag <<< "$spec"
  tsv=$out/$prof.$len.tsv
  /usr/bin/time -f "%e s" taskset -c 0 "$sim" --illumina_report "$pairs" --sequencer "$prof" --read_length "$len" \
      --fragment_mean "$frag" --fragment_stdev 50 > "$tsv" 2> "$out/time.txt" || { echo "FAILED $spec"; cat "$out/time.txt"; continue; }
  echo "== $spec ($(cat "$out/time.txt"))"; sed -n 2,4p "$tsv" | cut -c1-240
  awk -F'\t' -v L="$len" 'NR>5 && ($1==1 || $1==10 || $1==int(L/2) || $1==int(L*0.8) || $1==L) {printf "  cycle %s: R1 %.1f (pub %.1f)  R2 %.1f (pub %.1f)\n", $1, $2, $4, $3, $5}' "$tsv"
done
