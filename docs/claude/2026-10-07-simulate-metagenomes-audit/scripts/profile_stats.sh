#!/usr/bin/env bash
# The instruments' read statistics (simulate_metagenomes --illumina_report: mean written quality, Q30 share,
# substitutions, insertions, deletions and Ns per base, by read) of two binaries, 50,000 pairs each, to show that the
# per-pair changes (a ziggurat, geometric gaps for rare events, read 2 on its own stream) keep the model's statistics.
# Usage: profile_stats.sh SIM_BEFORE SIM_AFTER
set -uo pipefail
echo -e "instrument\tread\tbinary\tmean_Q\tQ30\tsubstitutions\tinsertions\tdeletions\tNs"
for setup in HS20:100:350 HS25:125:350 HSXt:150:350 NovaSeq:150:350 MSv3:250:550; do
  IFS=: read -r sequencer length fragment <<< "$setup"
  for label in before after; do
    sim=$1; [ "$label" = after ] && sim=$2
    taskset -c 0-3 nice -n 5 "$sim" --illumina_report 50000 --sequencer "$sequencer" --read_length "$length" \
        --fragment_mean "$fragment" --fragment_stdev 50 --seed 11 |
      awk -F'\t' -v s="$sequencer" -v b="$label" '$1 == "R1" || $1 == "R2" {
          printf "%s\t%s\t%s\t%.3f\t%.4f\t%.6f\t%.7f\t%.7f\t%.6f\n", s, $1, b, $2, $3, $4, $5, $6, $7 }'
  done
done
