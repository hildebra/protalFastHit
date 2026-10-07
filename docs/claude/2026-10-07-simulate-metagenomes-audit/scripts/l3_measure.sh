#!/usr/bin/env bash
# L3 at benchmark scale: one HiFi sample of 250 Mb from the 400 synthetic genomes (bench.sh's long_samples.tsv and
# long_genomes.tsv in WORK), templates cut at contigs' ends (SIM_CUT, 746ea66) against placed where they fit (SIM_PLACED):
# rounds, time and each genome's share of the bases against its weight's (long_shares.py).
# Usage: l3_measure.sh SIM_CUT SIM_PLACED WORK
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
w=$3
for run in "cut:$1" "placed:$2"; do
  label=${run%%:*} sim=${run#*:}
  out=$w/l3_$label; rm -rf "$out"; mkdir -p "$out"
  sed "s#@OUT@#$out#" "$w/long_samples.tsv" > "$w/l3_samples.tsv"
  /usr/bin/time -f "%e s wall, %U s user" -o "$w/l3_time.txt" taskset -c 0-3 nice -n 5 "$sim" --long_samples "$w/l3_samples.tsv" \
      --long_genomes "$w/long_genomes.tsv" --long_setup hifi:15000:3000:3 --long_stats "$w/l3_stats.tsv" -t 1 > "$w/l3.log" 2>&1
  echo "== $label: $(cat "$w/l3_time.txt"); $(tail -1 "$w/l3_stats.tsv")"
  python3 "$here/long_shares.py" "$out/long1.fq.zst" "$w/long_genomes.tsv" long1
done
