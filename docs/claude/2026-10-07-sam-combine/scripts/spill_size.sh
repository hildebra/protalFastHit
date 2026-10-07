#!/usr/bin/env bash
# spill_size.sh PROTAL N - the largest size the spill folder reaches during a --profile_only run over memory.sh's cohort
# of N samples (docs/claude/2026-10-07-sam-combine), sampled every 0.2 s, and the files it held then.
set -u
PROTAL=$1; N=$2
DB=${DB:-$HOME/protal-perf/db900n}
M=${M:-$HOME/samcombine/mem}
rm -rf "$M/spill_size" "$M/run_spill_size"
taskset -c 0-3 nice -n 5 "$PROTAL" --db "$DB" --profile_only "$M/cohort$N/*.sam.zst" -o "$M/run_spill_size" -t 4 --no_qcmsa \
  --strain_spill "$M/spill_size" > "$M/run_spill_size.log" 2>&1 &
pid=$!
max=0; files=0
while kill -0 $pid 2>/dev/null; do
  size=$(du -sb "$M/spill_size" 2>/dev/null | cut -f1)
  if [ -n "$size" ] && [ "$size" -gt "$max" ]; then max=$size; files=$(ls "$M/spill_size" | wc -l); fi
  sleep 0.2
done
wait $pid
echo "exit $?; largest: $((max / 1000000)) MB in $files file(s); left: $(ls "$M/spill_size" 2>/dev/null | wc -l)"
