#!/usr/bin/env bash
# Memory over one protal --build, by phase: trace_build.sh LABEL THREADS [protal args]
# BIN (default ~/buildmem/base/build/protal), PRISTINE (a converted database folder, copied first: the
# build packs database.protal and removes the separate files). Every 0.05 s VmRSS, RssAnon, VmHWM from
# /proc; a sample belongs to the phase of the last line protal printed before it. Output in
# ~/buildmem/runs/LABEL (phases.txt; stdout.ts and mem.tsv in seconds since the start).
set -u
BIN=${BIN:-$HOME/buildmem/base/build/protal}
PRISTINE=${PRISTINE:-$HOME/buildmem/db900.pristine}
label=$1; t=$2; shift 2
out=$HOME/buildmem/runs/$label; rm -rf $out; mkdir -p $out
db=$out/db; cp -a $PRISTINE $db
full=""; for f in $db/full_reference.fna $db/full_reference.fna.zst; do [ -e $f ] && full="--full_reference $f"; done
t0=$EPOCHREALTIME
( $BIN --build --no_profile -t $t --db $db --reference $db/reference.fna $full --compress_level 3 "$@" 2>$out/stderr.log |
  while IFS= read -r line; do printf '%s\t%s\n' "$EPOCHREALTIME" "$line"; done > $out/stdout.raw ) &
sleep 0.1
pid=$(pgrep -n -f "$BIN --build --no_profile -t $t --db $db")
while kill -0 $pid 2>/dev/null; do
  s=$(awk '/^VmRSS/{r=$2} /^RssAnon/{a=$2} /^VmHWM/{h=$2} END{print r"\t"a"\t"h}' /proc/$pid/status 2>/dev/null)
  [ -n "$s" ] && printf '%s\t%s\n' "$EPOCHREALTIME" "$s" >> $out/mem.raw
  sleep 0.05
done
wait
awk -F'\t' -v t0=$t0 '{ printf "%.3f\t%s\n", $1 - t0, $2 }' $out/stdout.raw > $out/stdout.ts
awk -F'\t' -v t0=$t0 -v OFS='\t' '{ $1 = sprintf("%.3f", $1 - t0); print }' $out/mem.raw > $out/mem.tsv
# Phases: protal's timers ("<phase> took ..."); a sample belongs to the first timer line printed after it,
# so the phases follow the build's own order whatever it is.
awk -F'\t' 'NR==FNR { if ($2 ~ / took /) { n++; ts[n]=$1; ph[n]=$2; sub(/ took .*/, "", ph[n]) } next }
  { while (i < n && ts[i+1] < $1) i++; q = (i < n ? ph[i+1] : "(after the last timer)"); if ($2 > peak[q]) peak[q] = $2; if ($3 > anon[q]) anon[q] = $3;
    if (!(q in first)) { first[q] = $1; order[++k] = q } if ($4 > hwm) hwm = $4 }
  END { for (j = 1; j <= k; j++) printf "%8.2f s  %7.0f MB  (anon %7.0f MB)  %s\n", first[order[j]], peak[order[j]]/1024, anon[order[j]]/1024, order[j];
        printf "VmHWM %.0f MB\n", hwm/1024 }' $out/stdout.ts $out/mem.tsv | tee $out/phases.txt
grep -E "took" $out/stdout.ts | cut -f2 | tr '\n' ';' | head -c 1500; echo
