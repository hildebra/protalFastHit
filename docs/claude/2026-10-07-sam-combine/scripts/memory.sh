#!/usr/bin/env bash
# memory.sh PROTAL TAG [N...] - what a --profile_only run over N copies of one dense sample keeps in memory
# (docs/claude/2026-10-07-sam-combine). db900n (900-species world) and w900 (1M pairs from 60 species, ~25x each) of
# the performance reports; the SAM is aligned once (by PROTAL, if missing). Each cohort is profiled with strain MSAs
# and with --no_strains; protal's own "Memory after profiling" line (what the strain stage starts from) and the peak
# RSS (/usr/bin/time) are recorded in $M/TAG.tsv. 4 cores.
set -u
PROTAL=$1; TAG=$2; shift 2
NS=${*:-1 2 4 8 16}
DB=${DB:-$HOME/protal-perf/db900n}
R=$HOME/protal-perf/reads/w900/reads
M=${M:-$HOME/samcombine/mem}
RUN="taskset -c 0-3 nice -n 5"
mkdir -p "$M"
if [ ! -f "$M/aln/w900.sam.zst" ]; then
  $RUN "$PROTAL" --db "$DB" -1 "$R/w900_1_R1.fq.gz" -2 "$R/w900_1_R2.fq.gz" --prefix w900 -o "$M/aln" -t 4 --no_profile \
    > "$M/align.log" 2>&1 || { echo "alignment failed"; tail -20 "$M/align.log"; exit 1; }
fi
out=$M/$TAG.tsv
echo -e "samples\tstrains\tafter_profiling_gb\tpeak_gb\twall_s" > "$out"
for n in $NS; do
  rm -rf "$M/cohort$n"; mkdir -p "$M/cohort$n"
  # Hard links: protal refuses a SAM given twice, and resolves symbolic links to tell.
  for i in $(seq 1 "$n"); do ln "$M/aln/w900.sam.zst" "$M/cohort$n/s$i.sam.zst"; done
  for strains in yes no; do
    o=$M/run_${TAG}_${n}_$strains
    rm -rf "$o"
    extra=(); [ "$strains" = no ] && extra=(--no_strains)
    /usr/bin/time -f "%M %e" -o "$o.time" $RUN "$PROTAL" --db "$DB" --profile_only "$M/cohort$n/*.sam.zst" -o "$o" -t 4 \
      --no_qcmsa "${extra[@]}" > "$o.log" 2>&1
    rc=$?
    after=$(grep -o "Memory after profiling: [0-9.]* GB resident" "$o.log" | grep -o "[0-9.]*" | head -1)
    read -r peak_kb wall < "$o.time"
    echo -e "$n\t$strains\t${after:-NA}\t$(awk -v k="$peak_kb" 'BEGIN{printf "%.3f", k/1048576}')\t$wall" | tee -a "$out"
    [ $rc -eq 0 ] || { echo "exit $rc"; grep -E "Error|rror:" "$o.log" | head -5; }
  done
done
