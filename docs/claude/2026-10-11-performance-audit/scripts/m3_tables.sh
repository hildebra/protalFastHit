#!/usr/bin/env bash
# The start-up's tables one by one: the experiment build (tim: 0.8.1 + AuditTimed) with --sequential_load, each table's
# wall time, resident memory before and after and its own peak (VmHWM reset before it); and 0.8.1 itself with the
# default concurrent start-up, for its "Loading the taxonomy, models and tables" bracket. The pe100k sample, profiled
# (so that every table loads), --no_strains --no_qcmsa, 4 threads on 4 pinned cores, REPS rounds.
set -uo pipefail
W=$HOME/perf9; S=$W/samples; R=${OUT:-$W/runs/tables}; reps=${REPS:-3}; DB=${DB:-$W/arc/protal_db}
mkdir -p "$R/logs"
printf 'rep\ttable\tseconds\trss_before_mb\trss_after_mb\tpeak_mb\n' > "$R/tables.tsv"
printf 'rep\tvariant\twall_s\tmax_rss_mb\tgene_tables_s\tpreload_s\ttables_s\tindex_s\n' > "$R/startup.tsv"
secs() { awk '{ t = 0; for (i = 1; i <= NF; i++) { if ($i ~ /^[0-9]+m$/) t += 60 * $i; else if ($i ~ /^[0-9]+s$/) t += $i; else if ($i ~ /^[0-9]+ms$/) t += $i / 1000 } printf "%.3f", t }'; }
for rep in $(seq 1 "$reps"); do
  for variant in tim_sequential v081_concurrent v081_sequential; do
    case $variant in
      tim_sequential) bin=$W/tim/build/protal; opts=(--sequential_load) ;;
      v081_concurrent) bin=$W/v081/build/protal; opts=() ;;
      v081_sequential) bin=$W/v081/build/protal; opts=(--sequential_load) ;;
    esac
    out=$R/work; rm -rf "$out"; mkdir -p "$out"
    log=$R/logs/$variant.$rep.log
    /usr/bin/time -f "%e %M" -o "$out/time.txt" taskset -c 0-3 nice -n 5 "$bin" --db "$DB" -1 "$S/pe100k_R1.fq.gz" -2 "$S/pe100k_R2.fq.gz" \
        -o "$out/s" -t 4 --no_strains --no_qcmsa "${opts[@]}" > "$log" 2>&1 || echo "FAILED $variant $rep"
    read -r wall rss < "$out/time.txt"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$rep" "$variant" "$wall" "$((rss / 1024))" \
      "$(grep -m1 '^Loading the gene tables took' "$log" | secs)" "$(grep -m1 '^Preload genomes took' "$log" | secs)" \
      "$(grep -m1 '^Loading the taxonomy, models and tables took' "$log" | secs)" "$(grep -m1 '^Load Index took' "$log" | secs)" >> "$R/startup.tsv"
    [ "$variant" = tim_sequential ] && grep '^AUDIT table' "$log" | awk -F'\t' -v r="$rep" '{ print r "\t" $2 "\t" $3 "\t" $6 "\t" $8 "\t" $10 }' >> "$R/tables.tsv"
    rm -rf "$out"
  done
done
# what the tables hold (counts), from the last tim log
grep -E '^(Gene conservation|Suspect copies|Species priors|Species neighbours|Congener gaps|Foreign rates|Column weights|Ancestry sites|Strain alleles|Gene neighbours):' "$R/logs/tim_sequential.$reps.log" > "$R/table_counts.txt"
# the tables' sizes in the database (unpacked by db078.sh)
ls -la "$W/db078" > "$R/table_files.txt" 2> /dev/null
echo DONE
