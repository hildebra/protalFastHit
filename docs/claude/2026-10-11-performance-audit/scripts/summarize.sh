#!/usr/bin/env bash
# Medians of the audit's measurements: m1 (alignment A/B), m3 (start-up tables), m5 (profiling), m4 (cohort).
W=$HOME/perf9/runs
median() { sort -g | awk '{ v[NR] = $1 } END { if (NR == 0) print "NA"; else if (NR % 2) print v[(NR + 1) / 2]; else print (v[NR / 2] + v[NR / 2 + 1]) / 2 }'; }
if [ -s "$W/align/runs.tsv" ]; then
  echo "== alignment A/B (medians over reps): sample variant n wall_s user_s max_rss_mb aligning_s"
  for key in $(awk -F'\t' 'NR > 1 { print $1 "|" $2 }' "$W/align/runs.tsv" | sort -u); do
    s=${key%%|*}; v=${key#*|}
    rows=$(awk -F'\t' -v s="$s" -v v="$v" 'NR > 1 && $1 == s && $2 == v' "$W/align/runs.tsv")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$v" "$(echo "$rows" | wc -l)" "$(echo "$rows" | cut -f4 | median)" \
      "$(echo "$rows" | cut -f5 | median)" "$(echo "$rows" | cut -f7 | median)" "$(echo "$rows" | cut -f8 | median)"
  done
  echo "== counts lines (rep 1)"
  awk -F'\t' 'NR > 1 && $3 == 1 { print $1 "\t" $2 "\t" substr($9, 1, 330) }' "$W/align/runs.tsv"
fi
if [ -s "$W/tables/tables.tsv" ]; then
  echo "== start-up tables (tim, sequential; medians): table seconds rss_after-rss_before peak-rss_before"
  for t in $(awk -F'\t' 'NR > 1 { print $2 }' "$W/tables/tables.tsv" | awk '!seen[$0]++'); do
    rows=$(awk -F'\t' -v t="$t" 'NR > 1 && $2 == t' "$W/tables/tables.tsv")
    printf '%s\t%s\t%s\t%s\n' "$t" "$(echo "$rows" | cut -f3 | median)" \
      "$(echo "$rows" | awk -F'\t' '{ print $5 - $4 }' | median)" "$(echo "$rows" | awk -F'\t' '{ print $6 - $4 }' | median)"
  done
  echo "== start-up brackets"; cat "$W/tables/startup.tsv"
  cat "$W/tables/table_counts.txt" 2> /dev/null | cut -c1-250
fi
if [ -s "$W/profiling/profiling.tsv" ]; then
  echo "== profiling (--profile_only; medians): sample variant wall_s max_rss_mb, and rep 1's line"
  for key in $(awk -F'\t' 'NR > 1 { print $1 "|" $2 }' "$W/profiling/profiling.tsv" | awk '!seen[$0]++'); do
    s=${key%%|*}; v=${key#*|}
    rows=$(awk -F'\t' -v s="$s" -v v="$v" 'NR > 1 && $1 == s && $2 == v' "$W/profiling/profiling.tsv")
    printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$v" "$(echo "$rows" | cut -f4 | median)" "$(echo "$rows" | cut -f5 | median)" \
      "$(echo "$rows" | head -1 | cut -f6 | cut -c1-260)"
  done
  echo "== SAM finishing"; cat "$W/profiling/sam_finish.tsv" | cut -c1-260
fi
if [ -s "$W/cohort/cohort.tsv" ]; then echo "== cohort"; cat "$W/cohort/cohort.tsv"; fi
