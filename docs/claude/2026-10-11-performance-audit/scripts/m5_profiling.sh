#!/usr/bin/env bash
# The profiling stage of 0.8.1 by itself: each sample aligned once (its SAM kept), then profiled again from it
# (--profile_only) by default and with the tables' groups off (--no_site_weights, --no_gene_neighbours,
# --no_allele_scores), REPS alternated rounds; and the SAM written as .sam.gz against .sam.zst (the header's room).
# 4 threads on 4 pinned cores, --no_strains --no_qcmsa. profiling.tsv: sample variant rep wall and protal's line
# "Profiling sample X took ...".
set -uo pipefail
W=$HOME/perf9; S=$W/samples; R=${OUT:-$W/runs/profiling}; reps=${REPS:-3}; DB=${DB:-$W/arc/protal_db}
mkdir -p "$R/logs" "$R/sams"
pin=(taskset -c 0-3 nice -n 5)
declare -A reads=( [pe]="-1 $S/pe1m/reads/pe1m_1_R1.fq.gz -2 $S/pe1m/reads/pe1m_1_R2.fq.gz"
                   [se]="-1 $S/pe1m/reads/pe1m_1_R1.fq.gz --read_type se"
                   [pb]="-1 $S/long/pb_100000000.fq.gz --read_type pb"
                   [ont]="-1 $S/long/ont_100000000.fq.gz --read_type ont" )
for s in pe se pb ont; do
  ls "$R/sams/$s/"*.sam.zst > /dev/null 2>&1 && continue
  # shellcheck disable=SC2086
  "${pin[@]}" "$W/v081/build/protal" --db "$DB" ${reads[$s]} -o "$R/sams/$s" --no_profile -t 4 > "$R/logs/$s.align.log" 2>&1
done
printf 'sample\tvariant\trep\twall_s\tmax_rss_mb\tprofiling_line\n' > "$R/profiling.tsv"
for rep in $(seq 1 "$reps"); do
  for s in pe se pb ont; do
    sam=$(ls "$R/sams/$s/"*.sam.zst | head -1)
    for variant in default no_site_weights no_gene_neighbours no_allele_scores; do
      case $variant in default) opts=() ;; *) opts=("--$variant") ;; esac
      out=$R/work; rm -rf "$out"; mkdir -p "$out"
      log=$R/logs/$s.$variant.$rep.log
      /usr/bin/time -f "%e %M" -o "$R/time.txt" "${pin[@]}" "$W/v081/build/protal" --db "$DB" --profile_only "$sam" --read_type "$s" \
          -o "$out" -t 4 --no_strains --no_qcmsa "${opts[@]}" > "$log" 2>&1 || echo "FAILED $s $variant"
      read -r wall rss < <(tail -1 "$R/time.txt")
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$variant" "$rep" "$wall" "$((rss / 1024))" \
        "$(grep -m1 '^Profiling sample .* took' "$log" | sed 's/^Profiling sample [^ ]* //')" >> "$R/profiling.tsv"
    done
  done
done
# the SAM's finishing: .sam.gz (records copied behind the header) against .sam.zst (written behind its room)
printf 'sample\tformat\trep\tfinish\tline\n' > "$R/sam_finish.tsv"
for rep in 1 2; do
  for s in pe ont; do
    for fmt in gz zst; do
      out=$R/work; rm -rf "$out"; mkdir -p "$out"
      log=$R/logs/$s.sam_$fmt.$rep.log
      # shellcheck disable=SC2086
      "${pin[@]}" "$W/v081/build/protal" --db "$DB" ${reads[$s]} -o "$out" --sam_format "$fmt" --no_profile -t 4 > "$log" 2>&1
      printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$fmt" "$rep" "$(grep -m1 '^Writing the SAM header and file took' "$log" | sed 's/.*took //')" \
        "$(grep -m1 '^SAM header:' "$log" | cut -c1-200)" >> "$R/sam_finish.tsv"
      rm -rf "$out"
    done
  done
done
echo DONE
