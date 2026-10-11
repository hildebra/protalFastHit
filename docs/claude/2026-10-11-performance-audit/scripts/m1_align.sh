#!/usr/bin/env bash
# Alignment A/B on the archaea world: 0.7.8 (its own database of the same genes) and 0.8.1 (the finished database) on the
# same reads, alignment only (--no_profile), 4 threads on 4 pinned cores, REPS rounds with the variants alternated in
# each; 0.8.1 also with its new alignment work switched off one at a time. One line per run in runs.tsv:
# sample variant rep wall_s user_s sys_s max_rss_mb aligning_s and protal's counts line.
set -uo pipefail
W=$HOME/perf9; S=$W/samples; R=${OUT:-$W/runs/align}; reps=${REPS:-3}
mkdir -p "$R/logs" "$R/misc"
pe1=$S/pe1m/reads/pe1m_1_R1.fq.gz pe2=$S/pe1m/reads/pe1m_1_R2.fq.gz
pb=$S/long/pb_100000000.fq.gz ont=$S/long/ont_100000000.fq.gz
# sample|variant|binary|database|options
configs=(
  "pe|v078|v078|$W/db078|"
  "pe|v081|v081|$W/arc/protal_db|"
  "pe|v081_noadapt|v081|$W/arc/protal_db|--adaptive_candidates 0"
  "pe|v081_noalleles|v081|$W/arc/protal_db|--no_allele_scores"
  "pe|v081_nonb|v081|$W/arc/protal_db|--no_gene_neighbours"
  "se|v078|v078|$W/db078|"
  "se|v081|v081|$W/arc/protal_db|"
  "pb|v078|v078|$W/db078|"
  "pb|v081|v081|$W/arc/protal_db|"
  "ont|v078|v078|$W/db078|"
  "ont|v081|v081|$W/arc/protal_db|"
  "ont|v081_nobound|v081|$W/arc/protal_db|--no_indel_bound"
)
[ -s "$R/runs.tsv" ] || printf 'sample\tvariant\trep\twall_s\tuser_s\tsys_s\tmax_rss_mb\taligning_s\tcounts\n' > "$R/runs.tsv"
secs() {  # protal's "took 1m 2s 345ms" -> seconds
  awk '{ t = 0; for (i = 1; i <= NF; i++) { if ($i ~ /^[0-9]+m$/) t += 60 * $i; else if ($i ~ /^[0-9]+s$/) t += $i; else if ($i ~ /^[0-9]+ms$/) t += $i / 1000 } printf "%.3f", t }'
}
for rep in $(seq 1 "$reps"); do
  for c in "${configs[@]}"; do
    IFS='|' read -r sample variant bin db opts <<< "$c"
    tag=$sample.$variant.$rep
    grep -q "^$sample	$variant	$rep	" "$R/runs.tsv" && continue
    out=$R/work/$tag; rm -rf "$out"; mkdir -p "$out"
    case $sample in
      pe) reads=(-1 "$pe1" -2 "$pe2") ;;
      se) reads=(-1 "$pe1" --read_type se) ;;
      pb) reads=(-1 "$pb" --read_type pb) ;;
      ont) reads=(-1 "$ont" --read_type ont) ;;
    esac
    # shellcheck disable=SC2086
    /usr/bin/time -f "%e %U %S %M" -o "$out/time.txt" taskset -c 0-3 nice -n 5 "$W/$bin/build/protal" --db "$db" "${reads[@]}" \
        -o "$out/s" --no_profile -t 4 $opts > "$R/logs/$tag.log" 2>&1
    status=$?
    read -r wall user sys rss < "$out/time.txt"
    aligning=$(grep -m1 '^Aligning reads took' "$R/logs/$tag.log" | secs)
    counts=$(grep -m1 -E '^Sample .*candidate alignments tried' "$R/logs/$tag.log" | cut -d: -f2- | sed 's/^ //')
    [ $status -eq 0 ] || counts="FAILED exit $status"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$variant" "$rep" "$wall" "$user" "$sys" "$((rss / 1024))" "$aligning" "$counts" >> "$R/runs.tsv"
    cp "$out"/s/misc/*_runtime.tsv "$R/misc/$tag.runtime.tsv" 2> /dev/null
    rm -rf "$out"
  done
done
echo DONE
