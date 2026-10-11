#!/usr/bin/env bash
# The tables' own costs, ref against work, after impl_after.sh: profile-only runs of the pe1m sample's SAM (no index:
# every table loads before the profiling, so the bracket and the peak are the tables'), REPS rounds on 4 pinned cores;
# then callgrind on one thread with --sequential_load, the three variants side by side on cores 1-3, and the loaders'
# inclusive instruction counts.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
until grep -q "AFTER DONE" "$A/impl_after.out" 2> /dev/null; do sleep 15; done
W=$HOME/impl/perf; S=$W/samples; R=$W/tables; reps=${REPS:-3}
REF=$HOME/impl/ref/build/protal; WORK=$HOME/impl/work/build/protal
rm -rf "$R"; mkdir -p "$R/logs"
pin=(taskset -c 0-3 nice -n 5)
secs() { awk '{ t = 0; for (i = 1; i <= NF; i++) { if ($i ~ /^[0-9]+m$/) t += 60 * $i; else if ($i ~ /^[0-9]+s$/) t += $i; else if ($i ~ /^[0-9]+ms$/) t += $i / 1000 } printf "%.3f", t }'; }
mem() { grep -m1 "^Memory after $1:" "$2" | sed -E 's/.*: ([0-9.]+) GB resident, peak ([0-9.]+) GB/\1\/\2/'; }
variant() {  # name: sets bin db
  case $1 in ref_text) bin=$REF db=$W/db_ref ;; work_text) bin=$WORK db=$W/db_ref ;; work_bin) bin=$WORK db=$W/db_bin ;; esac
}

echo "[$(date +%T)] the SAM"
"${pin[@]}" "$WORK" --db "$W/db_bin/database.protal" -1 "$S/pe1m/reads/pe1m_1_R1.fq.gz" -2 "$S/pe1m/reads/pe1m_1_R2.fq.gz" -o "$R/sam" \
    -t 4 --no_profile > "$R/logs/sam.log" 2>&1 || { echo "SAM FAIL"; tail "$R/logs/sam.log"; exit 1; }
sam=$(find "$R/sam" -maxdepth 2 -name "*.sam*" | head -1)
echo "SAM: $sam"

echo "[$(date +%T)] profile-only runs"
printf 'variant\trep\twall_s\tmax_rss_mb\ttables_s\tpreload_s\tprofiling_s\tafter_profiling\n' > "$R/profile_only.tsv"
for rep in $(seq 1 "$reps"); do
  for v in ref_text work_text work_bin; do
    variant "$v"
    out=$R/out; rm -rf "$out"
    log=$R/logs/po.$v.$rep.log
    /usr/bin/time -f "%e %M" -o "$R/time.txt" "${pin[@]}" "$bin" --db "$db/database.protal" --profile_only "$sam" -o "$out" -t 4 --no_qcmsa \
        > "$log" 2>&1 || echo "FAILED $v $rep"
    read -r wall rss < <(tail -1 "$R/time.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$v" "$rep" "$wall" "$((rss / 1024))" \
      "$(grep -m1 '^Loading the taxonomy, models and tables took' "$log" | secs)" "$(grep -m1 '^Preload genomes took' "$log" | secs)" \
      "$(grep -m1 '^Profiling took' "$log" | secs)" "$(mem profiling "$log")" >> "$R/profile_only.tsv"
    # the outputs of the three variants compared (all but misc/ and logs)
    if [ "$rep" = 1 ]; then mkdir -p "$R/profiles"; mv "$out" "$R/profiles/$v"; else rm -rf "$out"; fi
  done
done
cat "$R/profile_only.tsv"
for v in work_text work_bin; do
  a=$R/profiles/ref_text b=$R/profiles/$v bad=0
  while IFS= read -r f; do
    rel=${f#$a/}
    case $rel in misc/*|*.log) continue ;; esac
    cmp -s "$a/$rel" "$b/$rel" || { echo "  differs: $rel"; bad=1; }
  done < <(find "$a" -type f | sort)
  echo "$([ $bad = 0 ] && echo same || echo DIFFERENT): ref_text vs $v ($(find "$a" -type f | grep -v '/misc/' | wc -l) files)"
done

echo "[$(date +%T)] callgrind"
core=1
for v in ref_text work_text work_bin; do
  variant "$v"
  ( taskset -c $core nice -n 5 valgrind --tool=callgrind --callgrind-out-file="$R/callgrind.out.$v" \
        "$bin" --db "$db/database.protal" --profile_only "$sam" -o "$R/cg_$v" -t 1 --no_qcmsa --sequential_load > "$R/logs/cg.$v.log" 2>&1
    callgrind_annotate --inclusive=yes --threshold=100 "$R/callgrind.out.$v" 2> /dev/null > "$R/$v.all.txt" ) &
  core=$((core + 1))
done
wait
names=("PROGRAM TOTALS" "LoadAllGenomes(" "protal::LoadStrainAlleles(" "LoadColumnWeights(" "LoadAncestrySites(" "LoadCongenerGaps("
       "LoadSpeciesNeighbours(" "LoadSpeciesPriors(" "LoadGeneConservation(" "LoadSuspectCopies(" "LoadForeignRates(" "ReadGeneNeighbours("
       "LoadTaxonomy(" "ProfileSample(")
printf 'function' > "$R/callgrind.tsv"; for v in ref_text work_text work_bin; do printf '\t%s' "$v" >> "$R/callgrind.tsv"; done; echo >> "$R/callgrind.tsv"
for n in "${names[@]}"; do
  printf '%s' "$n" >> "$R/callgrind.tsv"
  for v in ref_text work_text work_bin; do
    printf '\t%s' "$(grep -F -m1 "$n" "$R/$v.all.txt" | awk '{ print $1 }')" >> "$R/callgrind.tsv"
  done
  echo >> "$R/callgrind.tsv"
done
cat "$R/callgrind.tsv"
rm -rf "$R"/cg_* "$R/sam"
echo "[$(date +%T)] TABLES DONE"
