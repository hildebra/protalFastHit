#!/usr/bin/env bash
# --flex_indels on ONT reads of the small real world (docs/claude/2026-10-13-flex-indels): fi_ont_sim.sh's samples (the
# first SAMPLES communities of its paired-end test point rl150_p100000, BASES read bases each, at each pbsim3 mean
# ACCURACY), one protal run of all of them per arm (a map; 4 threads on cores 0-3, --no_strains), every record scored
# against its read's true source. Arms (ARMS): plain (no gate), cC (--sparse_first C), cCs (and --sparse_seeding), a tT
# after either for --sparse_first_score T (else the database's, here 14), each with an f last for --flex_indels ont (flex:
# plain with it). Writes, per accuracy, WORK/<acc>/reads.tsv (the arms' read-level counts), WORK/<acc>/runs.tsv (time,
# the reads the gate gave up, the seeding's and alignment's counts) and the logs; TABLES (environment) names the two
# tables' suffix (reads<TABLES>.tsv), so that a second sweep keeps the first's.
#   reads.tsv: records (not secondary), those of the read's own species and genus, the reads, those with a record, the
#   reads of a species the database has (in_db), those with a record of their species (in_db_species).
#   bash fi_ont.sh PROTAL SIMULATE WORK
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
protal=$1 simulate=$2 work=$3
build=${BUILD:-$HOME/realworld_sp_s1}
db=$build/protal_db/database.protal
pe_sim=$build/work/test/points/rl150_p100000/sim
model=${MODEL:-$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model}
samples=${SAMPLES:-6} bases=${BASES:-20000000}
read -r -a accuracies <<< "${ACCURACIES:-0.97 0.90}"
read -r -a arms <<< "${ARMS:-plain flex c1 c1f c8 c8f c4s c4sf c8s c8sf}"
pin=(taskset -c 0-3 nice -n 5)
mkdir -p "$work"
arm_args() {  # plain, flex, cC, cCs, either with tT, and each with a trailing f
  local a=$1
  [ "$a" = flex ] && { echo "--sparse_first 0 --flex_indels ont"; return; }
  [ "$a" = plain ] && { echo "--sparse_first 0"; return; }
  [[ $a =~ ^c([0-9]+)(s?)(t([0-9]+))?(f?)$ ]] || { echo "unknown arm $a" >&2; exit 1; }
  echo "--sparse_first ${BASH_REMATCH[1]}${BASH_REMATCH[2]:+ --sparse_seeding}${BASH_REMATCH[4]:+ --sparse_first_score ${BASH_REMATCH[4]}}${BASH_REMATCH[5]:+ --flex_indels ont}"
}
tables=${TABLES:-}
for acc in "${accuracies[@]}"; do
  dir=$work/$acc; mkdir -p "$dir"
  rows=$dir/rows.tsv
  bash "$here/fi_ont_sim.sh" "$simulate" "$pe_sim" "$build/work/genome_store" "$model" "$dir/sim" "ont_a${acc#0.}" "$bases" \
      "$samples" "$acc" > "$rows" || { echo "simulation at $acc failed"; exit 1; }
  # Each read's genome (its sample's k-th genome) and whether the database has its species: read, genome, in_db.
  if [ ! -s "$dir/read_sources.tsv" ]; then
    while IFS=$'\t' read -r sample reads _; do
      zstdcat "$reads" | awk -v s="$sample" 'NR % 4 == 1 { print s "\t" substr($1, 2) }'
    done < "$rows" | awk -F'\t' '
      FILENAME == ARGV[1] { if (FNR > 1) place[$1, n[$1]++] = $2; next }
      FILENAME == ARGV[2] { if (FNR == 1) { for (i = 1; i <= NF; i++) col[$i] = i; next } truth[$col["genome"]] = $col["taxonomy"]; next }
      FILENAME == ARGV[3] { in_db[$4]; next }
      { k = $2; sub(/^g/, "", k); sub(/x_.*/, "", k); g = place[$1, k + 0]; print $2 "\t" g "\t" (truth[g] in in_db ? 1 : 0) }' \
      "$dir/sim/tmp/genomes.tsv" "$pe_sim/manifest.tsv" "$build/protal_db/genome2tiid.tsv" - > "$dir/read_sources.tsv"
  fi
  runs=$dir/runs$tables.tsv reads_table=$dir/reads$tables.tsv
  echo -e "arm\twall_s\tuser_s\tseeding_s\tgiven_up\tkmers\tfound\tcells\tseeds\tpaired\tdropped\tanchors\ttried\trecords" > "$runs"
  echo -e "arm\trecords\tspecies\tgenus\treads\treads_with_record\tin_db\tin_db_species" > "$reads_table"
  for arm in "${arms[@]}"; do
    out=$dir/arms/$arm; rm -rf "$out"; mkdir -p "$out"
    map=$dir/$arm.map
    { echo -e "#OUTPUT_DIR\t$out"; echo -e "#SAMPLEID\tPREFIX\tFIRST\tSECOND\tPROFILE_TRUTH\tREAD_TYPE"
      awk -F'\t' '{ printf "%s\t%s\t%s\t-\t%s\tont\n", $1, $1, $2, $3 }' "$rows"; } > "$map"
    log=$dir/$arm.log
    # shellcheck disable=SC2046
    "${pin[@]}" /usr/bin/time -f "TIME %e %U %S" "$protal" --db "$db" --map "$map" -t 4 --no_strains $(arm_args "$arm") > "$log" 2>&1
    echo "$acc $arm exit $?"
    seeding=$(cat "$out"/misc/*_runtime.tsv 2> /dev/null | awk -F'\t' '$1 == "Seeding" { s += $2 } END { printf "%.2f", s }')
    awk -v a="$arm" -v seeding="$seeding" '
      /^TIME / { wall = $2; user = $3 }
      / sparse first \(.*reads given up/ { for (i = 1; i <= NF; i++) if ($(i + 1) == "of" && $(i + 3) == "reads") g += $i }
      / seeding: [0-9]+ k-mers looked up/ { for (i = 1; i <= NF; i++) {
          if ($(i + 1) == "k-mers") k += $i; if ($(i + 1) == "in" && $(i + 2) == "the") f += $i
          if ($(i + 1) == "flex" && $(i + 2) == "cells") c += $i; if ($(i + 1) == "seeds,") s += $i
          if ($(i + 1) == "of" && $(i + 2) == "them") p += $i; if ($(i + 1) == "lookups" && $(i + 2) == "dropped") d += $i
          if ($(i + 1) == "anchors") an += $i } }
      /candidate alignments tried/ { for (i = 1; i <= NF; i++) if ($(i + 1) == "candidate" && $(i + 2) == "alignments") t += $i }
      /records written/ { for (i = 1; i <= NF; i++) if ($(i + 1) == "records" && $(i + 2) == "written") r += $i }
      END { printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", a, wall, user, seeding, g, k, f, c, s, p, d, an, t, r }' "$log" >> "$runs"
    # Every record (not secondary) against its read's source: its taxon's lineage (genome2tiid.tsv) against the
    # source genome's (the manifest's taxonomy).
    for sam in "$out"/*.sam.zst; do zstdcat "$sam" | awk -F'\t' '!/^@/ && $3 != "*" && int($2 / 256) % 2 == 0 { print $1 "\t" $3 }'; done |
      awk -F'\t' -v a="$arm" '
        FILENAME == ARGV[1] { src[$1] = $2; db[$1] = $3; reads++; in_db += $3; next }
        FILENAME == ARGV[2] { if (FNR == 1) { for (i = 1; i <= NF; i++) col[$i] = i; next } truth[$col["genome"]] = $col["taxonomy"]; next }
        FILENAME == ARGV[3] { lineage[$2] = $4; next }
        function genus(l) { sub(/;s__.*/, "", l); return l }
        { t = $2; sub(/_.*/, "", t); g = src[$1]; records++; with[$1]
          if (lineage[t] == truth[g]) { species++; right[$1] }
          if (genus(lineage[t]) == genus(truth[g])) genus_ok++ }
        END { nw = 0; for (r in with) nw++; ns = 0; for (r in right) ns += db[r]
              printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", a, records, species, genus_ok, reads, nw, in_db, ns }' \
        "$dir/read_sources.tsv" "$pe_sim/manifest.tsv" "$build/protal_db/genome2tiid.tsv" - >> "$reads_table"
  done
  echo "== accuracy $acc"; column -t "$runs"; column -t "$reads_table"
done
