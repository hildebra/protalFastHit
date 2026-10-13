#!/usr/bin/env bash
# What seeding the kept reads from the gate's cores alone (--sparse_seeding) changes against seeding them from all
# (docs/claude/2026-10-13-sparse-seeding), on one read file: the arms alternated for REPS rounds, 4 threads on cores 0-3,
# --no_strains, the profile kept. An arm is "plain" (no gate), "cC" (--sparse_first C) or "cCs" (--sparse_first C
# --sparse_seeding); ARMS (environment), default "plain c8 c4s c8s c16s"; REF (default c8) the arm the others are
# compared with read by read. Per run (WORK/runs.tsv): wall and CPU seconds, the aligning time, the reads given up, the
# k-mers the kept reads looked up, the candidate alignments tried, the records written. Per arm, from the last round
# (WORK/compare.tsv): the reads with a record (not secondary), the (read, gene) and (read, taxon) pairs of those records,
# each against REF's (shared, lost, gained), and the profile's species against REF's (lost, gained, Bray-Curtis).
#   bash ss_compare.sh PROTAL DATABASE WORK REPS -- READ ARGUMENTS (e.g. -1 HIFI.fq.zst --read_type pb)
#   COMPARE_ONLY=1 bash ss_compare.sh PROTAL DATABASE WORK 0     (the comparison again, from WORK's saved files)
set -uo pipefail
protal=$1 db=$2 work=$3 reps=$4; shift 4
[ "${1:-}" = "--" ] && shift
reads=("$@")
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
read -r -a order <<< "${ARMS:-plain c8 c4s c8s c16s}"
ref=${REF:-c8}
arm_args() {  # plain -> --sparse_first 0; cC -> --sparse_first C; cCs -> and --sparse_seeding
  case $1 in
    plain) echo "--sparse_first 0" ;;
    c*s) local c=${1#c}; echo "--sparse_first ${c%s} --sparse_seeding" ;;
    c*) echo "--sparse_first ${1#c}" ;;
  esac
}
seconds() {  # "1m 2s 345ms" or "2s 345ms" or "345ms" -> seconds
  awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) s += substr(v, 1, length(v) - 2) / 1000;
         else if (v ~ /m$/) s += 60 * substr(v, 1, length(v) - 1); else if (v ~ /s$/) s += substr(v, 1, length(v) - 1) } printf "%.3f", s }'
}
# COMPARE_ONLY=1 (environment): no runs, the comparison of the arms' saved records and profiles again.
[ -n "${COMPARE_ONLY:-}" ] && reps=0
[ "$reps" -gt 0 ] && echo -e "arm\trep\twall_s\tuser_s\taligning_s\tgiven_up\tkmers_kept\ttried\trecords" > "$work/runs.tsv"
for rep in $(seq 1 "$reps"); do
  for arm in "${order[@]}"; do
    dir=$work/run_$arm; rm -rf "$dir"; log=$work/$arm.$rep.log
    # shellcheck disable=SC2046
    "${pin[@]}" /usr/bin/time -f "TIME %e %U %S" "$protal" --db "$db" "${reads[@]}" --prefix s -o "$dir" -t 4 --no_strains \
        $(arm_args "$arm") > "$log" 2>&1
    read -r wall user < <(grep "^TIME" "$log" | awk '{ print $2, $3 }')
    aligning=$(grep -m1 "Aligning reads took" "$log" | sed 's/.*took //' | seconds)
    given=$(grep -m1 "^Sample s sparse first" "$log" | sed -E 's/.*: ([0-9]+) of [0-9]+ (pairs|reads) given up.*/\1/')
    kmers=$(grep -m1 "^Sample s seeding: " "$log" | sed -E 's/.*seeding: ([0-9]+) k-mers.*/\1/')
    tried=$(grep -m1 "candidate alignments tried" "$log" | sed -E 's/.* ([0-9]+) candidate alignments tried.*/\1/')
    written=$(grep -m1 "records written" "$log" | sed -E 's/.* ([0-9]+) records written.*/\1/')
    echo -e "$arm\t$rep\t$wall\t$user\t$aligning\t${given:-0}\t$kmers\t$tried\t$written" | tee -a "$work/runs.tsv"
    # The last round's records and profile, for the comparison.
    if [ "$rep" = "$reps" ]; then
      zstdcat "$dir/s.sam.zst" | awk -F'\t' '!/^@/ && int($2 / 256) % 2 == 0 { print $1 "\t" $3 }' | sort -u > "$work/$arm.pairs"
      cp "$dir/s.profile" "$work/$arm.profile" 2> /dev/null
    fi
  done
done
echo "== medians"
for arm in "${order[@]}"; do
  awk -F'\t' -v a="$arm" '$1 == a { w[++n] = $3; u[n] = $4; al[n] = $5 }
       END { asort(w); asort(u); asort(al); m = int((n + 1) / 2); printf "%s\twall %.2f\tuser %.2f\taligning %.3f\n", a, w[m], u[m], al[m] }' "$work/runs.tsv"
done
# Read by read against REF: (read, gene) pairs as written; reads; (read, taxon) with the taxon the reference's prefix
# before "_" (protal's references are <taxid>_<gene id>).
count() { wc -l < "$1" | tr -d ' '; }
echo -e "arm\treads\treads_shared\treads_lost\treads_gained\tgene_pairs\tgene_shared\tgene_lost\tgene_gained\ttaxon_pairs\ttaxon_shared\ttaxon_lost\ttaxon_gained\tspecies\tspecies_lost\tspecies_gained\tbray_curtis" > "$work/compare.tsv"
for arm in "${order[@]}"; do
  for kind in reads gene taxon; do
    case $kind in
      reads) cut -f1 "$work/$arm.pairs" | sort -u > "$work/$arm.reads"; cut -f1 "$work/$ref.pairs" | sort -u > "$work/ref.reads" ;;
      gene) cp "$work/$arm.pairs" "$work/$arm.gene"; cp "$work/$ref.pairs" "$work/ref.gene" ;;
      taxon) awk -F'\t' '{ split($2, t, "_"); print $1 "\t" t[1] }' "$work/$arm.pairs" | sort -u > "$work/$arm.taxon"
             awk -F'\t' '{ split($2, t, "_"); print $1 "\t" t[1] }' "$work/$ref.pairs" | sort -u > "$work/ref.taxon" ;;
    esac
    n=$(count "$work/$arm.$kind"); shared=$(comm -12 "$work/ref.$kind" "$work/$arm.$kind" | wc -l)
    lost=$(comm -23 "$work/ref.$kind" "$work/$arm.$kind" | wc -l); gained=$(comm -13 "$work/ref.$kind" "$work/$arm.$kind" | wc -l)
    printf '%s\t%s\t%s\t%s\t' "$n" "$shared" "$lost" "$gained"
  done > "$work/$arm.row"
  # The profiles' species (lineage column) and abundances; the "?" line left out. (Bray-Curtis computed before printf:
  # a ">" among its arguments would redirect its output.)
  species=$(awk -F'\t' 'FNR == 1 { f++ } NF >= 3 && $1 != "?" && $1 !~ /^#/ { a[f, $2] = $3 + 0; seen[$2] }
      END { for (s in seen) { in1 = ((1, s) in a); in2 = ((2, s) in a); x = in1 ? a[1, s] : 0; y = in2 ? a[2, s] : 0
              n2 += in2; lost += (in1 && !in2); gained += (in2 && !in1); num += (x < y ? x : y); s1 += x; s2 += y }
            bc = 0; if (s1 + s2 > 0) bc = 1 - 2 * num / (s1 + s2)
            printf "%d\t%d\t%d\t%.4f", n2, lost, gained, bc }' "$work/$ref.profile" "$work/$arm.profile")
  echo -e "$arm\t$(cat "$work/$arm.row")$species" >> "$work/compare.tsv"
done
column -t "$work/compare.tsv"
grep -h "^Sparse seeding" "$work"/*.1.log | head -1
