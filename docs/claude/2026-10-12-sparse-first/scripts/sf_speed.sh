#!/usr/bin/env bash
# The speed of --sparse_first on a synthetic metagenome (docs/claude/2026-10-12-sparse-first): the arms alternated,
# REPS rounds, aligning only (--no_profile), 4 threads on cores 0-3. An arm is "plain" or cC_tT (--sparse_first C
# --sparse_first_score T); ARMS (environment) lists them, default "plain c1_t14 c2_t14 c4_t14". Per run: wall and CPU
# seconds (/usr/bin/time), the "Aligning reads took" timer, the records written, the k-mers the kept reads looked up,
# the reads given up and the sampled lookups. Writes WORK/speed.tsv and the medians.
#   bash sf_speed.sh PROTAL DATABASE WORK REPS -- READ ARGUMENTS (e.g. -1 R1 -2 R2, or -1 HIFI --read_type pb)
set -uo pipefail
protal=$1 db=$2 work=$3 reps=$4; shift 4
[ "${1:-}" = "--" ] && shift
reads=("$@")
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
read -r -a order <<< "${ARMS:-plain c1_t14 c2_t14 c4_t14}"
arm_args() {  # cC_tT -> --sparse_first C --sparse_first_score T
  [ "$1" = plain ] && return
  local c=${1#c}; c=${c%%_*}; local t=${1##*_t}
  echo "--sparse_first $c --sparse_first_score $t"
}
out=$work/speed.tsv
echo -e "arm\trep\twall_s\tuser_s\tsys_s\taligning_s\trecords_written\tkept_kmers\tgiven_up\tsampled_lookups" > "$out"
seconds() {  # "1m 2s 345ms" or "2s 345ms" or "345ms" -> seconds
  awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) s += substr(v, 1, length(v) - 2) / 1000;
         else if (v ~ /m$/) s += 60 * substr(v, 1, length(v) - 1); else if (v ~ /s$/) s += substr(v, 1, length(v) - 1) } printf "%.3f", s }'
}
for rep in $(seq 1 "$reps"); do
  for arm in "${order[@]}"; do
    dir=$work/run_$arm; rm -rf "$dir"
    log=$work/$arm.$rep.log
    # shellcheck disable=SC2046
    "${pin[@]}" /usr/bin/time -f "TIME %e %U %S" "$protal" --db "$db" "${reads[@]}" --prefix s -o "$dir" -t 4 \
        --no_strains --no_profile $(arm_args "$arm") > "$log" 2>&1
    read -r wall user sys < <(grep "^TIME" "$log" | awk '{ print $2, $3, $4 }')
    aligning=$(grep -m1 "Aligning reads took" "$log" | sed 's/.*took //' | seconds)
    written=$(grep -m1 "records written" "$log" | sed -E 's/.* ([0-9]+) records written.*/\1/')
    kmers=$(grep -m1 " seeding: " "$log" | sed -E 's/.*seeding: ([0-9]+) k-mers.*/\1/')
    given=$(grep -m1 "^Sample .* sparse first" "$log" | sed -E 's/.*: ([0-9]+) of [0-9]+ (pairs|reads) given up.*/\1/')
    sampled=$(grep -m1 "^Sample .* sparse first" "$log" | sed -E 's/.* ([0-9]+) sampled lookups.*/\1/')
    echo -e "$arm\t$rep\t$wall\t$user\t$sys\t$aligning\t$written\t$kmers\t${given:-0}\t${sampled:-0}" | tee -a "$out"
  done
done
echo "== medians"
for arm in "${order[@]}"; do
  awk -F'\t' -v a="$arm" '$1 == a { w[++n] = $3; u[n] = $4; al[n] = $6; r = $7; k = $8; g = $9; s = $10 }
       END { asort(w); asort(u); asort(al); m = int((n + 1) / 2)
             printf "%s\twall %.2f\tuser %.2f\taligning %.2f\trecords %s\tkept k-mers %s\tgiven up %s\tsampled %s\n", a, w[m], u[m], al[m], r, k, g, s }' "$out"
done
