#!/usr/bin/env bash
# The speed of KeepSparseLookups: base (5334f2c) and reuse (+ the patch) alternated, REPS rounds, at --sparse_first 1
# (and plain, no gate, from the base binary, for scale), aligning only (--no_profile), 4 threads on cores 0-3. Per run:
# wall and CPU seconds (/usr/bin/time), the "Aligning reads took" timer and, from misc/<sample>_runtime.tsv, the
# "Seed- and Anchor-finding", "Seeding" and "Sparse first" stages (seconds summed over the threads). Writes
# WORK/speed.tsv and the medians.
#   bash sfr_speed.sh BASE_PROTAL REUSE_PROTAL DATABASE WORK REPS -- READ ARGUMENTS (e.g. -1 R1 -2 R2)
set -uo pipefail
base=$1 reuse=$2 db=$3 work=$4 reps=$5; shift 5
[ "${1:-}" = "--" ] && shift
reads=("$@")
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
read -r -a order <<< "${ARMS:-base_c1 reuse_c1 plain}"
out=$work/speed.tsv
echo -e "arm\trep\twall_s\tuser_s\tsys_s\taligning_s\tseed_anchor_s\tseeding_s\tsparse_s\trecords_written" > "$out"
seconds() {  # "1m 2s 345ms" or "2s 345ms" or "345ms" -> seconds
  awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) s += substr(v, 1, length(v) - 2) / 1000;
         else if (v ~ /m$/) s += 60 * substr(v, 1, length(v) - 1); else if (v ~ /s$/) s += substr(v, 1, length(v) - 1) } printf "%.3f", s }'
}
for rep in $(seq 1 "$reps"); do
  for arm in "${order[@]}"; do
    case $arm in
      base_c1) bin=$base; args=(--sparse_first 1) ;;
      reuse_c1) bin=$reuse; args=(--sparse_first 1) ;;
      plain) bin=$base; args=() ;;
    esac
    dir=$work/run_$arm; rm -rf "$dir"
    log=$work/$arm.$rep.log
    "${pin[@]}" /usr/bin/time -f "TIME %e %U %S" "$bin" --db "$db" "${reads[@]}" --prefix s -o "$dir" -t 4 \
        --no_strains --no_profile "${args[@]}" > "$log" 2>&1
    read -r wall user sys < <(grep "^TIME" "$log" | awk '{ print $2, $3, $4 }')
    aligning=$(grep -m1 "Aligning reads took" "$log" | sed 's/.*took //' | seconds)
    runtime=$(ls "$dir"/misc/*_runtime.tsv 2> /dev/null | head -1)
    stage() { awk -F'\t' -v s="$1" '$1 == s { printf "%.3f", $2 }' "$runtime"; }
    seed_anchor=$(stage "Seed- and Anchor-finding"); seeding=$(stage "Seeding"); sparse=$(stage "Sparse first")
    written=$(grep -m1 "records written" "$log" | sed -E 's/.* ([0-9]+) records written.*/\1/')
    echo -e "$arm\t$rep\t$wall\t$user\t$sys\t$aligning\t$seed_anchor\t$seeding\t${sparse:-0}\t$written" | tee -a "$out"
  done
done
echo "== medians"
for arm in "${order[@]}"; do
  awk -F'\t' -v a="$arm" '$1 == a { w[++n] = $3; u[n] = $4; al[n] = $6; sa[n] = $7; se[n] = $8; sp[n] = $9; r = $10 }
       END { asort(w); asort(u); asort(al); asort(sa); asort(se); asort(sp); m = int((n + 1) / 2)
             printf "%s\twall %.2f\tuser %.2f\taligning %.3f\tseed+anchor %.3f\tseeding %.3f\tsparse %.3f\trecords %s\n",
                    a, w[m], u[m], al[m], sa[m], se[m], sp[m], r }' "$out"
done
