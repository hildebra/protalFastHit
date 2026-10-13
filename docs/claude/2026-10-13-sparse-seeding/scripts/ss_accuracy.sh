#!/usr/bin/env bash
# The calls of HiFi test samples against their truth with and without --sparse_seeding
# (docs/claude/2026-10-13-sparse-seeding): the small real world build's two HiFi points made again as the sparse-first
# report did (its sf_hifi_sim.sh: each sample replays a regenerated paired-end community; pb_b3000000 from
# rl150_p10000, pb_b30000000 from rl150_p100000, 12 samples each), profiled with the build's database and models
# (trained on runs without a gate), one run of all 24 per arm, 4 threads on cores 0-3. Arms as ss_compare.sh's
# (plain, cC, cCs); ARMS (environment), default "plain c8 c4s c8s c16s". Writes WORK/calls.tsv (TP, FP, FN and records
# per sample), WORK/f1.tsv (by point and over all) and WORK/arms.tsv (wall and CPU seconds).
#   bash ss_accuracy.sh PROTAL SIMULATE BUILD WORK
set -uo pipefail
protal=$1 simulate=$2 build=$3 work=$4
here=$(cd "$(dirname "$0")" && pwd)
hifi_sim=$here/../../2026-10-12-sparse-first/scripts/sf_hifi_sim.sh
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
db=$build/protal_db/database.protal
read -r -a order <<< "${ARMS:-plain c8 c4s c8s c16s}"
arm_args() {
  case $1 in
    plain) echo "--sparse_first 0" ;;
    c*s) local c=${1#c}; echo "--sparse_first ${c%s} --sparse_seeding" ;;
    c*) echo "--sparse_first ${1#c}" ;;
  esac
}

# 1. The paired-end points the HiFi points replay: the recorded command with this simulator and new output folders.
for point in rl150_p10000 rl150_p100000; do
  dir=$work/sim/$point
  [ -s "$dir/done" ] && continue
  rm -rf "$dir"; mkdir -p "$dir"
  mapfile -t cmd < <(sed -n '/"command": \[/,/\]/p' "$build/work/test/points/$point/simulated.json" | sed '1d;$d' |
                     sed -E 's/^ *"//; s/",?$//')
  cmd[0]=$simulate
  for i in "${!cmd[@]}"; do
    [ "${cmd[$i]}" = "-o" ] && cmd[$((i + 1))]=$dir/sim
    [ "${cmd[$i]}" = "--protal_metafile" ] && cmd[$((i + 1))]=$dir/protal
  done
  "${pin[@]}" "${cmd[@]}" > "$dir/simulate.log" 2>&1 || { echo "simulation of $point failed"; tail "$dir/simulate.log"; exit 1; }
  touch "$dir/done"
done

# 2. The HiFi samples and their map rows (sample, reads, truth).
rows=$work/hifi_rows.tsv
{
  bash "$hifi_sim" "$simulate" "$work/sim/rl150_p10000/sim" "$build/work/genome_store" "$work/hifi/pb_b3000000" pb_b3000000 3000000 12
  bash "$hifi_sim" "$simulate" "$work/sim/rl150_p100000/sim" "$build/work/genome_store" "$work/hifi/pb_b30000000" pb_b30000000 30000000 12
} > "$rows" || { echo "HiFi simulation failed"; exit 1; }
echo "$(wc -l < "$rows") HiFi samples"

# 3. The arms: one run of all samples each (a map with the arm's output folder).
echo -e "arm\twall_s\tuser_s" > "$work/arms.tsv"
echo -e "arm\tsample\tTP\tFP\tFN\trecords\tgiven_up" > "$work/calls.tsv"
for arm in "${order[@]}"; do
  out=$work/arms/$arm; rm -rf "$out"; mkdir -p "$out"
  map=$work/$arm.map
  { echo -e "#OUTPUT_DIR\t$out"; echo -e "#SAMPLEID\tPREFIX\tFIRST\tSECOND\tPROFILE_TRUTH\tREAD_TYPE"
    awk -F'\t' '{ printf "%s\t%s\t%s\t-\t%s\tpb\n", $1, $1, $2, $3 }' "$rows"; } > "$map"
  # shellcheck disable=SC2046
  "${pin[@]}" /usr/bin/time -f "TIME %e %U %S" "$protal" --db "$db" --map "$map" -t 4 --no_strains $(arm_args "$arm") \
      > "$work/$arm.log" 2>&1
  echo "$arm exit $?"
  grep "^TIME" "$work/$arm.log" | awk -v a="$arm" '{ printf "%s\t%s\t%s\n", a, $2, $3 }' >> "$work/arms.tsv"
  awk -v a="$arm" '
    /^Sample [^ ]+: [0-9]+ (read pairs|reads),/ { s = $2; sub(/:$/, "", s); n = split($0, w, " "); for (i = 1; i < n; i++) if (w[i + 1] == "records" && w[i + 2] == "written") rec[s] = w[i] }
    /^Sample [^ ]+ sparse first/ { s = $2; for (i = 3; i < NF; i++) if ($(i + 1) == "of") { g[s] = $i; break } }
    /^Sample [^ ]+: TP [0-9]+, FP [0-9]+, FN [0-9]+/ { s = $2; sub(/:$/, "", s); tp = $4; fp = $6; fn = $8; gsub(/,/, "", tp); gsub(/,/, "", fp); gsub(/[,(]/, "", fn); row[s] = tp "\t" fp "\t" fn }
    END { for (s in row) printf "%s\t%s\t%s\t%s\t%s\n", a, s, row[s], rec[s], (s in g) ? g[s] : 0 }' "$work/$arm.log" >> "$work/calls.tsv"
done

# 4. F1 by arm and point, and over all samples.
awk -F'\t' 'NR > 1 { d = $2; sub(/_s_[0-9]+$/, "", d); tp[$1, d] += $3; fp[$1, d] += $4; fn[$1, d] += $5; rec[$1, d] += $6
              tp[$1, "all"] += $3; fp[$1, "all"] += $4; fn[$1, "all"] += $5; rec[$1, "all"] += $6; arms[$1]; ds[d]; ds["all"] }
     END { printf "arm\tpoint\tTP\tFP\tFN\tprecision\trecall\tF1\trecords\n"
           for (a in arms) for (d in ds) { p = tp[a, d] / (tp[a, d] + fp[a, d]); r = tp[a, d] / (tp[a, d] + fn[a, d])
             printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\n", a, d, tp[a, d], fp[a, d], fn[a, d], p, r, 2 * p * r / (p + r), rec[a, d] } }' \
    "$work/calls.tsv" | sort -k2,2 -k1,1 | tee "$work/f1.tsv"
cat "$work/arms.tsv"
