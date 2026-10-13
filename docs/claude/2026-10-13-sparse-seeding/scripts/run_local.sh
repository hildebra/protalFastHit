#!/usr/bin/env bash
# The local runs of docs/claude/2026-10-13-sparse-seeding (WSL, cores 0-3):
#   1. --sparse_seeding changes nothing where the gate looks up every core (paired-end, C 1);
#   2. 300 Mb of HiFi reads of the archaea world's community (~/seed), arms compared read by read (ss_compare.sh);
#   3. the small real world's 24 HiFi test samples against their truth (ss_accuracy.sh).
#   bash run_local.sh PROTAL SIMULATE WORK
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
protal=$1 simulate=$2 work=$3
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
db=~/seed/db/database.protal
records() { zstdcat "$1" | grep -v '^@' | sort | md5sum | cut -c1-16; }

echo "== 1. paired-end C 1 with and without --sparse_seeding $(date +%T)"
for arm in c1 c1s; do
  extra=(); [ $arm = c1s ] && extra=(--sparse_seeding)
  rm -rf "$work/exact_$arm"
  "${pin[@]}" "$protal" --db "$db" -1 ~/seed/samples/pe100k_R1.fq.gz -2 ~/seed/samples/pe100k_R2.fq.gz --prefix s \
      -o "$work/exact_$arm" -t 4 --no_strains --sparse_first 1 "${extra[@]}" > "$work/exact_$arm.log" 2>&1
  echo "$arm exit $?"
done
a=$(records "$work/exact_c1/s.sam.zst"); b=$(records "$work/exact_c1s/s.sam.zst")
echo "SAM records: $([ "$a" = "$b" ] && echo same || echo DIFFER)"
cmp -s "$work/exact_c1/s.profile" "$work/exact_c1s/s.profile" && echo "profile: same" || echo "profile: DIFFERS"
diff <(grep -E "^Sample s( seeding: | sparse first \(|: .*records written)" "$work/exact_c1.log") \
     <(grep -E "^Sample s( seeding: | sparse first \(|: .*records written)" "$work/exact_c1s.log") > /dev/null && echo "count lines: same" || echo "count lines: DIFFER"
grep -h "^Sparse seeding" "$work/exact_c1s.log"

echo "== 2. archaea world, 300 Mb of HiFi $(date +%T)"
hifi=$work/hifi300
COMMUNITIES=pe1m_1 bash "$here/../../2026-10-12-sparse-first/scripts/sf_hifi_sim.sh" "$simulate" ~/seed/samples/pe1m - "$hifi" pb300m \
    300000000 1 > /dev/null || { echo "HiFi simulation failed"; exit 1; }
bash "$here/ss_compare.sh" "$protal" "$db" "$work/compare_archaea" 3 -- -1 "$hifi/reads/pb300m_s_1.fq.zst" --read_type pb

echo "== 3. small real world, 24 HiFi samples against their truth $(date +%T)"
bash "$here/ss_accuracy.sh" "$protal" "$simulate" ~/realworld_sp_s1 "$work/accuracy"
echo "== done $(date +%T)"
