#!/usr/bin/env bash
# After impl_final.sh (the final build): (1) a --no_profile run's SAM against a profiled run's, ref and work (ref's
# --no_profile left the column weights unread, so the weighted site shift had none); (2) the loaders' instructions of the
# final build (binary headers read without their whole frame), profile-only under callgrind on one thread, beside the
# ASan job (cores 2 and 3; instruction counts do not mind).
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
until grep -q "FINAL DONE (round 2)" "$A/impl_round2.out" 2> /dev/null; do sleep 15; done
W=$HOME/impl/perf; S=$W/samples; R=$W/post
REF=$HOME/impl/ref/build/protal; WORK=$HOME/impl/work/build/protal
rm -rf "$R"; mkdir -p "$R"
pin=(taskset -c 0-3 nice -n 5)
for e in 1 2; do zcat "$S/pe1m/reads/pe1m_1_R${e}.fq.gz" | head -400000 | gzip -1 > "$R/pe100k_R$e.fq.gz"; done
echo "[$(date +%T)] SAMs"
run() {  # name bin db opts...
  local name=$1 bin=$2 db=$3; shift 3
  "${pin[@]}" "$bin" --db "$db/database.protal" -1 "$R/pe100k_R1.fq.gz" -2 "$R/pe100k_R2.fq.gz" -o "$R/$name" -t 4 --no_qcmsa "$@" > "$R/$name.log" 2>&1 \
      || echo "FAILED $name"
}
run ref_full "$REF" "$W/db_ref"
run ref_noprofile "$REF" "$W/db_ref" --no_profile
run work_full "$WORK" "$W/db_bin"
run work_noprofile "$WORK" "$W/db_bin" --no_profile
records() { zstd -dc "$(find "$R/$1" -name "*.sam.zst" | head -1)" | grep -v '^@' | sort; }
for pair in "ref_full work_full" "ref_full work_noprofile" "ref_full ref_noprofile" "work_full work_noprofile"; do
  set -- $pair
  if cmp -s <(records "$1") <(records "$2"); then echo "same SAM records: $1 vs $2"; else echo "different SAM records: $1 vs $2 ($(diff <(records "$1") <(records "$2") | grep -c '^[<>]') lines differ)"; fi
done
grep -h "^Column weights" "$R"/*.log | cut -c1-200

echo "[$(date +%T)] callgrind, the final build"
T=$W/tables; sam=$(find "$R/ref_full" -name "*.sam.zst" | head -1)
core=2
for v in work_text work_bin; do
  case $v in work_text) db=$W/db_ref ;; work_bin) db=$W/db_bin ;; esac
  ( taskset -c $core nice -n 5 valgrind --tool=callgrind --callgrind-out-file="$R/callgrind.out.$v" \
        "$WORK" --db "$db/database.protal" --profile_only "$sam" -o "$R/cg_$v" -t 1 --no_qcmsa --sequential_load > "$R/cg.$v.log" 2>&1
    callgrind_annotate --inclusive=yes --threshold=100 "$R/callgrind.out.$v" 2> /dev/null > "$R/$v.all.txt" ) &
  core=$((core + 1))
done
wait
names=("PROGRAM TOTALS" "protal::LoadStrainAlleles(" "LoadColumnWeights(" "binary_table::ReadHeader(" "binary_table::ReadArrays(" "ReadFramesPrefix("
       "ZSTD_decompressDCtx" "LoadAncestrySites(" "LoadCongenerGaps(")
printf 'function\twork_text\twork_bin\n' > "$R/callgrind.tsv"
for n in "${names[@]}"; do
  printf '%s\t%s\t%s\n' "$n" "$(grep -F -m1 "$n" "$R/work_text.all.txt" | awk '{ print $1 }')" "$(grep -F -m1 "$n" "$R/work_bin.all.txt" | awk '{ print $1 }')" >> "$R/callgrind.tsv"
done
cat "$R/callgrind.tsv"
rm -rf "$R"/cg_* "$R"/pe100k_R*.fq.gz
echo "[$(date +%T)] POST DONE"
