#!/usr/bin/env bash
# Without --flex_indels nothing changes (docs/claude/2026-10-13-flex-indels): HEAD's protal (REF) and the changed one
# (NEW), without the option and with --flex_indels none, on a paired-end sample (with its database) and an ONT sample
# (with its); the SAM records (sorted, md5) and the profile compared. 4 threads on cores 0-3.
#   bash fi_exact.sh REF NEW WORK PE_DATABASE PE_R1 PE_R2 ONT_DATABASE ONT_READS
set -uo pipefail
ref=$1 new=$2 work=$3 pe_db=$4 r1=$5 r2=$6 ont_db=$7 ont=$8
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
records() { zstdcat "$1" | grep -v '^@' | sort | md5sum | cut -c1-16; }
for sample in pe ont; do
  reads=(-1 "$r1" -2 "$r2") db=$pe_db; [ $sample = ont ] && { reads=(-1 "$ont" --read_type ont); db=$ont_db; }
  for arm in ref new none; do
    bin=$new extra=(); [ $arm = ref ] && bin=$ref; [ $arm = none ] && extra=(--flex_indels none)
    out=$work/$sample.$arm; rm -rf "$out"
    "${pin[@]}" "$bin" --db "$db" "${reads[@]}" --prefix s -o "$out" -t 4 --no_strains "${extra[@]}" > "$out.log" 2>&1
    echo -e "$sample\t$arm\texit $?\trecords $(records "$out/s.sam.zst")\tprofile $(md5sum < "$out/s.profile" | cut -c1-16)"
  done
  for arm in new none; do
    cmp -s "$work/$sample.ref/s.profile" "$work/$sample.$arm/s.profile" && echo "$sample $arm: profile identical to HEAD's" ||
      echo "$sample $arm: profile DIFFERS from HEAD's"
  done
done
