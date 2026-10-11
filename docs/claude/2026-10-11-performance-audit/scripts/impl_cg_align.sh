#!/usr/bin/env bash
# Exact: the alignment's instructions, ref against work (the weighted site shift reads the packed column weights now),
# after impl_tables.sh. The first 100,000 pairs of pe1m, aligned only (--no_profile), one thread, under callgrind; the
# two variants side by side on cores 1-2.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
until grep -q "TABLES DONE" "$A/impl_tables.out" 2> /dev/null; do sleep 15; done
W=$HOME/impl/perf; S=$W/samples; R=$W/cg_align
REF=$HOME/impl/ref/build/protal; WORK=$HOME/impl/work/build/protal
rm -rf "$R"; mkdir -p "$R"
for e in 1 2; do zcat "$S/pe1m/reads/pe1m_1_R${e}.fq.gz" | head -400000 | gzip -1 > "$R/pe100k_R$e.fq.gz"; done
echo "[$(date +%T)] callgrind"
core=1
for v in ref_text work_bin; do
  case $v in ref_text) bin=$REF db=$W/db_ref ;; work_bin) bin=$WORK db=$W/db_bin ;; esac
  ( taskset -c $core nice -n 5 valgrind --tool=callgrind --callgrind-out-file="$R/callgrind.out.$v" \
        "$bin" --db "$db/database.protal" -1 "$R/pe100k_R1.fq.gz" -2 "$R/pe100k_R2.fq.gz" -o "$R/out_$v" -t 1 --no_profile > "$R/$v.log" 2>&1
    callgrind_annotate --inclusive=yes --threshold=100 "$R/callgrind.out.$v" 2> /dev/null > "$R/$v.all.txt" ) &
  core=$((core + 1))
done
wait
names=("PROGRAM TOTALS" "protal::LoadIndex(" "LoadAllGenomes(" "protal::LoadStrainAlleles(" "LoadColumnWeights(" "classify::RunPairedEnd<"
       "SimpleAlignmentHandler::operator()" "SimpleAlignmentHandler::ScoreAlleles" "ColumnWeightsOf" "column_weights::Table::Expand"
       "column_weights::Table::Find")
printf 'function\tref_text\twork_bin\n' > "$R/callgrind.tsv"
for n in "${names[@]}"; do
  printf '%s\t%s\t%s\n' "$n" "$(grep -F -m1 "$n" "$R/ref_text.all.txt" | awk '{ print $1 }')" "$(grep -F -m1 "$n" "$R/work_bin.all.txt" | awk '{ print $1 }')" >> "$R/callgrind.tsv"
done
cat "$R/callgrind.tsv"
# the SAM records the same
a=$(find "$R/out_ref_text" -name "*.sam.zst" | head -1); b=$(find "$R/out_work_bin" -name "*.sam.zst" | head -1)
if cmp -s <(zstd -dc "$a" | grep -v '^@' | sort) <(zstd -dc "$b" | grep -v '^@' | sort); then echo "SAM records same"; else echo "SAM RECORDS DIFFER"; fi
rm -rf "$R"/out_* "$R"/pe100k_R*.fq.gz
echo "[$(date +%T)] CG ALIGN DONE"
