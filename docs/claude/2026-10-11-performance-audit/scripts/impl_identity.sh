#!/usr/bin/env bash
# Same outputs: HEAD (ref) against the changes (work) on the small real world's database (~/realworld_s1: column weights
# and strain alleles as text, no binary forms), as it is (text) and with the binary forms --compress_db adds (bin);
# pe, se, HiFi and ONT samples simulated from its genomes; each sample alone, by default and with --sequential_load, and
# a cohort of the four by default and with --profile_ahead. Compares every output file but the timers (misc/) and the
# logs, and the SAMs' records. 4 pinned cores.
set -uo pipefail
W=$HOME/impl; I=$W/id; R=$W/ref/build/protal; N=$W/work/build/protal; SIM=$W/work/build/simulate_metagenomes
pin=(taskset -c 0-3 nice -n 5)
mkdir -p "$I"
# The databases.
if [ ! -s "$I/db_bin/database.protal" ]; then
  rm -rf "$I/db_text" "$I/db_bin"; mkdir -p "$I/db_text" "$I/db_bin"
  cp "$HOME/realworld_s1/protal_db/database.protal" "$I/db_text/"
  cp "$HOME/realworld_s1/protal_db/database.protal" "$I/db_bin/"
  "${pin[@]}" "$N" --compress_db --db "$I/db_bin/database.protal" -t 4 > "$I/compress.log" 2>&1 || { echo "compress FAIL"; tail "$I/compress.log"; exit 1; }
  grep -E "^Wrote|^Rewrite" "$I/compress.log"
fi
# The samples.
S=$I/samples
if [ ! -s "$S/long/ont.fq.gz" ]; then
  rm -rf "$S"; mkdir -p "$S/long"
  "${pin[@]}" "$SIM" --genome_table "$HOME/realworld_s1/work/genomes.tsv" -o "$S/pe" -n 1 --sample_prefix pe --total_read_pairs 200000 \
      --species_per_sample 40 --strains_per_species 0.4 --congener_groups 0.5:2-5 --read_length 150 --sequencer HSXt --seed 5 -t 4 > "$S/pe.log" 2>&1
  M=$S/pe/manifests/pe_1.tsv
  for kind in pb ont; do
    printf 'sample\tout\tbases\tseed\n%s\t%s\t%s\t%s\n' "$kind" "$S/long/$kind.fq.gz" 20000000 9 > "$S/long/${kind}_samples.tsv"
    { printf 'sample\tgenome\tfasta\tweight\thost\n'; awk -F'\t' -v s="$kind" 'NR > 1 { print s "\t" $2 "\t" $11 "\t" $6 "\t0" }' "$M"; } > "$S/long/${kind}_genomes.tsv"
  done
  "${pin[@]}" "$SIM" --long_samples "$S/long/pb_samples.tsv" --long_genomes "$S/long/pb_genomes.tsv" --long_setup hifi:15000:3000:3 -t 4 > "$S/pb.log" 2>&1
  "${pin[@]}" "$SIM" --long_samples "$S/long/ont_samples.tsv" --long_genomes "$S/long/ont_genomes.tsv" \
      --long_setup qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36 --long_model "$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model" -t 4 > "$S/ont.log" 2>&1
fi
pe1=$S/pe/reads/pe_1_R1.fq.gz pe2=$S/pe/reads/pe_1_R2.fq.gz
reads() {
  case $1 in
    pe) echo "-1 $pe1 -2 $pe2" ;; se) echo "-1 $pe1 --read_type se" ;;
    pb) echo "-1 $S/long/pb.fq.gz --read_type pb" ;; ont) echo "-1 $S/long/ont.fq.gz --read_type ont" ;;
  esac
}
run() {  # name binary db extra...
  local name=$1 bin=$2 db=$3; shift 3
  rm -rf "$I/runs/$name"; mkdir -p "$I/runs"
  # shellcheck disable=SC2068
  "${pin[@]}" "$bin" --db "$db/database.protal" -o "$I/runs/$name" -t 4 --no_qcmsa $@ > "$I/runs/$name.log" 2>&1 || echo "FAILED $name"
}
for s in pe se pb ont; do
  # shellcheck disable=SC2046
  run "ref_text_$s" "$R" "$I/db_text" $(reads $s)
  run "work_text_$s" "$N" "$I/db_text" $(reads $s)
  run "work_bin_$s" "$N" "$I/db_bin" $(reads $s)
  run "work_bin_seq_$s" "$N" "$I/db_bin" $(reads $s) --sequential_load
done
map() {  # out
  { printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\tREAD_TYPE\n' "$1"
    printf 'pe\t%s\t%s\tpe\tpe\nse\t%s\t-\tse\tse\npb\t%s\t-\tpb\tpb\nont\t%s\t-\tont\tont\n' "$pe1" "$pe2" "$pe1" "$S/long/pb.fq.gz" "$S/long/ont.fq.gz"; } > "$I/$(basename "$1").map"
}
for v in ref_text work_bin work_bin_ahead; do
  rm -rf "$I/runs/cohort_$v"; map "$I/runs/cohort_$v"
  case $v in ref_text) bin=$R db=$I/db_text opts= ;; work_bin) bin=$N db=$I/db_bin opts= ;; work_bin_ahead) bin=$N db=$I/db_bin opts=--profile_ahead ;; esac
  "${pin[@]}" "$bin" --db "$db/database.protal" --map "$I/cohort_$v.map" -t 4 --no_qcmsa $opts > "$I/runs/cohort_$v.log" 2>&1 || echo "FAILED cohort $v"
done
# Compare: every file but misc/ (timers) and the SAMs, whose records are compared sorted (the header names the run).
compare() {  # a b
  local a=$1 b=$2 bad=0
  while IFS= read -r f; do
    rel=${f#$a/}
    case $rel in misc/*|*.log) continue ;; esac
    if [[ $rel == *.sam.zst ]]; then
      if ! cmp -s <(zstd -dc "$a/$rel" | grep -v '^@' | sort) <(zstd -dc "$b/$rel" 2> /dev/null | grep -v '^@' | sort); then echo "  differs: $rel"; bad=1; fi
    elif ! cmp -s "$a/$rel" "$b/$rel"; then echo "  differs: $rel"; bad=1; fi
  done < <(find "$a" -type f | sort)
  echo "$([ $bad = 0 ] && echo same || echo DIFFERENT): $(basename "$a") vs $(basename "$b") ($(find "$a" -type f | grep -v '/misc/' | wc -l) files)"
}
for s in pe se pb ont; do
  for v in work_text work_bin work_bin_seq; do
    # the run folders are named after the run, so the outputs' names are the same in each (prefix from the reads)
    compare "$I/runs/ref_text_$s" "$I/runs/${v}_$s"
  done
done
compare "$I/runs/cohort_ref_text" "$I/runs/cohort_work_bin"
compare "$I/runs/cohort_ref_text" "$I/runs/cohort_work_bin_ahead"
grep -h -E "^(Column weights|Strain alleles):|^Loading the (taxonomy|profiling)" "$I/runs/work_bin_pe.log" "$I/runs/work_text_pe.log" | cut -c1-250
