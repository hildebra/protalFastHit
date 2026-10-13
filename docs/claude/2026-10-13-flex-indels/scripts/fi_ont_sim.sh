#!/usr/bin/env bash
# ONT test samples for --flex_indels (docs/claude/2026-10-13-flex-indels), made as the build's collector makes a
# long-read point and as the sparse-first report's sf_hifi_sim.sh made HiFi ones: each sample replays the community of a
# paired-end test sample (its manifest's genomes, weighted by relative abundance times genome length), read by
# simulate_metagenomes --long_samples with pbsim3's qshmm Nanopore model (the GTDB build's --ont-setup,
# qshmm:QSHMM-ONT-HQ:8000:6000:ACCURACY:39/24/36), and its truth is that community's. Prints map rows (sample, reads,
# truth) to stdout; OUT/tmp/genomes.tsv lists each sample's genomes in the order the read names number them (g<k>x_<n>).
#   bash fi_ont_sim.sh SIMULATE PE_SIM GENOME_STORE MODEL OUT NAME BASES SAMPLES ACCURACY
set -uo pipefail
simulate=$1 pe_sim=$2 store=$3 model=$4 out=$5 name=$6 bases=$7 n=$8 accuracy=$9
mkdir -p "$out/tmp" "$out/reads"
pe_point=$(awk -F'\t' 'NR == 2 { s = $1; sub(/_s_[0-9]+$/, "", s); print s }' "$pe_sim/manifest.tsv")
echo -e "sample\tout\tbases\tseed" > "$out/tmp/samples.tsv"
echo -e "sample\tgenome\tfasta\tweight\thost" > "$out/tmp/genomes.tsv"
for s in $(seq 1 "$n"); do
  sample=${name}_s_$s community=${pe_point}_s_$s
  echo -e "$sample\t$out/reads/$sample.fq.zst\t$bases\t$(( (2000003 + s) * 101 ))" >> "$out/tmp/samples.tsv"
  awk -F'\t' -v c="$community" -v s="$sample" 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
       $1 == c { printf "%s\t%s\t%s\t%.17g\t0\n", s, $col["genome"], $col["fasta_path"], $col["relative_abundance"] * $col["genome_length"] }' \
      "$pe_sim/manifest.tsv" >> "$out/tmp/genomes.tsv"
done
if [ ! -s "$out/done" ]; then
  taskset -c 0-3 nice -n 5 "$simulate" --long_samples "$out/tmp/samples.tsv" --long_genomes "$out/tmp/genomes.tsv" \
      --long_setup "qshmm:QSHMM-ONT-HQ:8000:6000:$accuracy:39/24/36" --long_model "$model" --long_stats "$out/tmp/stats.tsv" \
      -t 4 --genome_store "$store" > "$out/simulate.log" 2>&1 \
      || { echo "simulation failed" >&2; tail "$out/simulate.log" >&2; exit 1; }
  touch "$out/done"
fi
for s in $(seq 1 "$n"); do
  echo -e "${name}_s_$s\t$out/reads/${name}_s_$s.fq.zst\t$pe_sim/protal_goldstd/${pe_point}_s_$s.profile_truth"
done
