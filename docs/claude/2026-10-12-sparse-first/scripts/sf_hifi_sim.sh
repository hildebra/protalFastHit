#!/usr/bin/env bash
# HiFi test samples for the --sparse_first accuracy check (docs/claude/2026-10-12-sparse-first), made as the build's
# collector makes a long-read point (collect_training_data.py long_unit_plan): each sample replays the community of a
# paired-end sample (its manifest's genomes, weighted by relative abundance times genome length), read by
# simulate_metagenomes --long_samples with the build's HiFi setup (hifi:15000:3000:3), and its truth is that
# community's truth file. Prints map rows (sample, reads, truth) to stdout.
#   bash sf_hifi_sim.sh SIMULATE PE_SIM GENOME_STORE OUT NAME BASES SAMPLES
# PE_SIM: a regenerated paired-end point's sim folder (manifest.tsv, protal_goldstd/); GENOME_STORE: the build's
# (work/genome_store), or - to read the manifest's FASTA files; NAME: the point's name (pb_b3000000); BASES: read bases
# per sample; SAMPLES: the first SAMPLES communities (<pe point>_s_1 ...), or those COMMUNITIES (environment) lists.
set -uo pipefail
simulate=$1 pe_sim=$2 store=$3 out=$4 name=$5 bases=$6 n=$7
mkdir -p "$out/tmp" "$out/reads"
pe_point=$(awk -F'\t' 'NR == 2 { s = $1; sub(/_s_[0-9]+$/, "", s); print s }' "$pe_sim/manifest.tsv")
read -r -a communities <<< "${COMMUNITIES:-}"
echo -e "sample\tout\tbases\tseed" > "$out/tmp/samples.tsv"
echo -e "sample\tgenome\tfasta\tweight\thost" > "$out/tmp/genomes.tsv"
for s in $(seq 1 "$n"); do
  sample=${name}_s_$s community=${communities[$((s - 1))]:-${pe_point}_s_$s}
  echo -e "$sample\t$out/reads/$sample.fq.zst\t$bases\t$(( (1000003 + s) * 101 ))" >> "$out/tmp/samples.tsv"
  awk -F'\t' -v c="$community" -v s="$sample" 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
       $1 == c { printf "%s\t%s\t%s\t%.17g\t0\n", s, $col["genome"], $col["fasta_path"], $col["relative_abundance"] * $col["genome_length"] }' \
      "$pe_sim/manifest.tsv" >> "$out/tmp/genomes.tsv"
done
store_args=(); [ "$store" != "-" ] && store_args=(--genome_store "$store")
if [ ! -s "$out/done" ]; then
  taskset -c 0-3 nice -n 5 "$simulate" --long_samples "$out/tmp/samples.tsv" --long_genomes "$out/tmp/genomes.tsv" \
      --long_setup hifi:15000:3000:3 --long_stats "$out/tmp/stats.tsv" -t 4 "${store_args[@]}" > "$out/simulate.log" 2>&1 \
      || { echo "simulation failed" >&2; tail "$out/simulate.log" >&2; exit 1; }
  touch "$out/done"
fi
for s in $(seq 1 "$n"); do
  echo -e "${name}_s_$s\t$out/reads/${name}_s_$s.fq.zst\t$pe_sim/protal_goldstd/${pe_point}_s_$s.profile_truth"
done
