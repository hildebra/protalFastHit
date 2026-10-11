#!/usr/bin/env bash
# The audit's long-read samples: the pe1m community (its genomes weighted by their read pairs) as PacBio HiFi and ONT
# reads (the build's default setups), 100 Mb each and 5 Mb for callgrind. 4 pinned cores.
set -uo pipefail
S=$HOME/perf9/samples
SIM=$HOME/perf9/v081/build/simulate_metagenomes
M=$S/pe1m/manifests/pe1m_1.tsv
L=$S/long; mkdir -p "$L"
pin=(taskset -c 0-3 nice -n 5)
make_tables() {  # kind bases...
    local kind=$1; shift
    printf 'sample\tout\tbases\tseed\n' > "$L/${kind}_samples.tsv"
    printf 'sample\tgenome\tfasta\tweight\thost\n' > "$L/${kind}_genomes.tsv"
    for b in "$@"; do
        local name=${kind}_${b}
        printf '%s\t%s\t%s\t%s\n' "$name" "$L/$name.fq.gz" "$b" 31 >> "$L/${kind}_samples.tsv"
        awk -F'\t' -v s="$name" 'NR > 1 { print s "\t" $2 "\t" $11 "\t" $6 "\t0" }' "$M" >> "$L/${kind}_genomes.tsv"
    done
}
make_tables pb 100000000 5000000
make_tables ont 100000000 5000000
"${pin[@]}" "$SIM" --long_samples "$L/pb_samples.tsv" --long_genomes "$L/pb_genomes.tsv" --long_setup hifi:15000:3000:3 \
    --long_stats "$L/pb_stats.tsv" -t 4 > "$L/pb.log" 2>&1 || { echo "FAIL pb"; tail -5 "$L/pb.log"; }
"${pin[@]}" "$SIM" --long_samples "$L/ont_samples.tsv" --long_genomes "$L/ont_genomes.tsv" \
    --long_setup qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36 --long_model "$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model" \
    --long_stats "$L/ont_stats.tsv" -t 4 > "$L/ont.log" 2>&1 || { echo "FAIL ont"; tail -5 "$L/ont.log"; }
cat "$L/pb_stats.tsv" "$L/ont_stats.tsv"
ls -la "$L"
# the callgrind pe subset: the first 100,000 pairs
for r in 1 2; do zcat "$S/pe1m/reads/pe1m_1_R$r.fq.gz" | head -400000 | gzip -1 > "$S/pe100k_R$r.fq.gz"; done
ls -la "$S"
