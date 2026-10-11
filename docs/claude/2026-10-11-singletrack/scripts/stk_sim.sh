#!/usr/bin/env bash
# stk_sim.sh: test samples for the Singletrack check, from the realworld_s1 genomes (6 GTDB r226 families, 487 genomes:
# representatives and strains) with the reference build's simulator: one paired-end sample (300k pairs, 40 species),
# one PacBio HiFi and one Nanopore sample (30 Mb each) of the same 40 genomes, weighted by length. On cores 0-3.
set -uo pipefail
W=$HOME/stk; D=$W/data; SIM=$W/ref/build/simulate_metagenomes
TABLE=$HOME/realworld_s1/work/genomes_simulated.tsv
ONT_MODEL=$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model
mkdir -p $D
run() { taskset -c 0-3 nice -n 5 "$@"; }
echo "== paired-end ($(date +%T))"
run $SIM --genome_table $TABLE -o $D/pe -n 1 --sample_prefix pe_s --total_read_pairs 300000 --species_per_sample 40 \
  --read_length 150 --seed 7 -t 4 --protal_metafile $D/pe/protal.meta > $D/pe.log 2>&1; echo "exit $?"
ls $D/pe/reads 2>/dev/null | head
M=$(ls $D/pe/manifest.tsv $D/pe/manifests/*.tsv 2>/dev/null | head -1); echo "manifest $M"; head -2 $M
# The long-read samples: the pe sample's genomes, weighted by relative abundance times length (as the build does).
mkdir -p $D/long
python3 - "$M" $D/long <<'EOF'
import csv, sys
manifest, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(manifest), delimiter="\t"))
with open(f"{out}/samples.tsv", "w") as fh:
    fh.write("sample\tout\tbases\tseed\n")
    fh.write(f"pb\t{out}/pb.fq.gz\t30000000\t101\n")
    fh.write(f"ont\t{out}/ont.fq.gz\t30000000\t202\n")
for name in ("pb", "ont"):
    with open(f"{out}/genomes_{name}.tsv", "w") as fh:
        fh.write("sample\tgenome\tfasta\tweight\thost\n")
        for r in rows:
            fh.write(f"{name}\t{r['genome']}\t{r['fasta_path']}\t{float(r['relative_abundance']) * float(r['genome_length'])!r}\t0\n")
print(len(rows), "genomes")
EOF
for name in pb ont; do
  echo "== $name ($(date +%T))"
  grep -P "^sample|^$name\t" $D/long/samples.tsv > $D/long/samples_$name.tsv
  if [ $name = pb ]; then setup="hifi:15000:3000:3"; model=""; else setup="qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36"; model="--long_model $ONT_MODEL"; fi
  run $SIM --long_samples $D/long/samples_$name.tsv --long_genomes $D/long/genomes_$name.tsv --long_setup "$setup" $model \
    --long_stats $D/long/stats_$name.tsv -t 4 > $D/long/$name.log 2>&1; echo "exit $?"
  ls -la $D/long/$name.fq.gz
done
echo "SIMDONE $(date +%T)"
