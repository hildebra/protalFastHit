#!/bin/bash
# Fix 3: the genome lengths build_gtdb_database.genome_length counts are the simulator's, and a table with them
# gives the same samples as one without. All genomes of the padded world are counted in Python and compared with
# the padding plan's sizes and with the simulator's own counts (simulate_metagenomes --test, every species in a
# sample); then one sample is simulated from each table with the same seed and the outputs compared.
# usage: check_lengths.sh SCRIPTS_DIR OUT
set -euo pipefail
scripts=$1 out=$2
sim=${SIM:-$HOME/build-prof/build/simulate_metagenomes}
rm -rf "$out" && mkdir -p "$out"
cut -f1-3 ~/bprof/world/genomes.tsv > "$out/three.tsv"
/usr/bin/time -o "$out/time.txt" -f "lengths of 2295 genomes on 4 threads: %e s wall, %U s user, %S s system" \
  python3 - "$scripts" "$out" > "$out/python.txt" <<'EOF'
import os, sys
sys.path.insert(0, sys.argv[1])
import build_gtdb_database as build
out = sys.argv[2]
four = build.with_lengths(os.path.join(out, "three.tsv"), os.path.join(out, "four.tsv"), 4)
rows = [l.rstrip("\n").split("\t") for l in open(four)]
lengths = [int(r[3]) for r in rows]
plan = {l.split("\t")[0]: int(l.split("\t")[2]) for l in open(os.path.expanduser("~/bprof/world/plan.tsv"))}
print("differ from the plan:", sum(plan[r[0]] != n for r, n in zip(rows, lengths)), "of", len(rows))
EOF
cat "$out/python.txt" "$out/time.txt"
# The simulator's own counts: every species once, three samples (--test: no reads).
"$sim" --genome_table "$out/three.tsv" -o "$out/test3" -n 3 --total_read_pairs 1000 --species_per_sample 765 --seed 2 --test > "$out/test3.log" 2>&1
awk -F'\t' 'NR == FNR { len[$1] = $4; next } FNR > 1 { n++; if (len[$2] != $5) bad++ } END { print n " genomes counted by the simulator, " bad + 0 " differ" }' \
  "$out/four.tsv" "$out/test3/manifest.tsv"
# The same sample from both tables.
for t in three four; do
  "$sim" --genome_table "$out/$t.tsv" -o "$out/sim_$t" -n 1 --total_read_pairs 20000 --species_per_sample 60 --seed 7 \
    --strains_per_species 0.3,0.1 -t 2 > "$out/sim_$t.log" 2>&1
done
for f in reads/sample_1_R1.fq.gz reads/sample_1_R2.fq.gz; do
  cmp -s "$out/sim_three/$f" "$out/sim_four/$f" && echo "$f identical" || echo "$f DIFFERS"
done
cmp -s <(sed 's/sim_three/SIM/g' "$out/sim_three/manifest.tsv") <(sed 's/sim_four/SIM/g' "$out/sim_four/manifest.tsv") \
  && echo "manifest.tsv identical but for its output folder" || echo "manifest.tsv DIFFERS"
grep -h "took" "$out"/sim_*.log
