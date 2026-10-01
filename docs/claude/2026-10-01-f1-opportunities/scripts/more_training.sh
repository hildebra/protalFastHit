#!/usr/bin/env bash
# A second training collection for 0.7.1's pipeline, as the pipeline collects its own (the same design, database,
# held-out species and protal) but with seed 101: twice the training data, to see whether F1 still grows with it.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
V=$B/V071
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
$PY $B/src/0.7.1/scripts/collect_training_data.py --db $V/training_db --genome_table $V/genomes.tsv -o $B/V071_training2 \
  --protal $B/bin/protal-0.7.1 --simulator $B/bin/simulate_metagenomes --samples 4 --read_pairs 1000,20000,200000 \
  --read_setups 100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50 --archaea 2 --species_per_sample 20-200 --seed 101 -t 6 \
  --taxonomy $V/internal_taxonomy.dmp --congeners 0 --read_types pe,se,pb,ont --long_read_bases 300000,6000000,60000000 \
  --pb_setup errhmm:ERRHMM-SEQUEL:15000:3000:0.999 --ont_setup qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36 \
  --pbsim $HOME/micromamba/envs/protal-db-build/bin/pbsim --strains_per_species 0.3,0.1 \
  --novel_species $V/heldout_species.txt --novel_clades 1 > $B/logs/training2.log 2>&1 || { tail -20 $B/logs/training2.log; exit 1; }
grep -E "^[0-9]+ taxa in" $B/logs/training2.log
