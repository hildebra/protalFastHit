#!/bin/bash
set -eu
B=/mnt/c/Users/hildebra/Documents/locDev/protal/scripts/mock_benchmark
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/a049aa2e-d238-4e0a-83b9-f87b6acf1fe4/scratchpad/samples_plus.tsv
PY=$HOME/micromamba/envs/protal-db-build/bin/python
taskset -c 0-3 $PY $B/download_mocks.py -o ~/Mocks --samples $S --sets d6331_pe_labs --max-pairs 1e6 --jobs 2 2>&1 | grep -v "^ENA" | tail -3
taskset -c 0-3 nice -n 5 $PY $B/score_mocks.py run --reads ~/Mocks --samples $S --sets d6331_pe_labs \
    --protal ~/mocktest/protal/build/protal --db ~/mocktest/build/protal_db -t 4 2>&1 | grep -v "^  /"
$PY $B/score_mocks.py score --reads ~/Mocks --samples $S --runs 0.8.3_4dfc549 > /dev/null
for s in d6331_dna_SRR38158090 d6331_hku_SRR23926922 d6331_auckland_SRR33710519; do
  f=$(find ~/Mocks/0.8.3_4dfc549 -name "$s.profile.log")
  echo "== $s"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $c["Lineage"] ~ /s__(Veillonella rogosae|Bacteroides fragilis|Escherichia coli|Prevotella corporis|Roseburia hominis|Faecalibacterium prausnitzii_E)$/ {n=split($c["Lineage"],a,";"); match($c["Summary"], /MeanANI: [0-9.]+/); printf "%-36s P=%.3f VCov=%7.3f ANI=%s\n", a[n], $c["Probability"], $c["VCov"], substr($c["Summary"],RSTART+9,6)}' $f
done
awk -F'\t' 'NR==1 || $1=="d6331_pe_labs" || $1=="d6331_dna"' ~/Mocks/0.8.3_4dfc549/scores/samples.tsv | cut -f2,10-19 | column -t -s $'\t'
