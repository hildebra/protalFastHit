#!/bin/bash
# Were archaea simulated (--archaea 1)? Which archaeal species are held out?
O=~/audit6/gtdb/${1:-b2}
echo "archaeal species in the genome table:"; grep -o 'd__Archaea;[^	]*' $O/genomes.tsv | awk -F';' '{print $NF}' | sort -u
echo "held out archaeal species:"; grep -F -f <(grep -o 'd__Archaea;[^	]*' $O/genomes.tsv | awk -F';' '{print $NF}' | sort -u) $O/heldout_species.txt
echo "archaeal species in the whole taxonomy:"; awk -F'\t' '$5=="species"' $O/internal_taxonomy.dmp | wc -l
for p in $O/training/points/*; do
  echo "== $(basename $p)"
  grep -o -- '--taxon [^ ]*' $p/sim/run_params.tsv 2>/dev/null || grep -i taxon $p/sim/run_params.tsv | head -3
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} {split($h["taxonomy"], t, ";"); print $h["sample"], t[1], t[7]}' $p/sim/manifest.tsv | sort -u | awk '{print $1, $2}' | sort | uniq -c
  grep -i -E "archaea|demand|taxon" $p/simulate.log | head -5
done
