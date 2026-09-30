cd ~/audit5/smoke
cut -f1-3,7 sim/manifest.tsv | sort -k3,3 -k1,1 | column -t
for f in protal/strains/*.raw.msa.fna; do
  echo "$f: $(grep -c '>' $f) seqs"
  grep '>' $f | tr '\n' ' '; echo
  awk 'NR==2{print length($0)}' $f
done
