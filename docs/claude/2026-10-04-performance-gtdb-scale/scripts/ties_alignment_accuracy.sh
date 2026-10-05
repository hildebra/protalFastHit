#!/usr/bin/env bash
# For the kept SAMs (paired-end 500k pairs, full database): per run, the primary records (not 0x100/0x800) and the share
# whose gene's taxon is the species of the genome the read came from (the read name's accession; genome2tiid.tsv);
# and the counts line. Both versions side by side.
F=${TIES:-$HOME/tiebench}; G=$HOME/bench071/V075/protal_db/genome2tiid.tsv
for d in $F/runs_v075/v073.full.pe.*_p500000_*; do
  [ -d $d ] || continue
  id=${d##*.full.pe.}; for v in v073 v075; do
    r=$F/runs_v075/$v.full.pe.$id
    acc=$(zstd -dc $r/$id.sam.zst | awk -F'\t' -v G=$G 'BEGIN { while ((getline l < G) > 0) { split(l, a, "\t"); tax[a[1]] = a[2] } }
      /^@/ { next } and($2, 0x904) { next } { n++; split($3, g, "_"); match($1, /^GC[AF]_[0-9]+\.[0-9]+/); acc = substr($1, 1, RLENGTH)
        if (acc in tax) { known++; ok += (tax[acc] == g[1]) } }
      END { printf "%d primary records, %d of known genomes, %.5f on the true species", n, known, ok / (known ? known : 1) }')
    echo "$v $id: $acc; $(grep -h "candidate alignments tried" $r.log | sed 's/^Sample [^:]*: //' | cut -c1-230); $(grep -h 'Elapsed' $r.time | awk '{print $NF}') wall"
  done
done
