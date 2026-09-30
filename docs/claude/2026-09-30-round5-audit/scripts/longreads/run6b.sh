#!/bin/bash
# Gene 48 of Mockella alpha: covered by paired-end reads, half (pb) or not at all (ont) by long reads.
S=$(dirname "$0")
W=~/audit6/longreads
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
OUT=$S/out_e6b.txt
cd $W/e6
{
m=$(awk '$2==48{print $1}' $D/gene2geneid.tsv); echo "gene 48 = $m"
grep -P "\t$m\t" $G/simulation/marker_positions.tsv
for acc in GCF_999001001.1 GCA_999001002.1 GCA_999001003.1; do python3 -c "
import sys; sys.path.insert(0,'$S'); from lrlib import read_fasta
g=read_fasta([l.split('\t')[2] for l in open('$G/simulation/genomes.tsv') if l.startswith('$acc')][0]); print('$acc', {k:len(v) for k,v in g.items()})"; done
for s in pe_a pb_a ont_a; do
  echo "-- $s records on *_48 (not secondary): $(grep -v '^@' pe_model/$s.sam | awk '!and($2,256) && $3 ~ /_48$/' | wc -l)"
  grep -v '^@' pe_model/$s.sam | awk '!and($2,256) && $3 ~ /_48$/ {h=$6; gsub(/[0-9]+[MXID]/,"",h); print $1,$2,$3,$4,$5,length($10),h}' | head -8
done
echo "== truth: long reads over gene 48 of GCA_999001002.1"
python3 - <<EOF
import pickle, sys
sys.path.insert(0, "$S")
from lrlib import World
w = World("$G", "$D")
acc = "GCA_999001002.1"
offs = pickle.load(open("pb_a.truth.pkl","rb"))["contig_offsets"][acc]
for mk in w.markers[acc]:
    if w.geneid[mk[0]] == 48:
        gs, ge = offs[mk[1]] + mk[2], offs[mk[1]] + mk[3]
        print("gene48 chrom", gs, ge, mk)
for f in ("pb_a", "ont_a"):
    T = pickle.load(open(f + ".truth.pkl", "rb"))
    n = 0
    for t in T["reads"]:
        if t["start"] <= gs and t["start"] + t["length"] >= ge:
            n += 1
            if n <= 3: print(f, t["name"], t["start"], t["length"], t["strand"])
    print(f, "reads spanning gene 48 wholly:", n)
EOF
} > $OUT 2>&1
echo done
