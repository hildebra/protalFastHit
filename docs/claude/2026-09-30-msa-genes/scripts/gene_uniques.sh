#!/bin/bash
S=$(cd "$(dirname "$0")" && pwd)
head -2 $HOME/genes_study/world/reference.map; head -2 $HOME/genes_study/world/unique_kmers.tsv
python3 $S/gene_uniques.py $HOME/genes_study/world $HOME/genes_study/tuneH2 $HOME/genes_study/db900
