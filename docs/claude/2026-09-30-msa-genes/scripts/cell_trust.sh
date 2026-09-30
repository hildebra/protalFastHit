#!/bin/bash
S=$(cd "$(dirname "$0")" && pwd)
python3 $S/cell_trust.py $HOME/genes_study/cong/out/strains $HOME/genes_study/tuneH2 $HOME/genes_study/cong/manifest.tsv
