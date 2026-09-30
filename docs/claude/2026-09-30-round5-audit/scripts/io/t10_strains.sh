#!/bin/bash
# T10: two samples with strain MSAs: identical outputs for -t 1 plain SAM vs -t 2 zst vs -t 2 gz,
# and after --profile_only of the zst SAMs; rerun with another --sam_format reuses an existing plain SAM.
set -u
source $(dirname "$0")/lib.sh
T=$W/t10; rm -rf $T; mkdir -p $T; cd $T
args="--db $DB -1 $R/sa_R1.fq,$R/sb_R1.fq -2 $R/sa_R2.fq,$R/sb_R2.fq --prefix sa,sb --no_qcmsa --msa_min_hcov 0"
prun a.log $args -o a -t 1 --sam_format sam
prun b.log $args -o b -t 2 --sam_format zst
prun c.log $args -o c -t 2 --sam_format gz
prun d.log --db $DB --profile_only b/sa.sam.zst,b/sb.sam.zst --prefix sa,sb -o d -t 2 --no_qcmsa --msa_min_hcov 0
sig() { (cd $1 && find . -type f ! -name '*.sam*' ! -name '*runtime*' ! -name '*.log' | sort | xargs md5sum); }
for d in a b c d; do sig $d > $d.sig; echo "$d: $(wc -l < $d.sig) files, strains files: $(find $d/strains -type f | wc -l)"; done
for d in b c d; do echo "a vs $d: $(diff a.sig $d.sig | grep -c '^[<>]') differing lines"; diff a.sig $d.sig | head -6; done
echo "--- rerun of a (plain SAMs) with the default --sam_format zst"
prun a2.log $args -o a -t 2
grep -a -E 'Skip|All alignments' a2.log | head -3
ls a/*.sam*
sig a > a2.sig; echo "a rerun vs a: $(diff a.sig a2.sig | grep -c '^[<>]') differing lines"
