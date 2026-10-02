#!/usr/bin/env bash
# Flank sizes after re-seeding (PROTAL_LR_DEBUG2), Nanopore 3 Mb: per alignment the read and gene bases left of
# the first link and right of the last, and the free ends.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/perf3/build_exp.sh lr || exit 1
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
for v in reseed noreseed; do env="PROTAL_LR_DEBUG2=1"; [ $v = noreseed ] && env="$env PROTAL_LR_NORESEED=1"
  rm -rf $W/f.$v
  env $env $HOME/mt-work/perf3-lr/src/build/protal --db $DB -1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont --no_profile --prefix s -o $W/f.$v -t 1 --no_qcmsa 2> $W/f.$v.err > /dev/null
  grep '^LRFLK' $W/f.$v.err > $W/f.$v.tsv
  echo "== $v: $(wc -l < $W/f.$v.tsv) anchored alignments"
  awk -F'\t' '{ n++; links += $2; rl += $4; gl += $5; rr += $8; gr += $9; win += $13; fr += $7; fe += $11
       lg[n] = $5; rg[n] = $9 }
    END { printf "links %.1f; left flank: read %.0f gene %.0f (free read %.0f); right flank: read %.0f gene %.0f (free read %.0f); window gene %.0f (means)\n", links/n, rl/n, gl/n, fr/n, rr/n, gr/n, fe/n, win/n }' $W/f.$v.tsv
  awk -F'\t' '{ print $5; print $9 }' $W/f.$v.tsv | sort -n | awk '{ a[NR] = $1 } END { printf "gene bases in a flank: median %d, p75 %d, p90 %d, p99 %d\n", a[int(NR/2)], a[int(NR*0.75)], a[int(NR*0.9)], a[int(NR*0.99)] }'
done
