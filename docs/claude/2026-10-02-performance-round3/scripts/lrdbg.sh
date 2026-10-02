#!/usr/bin/env bash
# The prototype with the debug print: per aligned candidate its gene length, window, links, the gene span of
# its links and the bases they cover; summaries per read type.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/perf3/build_exp.sh lrdbg || exit 1
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
for s in ont_b3000000 pb_b3000000; do t=ont; [ $s = pb_b3000000 ] && t=pb
  rm -rf $W/d.$s
  PROTAL_LR_DEBUG=1 $HOME/mt-work/perf3-lrdbg/src/build/protal --db $DB -1 $P/$s/sim/reads/${s}_s_1.fq.gz --read_type $t --no_profile --prefix s -o $W/d.$s -t 1 --no_qcmsa 2> $W/d.$s.err > /dev/null
  grep '^LRDBG' $W/d.$s.err > $W/d.$s.tsv
  echo "== $s: $(wc -l < $W/d.$s.tsv) aligned candidates"
  awk -F'\t' '{ n++; gl = $5; links += $7; left = $8; right = gl - $9; if (right < 0) right = 0; L += left; R += right; cov += $10 / gl; span += ($9 - $8) / gl
       k = $2 "|" $3 "|" $4; seen[k]++ }
    END { dup = 0; for (k in seen) if (seen[k] > 1) dup += seen[k] - 1
          printf "links per candidate %.1f; gene before the first link %.0f bases, after the last %.0f (means); links span %.0f%% of the gene, cover %.0f%%; same read+gene+strand aligned again: %d of %d\n", links / n, L / n, R / n, 100 * span / n, 100 * cov / n, dup, n }' $W/d.$s.tsv
  awk -F'\t' '{ left = $8; right = $5 - $9; if (right < 0) right = 0; print (left > right ? left : right) }' $W/d.$s.tsv | sort -n | awk '{ a[NR] = $1 } END { printf "longer flank (gene part): median %d, p75 %d, p90 %d\n", a[int(NR/2)], a[int(NR*0.75)], a[int(NR*0.9)] }'
done
# all anchors (candidates) of a read: per read+gene+strand, how many and the union of their links' gene spans
for s in ont_b3000000 pb_b3000000; do
  grep '^LRANC' $W/d.$s.err > $W/d.$s.anc.tsv
  echo "== $s anchors: $(wc -l < $W/d.$s.anc.tsv)"
  awk -F'\t' '{ k = $2 "|" $3 "|" $4; n[k]++; gl[k] = $5; if (!(k in lo) || $7 < lo[k]) lo[k] = $7; if ($8 > hi[k]) hi[k] = $8; sp = ($8 - $7) / $5; if (sp > best[k]) best[k] = sp }
    END { for (k in n) { g++; a += n[k]; multi += n[k] > 1; u += (hi[k] - lo[k]) / gl[k]; b += best[k] }
          printf "read+gene+strand: %d, anchors each %.2f, with more than one %.0f%%; span of the longest anchor %.0f%%, of all anchors together %.0f%% of the gene\n", g, a / g, 100 * multi / g, 100 * b / g, 100 * u / g }' $W/d.$s.anc.tsv
done
