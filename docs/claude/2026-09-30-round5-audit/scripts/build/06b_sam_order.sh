#!/bin/bash
# Do the -t 2 SAMs differ only in record order? Compare sorted records and headers.
set -u
C=~/audit6/build/compare
for s in s1 s2 s3; do
  for r in baseline_t1 baseline_t2 avx2_t1 avx2_t2; do
    zstd -dc $C/out_$r/$s.sam.zst > $C/$r.$s.sam
  done
  echo "== $s: records per run: $(for r in baseline_t1 baseline_t2 avx2_t1 avx2_t2; do grep -vc '^@' $C/$r.$s.sam; done | tr '\n' ' ')"
  for r in baseline_t2 avx2_t1 avx2_t2; do
    a=$(grep -v '^@' $C/baseline_t1.$s.sam | sort | md5sum | cut -c1-12); b=$(grep -v '^@' $C/$r.$s.sam | sort | md5sum | cut -c1-12)
    h1=$(grep '^@' $C/baseline_t1.$s.sam | md5sum | cut -c1-12); h2=$(grep '^@' $C/$r.$s.sam | md5sum | cut -c1-12)
    echo "  baseline_t1 vs $r: sorted records $( [ $a = $b ] && echo equal || echo DIFFER ), header $( [ $h1 = $h2 ] && echo equal || echo DIFFER ), raw $(cmp -s $C/baseline_t1.$s.sam $C/$r.$s.sam && echo equal || echo differ)"
  done
done
rm -f $C/*.s?.sam
echo "== static build status"; cat ~/audit6/build/04_static.rc 2>&1; tail -3 ~/audit6/build/04_static.out
echo DONE
