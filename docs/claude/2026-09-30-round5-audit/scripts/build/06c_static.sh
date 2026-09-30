#!/bin/bash
# The static binary on the 06_compare samples (-t 1) vs the dynamic baseline; RUNPATH of the
# dynamic binaries (build tree and installed copies).
set -u
A=~/audit6/build
C=$A/compare
S=$A/src
R1=$C/reads/s1_R1.fq,$C/reads/s2_R1.fq,$C/reads/s3_R1.fq
R2=$C/reads/s1_R2.fq,$C/reads/s2_R2.fq,$C/reads/s3_R2.fq
o=$C/out_static_t1; rm -rf $o
( cd $C && taskset -c 0,1 $S/build/protal_0.6.0_static --db $C/db -1 $R1 -2 $R2 --prefix s1,s2,s3 -o $o -t 1 --qcmsa_script $S/scripts/qcmsa.py > $C/log_static_t1.txt 2>&1 ); echo "static t1 rc=$?; $(tail -1 $C/log_static_t1.txt)"
( cd $o && find . -type f | sort | while read f; do case "$f" in *.zst) echo "$(zstd -dc "$f" | md5sum | cut -c1-32)  ${f%.zst} (zst, content)";; *) echo "$(md5sum < "$f" | cut -c1-32)  $f";; esac; done ) > $C/md5_static_t1.txt
echo "static_t1 vs baseline_t1: $(diff $C/md5_baseline_t1.txt $C/md5_static_t1.txt | grep -c '^[<>]') differing lines"; diff $C/md5_baseline_t1.txt $C/md5_static_t1.txt | grep '^[<>]' | grep -v runtime.tsv
echo "== RUNPATH/RPATH of dynamic binaries"
for b in $S/build/protal $S/build/protal_avx2 $A/prefix/bin/protal_baseline $A/prefix/bin/protal_avx2; do echo "-- $b"; readelf -d $b | grep -E 'RPATH|RUNPATH' | sed 's/:/\n    /g' | head -16; done
echo "== NEEDED of the installed baseline"; readelf -d $A/prefix/bin/protal_baseline | grep NEEDED
echo "== who pulls dlopen into the static binary"; nm $S/build/protal_0.6.0_static | grep -E ' (T|U|W) (dlopen|__libc_dlopen_mode|_dl_open)$' | head; grep -rln 'dlopen' $S/src 2>/dev/null | head -3
echo DONE
