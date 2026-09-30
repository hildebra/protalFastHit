#!/bin/bash
# The installed launcher on real work: (1) `just example` (examples/mini_db/run.sh) with
# PROTAL=<prefix>/bin/protal, (2) a 3-sample run through the launcher with qcmsa found next to the
# binary, (3) which qcmsa wins when PROTAL_QCMSA_SCRIPT is also set (docs/installation.md:98 says
# next-to-binary and $PATH come first), (4) the same run with the launcher forced to baseline.
set -u
A=~/audit6/build
P=$A/prefix
C=$A/compare
S=$A/src
cd $S
echo "== (1) examples/mini_db/run.sh through the installed launcher"
rm -rf $A/example_work
{ time PROTAL=$P/bin/protal THREADS=2 taskset -c 0,1 bash examples/mini_db/run.sh $A/example_work ; } > $A/example.log 2>&1; echo "rc=$?"
tail -12 $A/example.log
grep -E '^protal' $A/example_work/run_info.txt
R1=$C/reads/s1_R1.fq,$C/reads/s2_R1.fq,$C/reads/s3_R1.fq
R2=$C/reads/s1_R2.fq,$C/reads/s2_R2.fq,$C/reads/s3_R2.fq
echo "== (2) 3 samples through the launcher, qcmsa next to the binary"
rm -rf $A/inst_run; mkdir -p $A/inst_run
( cd $A/inst_run && taskset -c 0,1 $P/bin/protal --db $C/db -1 $R1 -2 $R2 --prefix s1,s2,s3 -o out -t 2 > log.txt 2>&1 ); echo "rc=$?"
grep -iE 'qcmsa' $A/inst_run/log.txt | head -5; ls $A/inst_run/out/strains 2>/dev/null | head; find $A/inst_run/out -name '*.msa.fna' | head
echo "== (3) PROTAL_QCMSA_SCRIPT set as well: which qcmsa runs?"
mkdir -p $A/fakeq; printf '#!/bin/sh\necho "fake qcmsa from PROTAL_QCMSA_SCRIPT ran" >> %s/fakeq/marker\nexit 1\n' $A > $A/fakeq/qcmsa_env; chmod +x $A/fakeq/qcmsa_env
rm -f $A/fakeq/marker; rm -rf $A/inst_run2; mkdir -p $A/inst_run2
( cd $A/inst_run2 && PROTAL_QCMSA_SCRIPT=$A/fakeq/qcmsa_env taskset -c 0,1 $P/bin/protal --db $C/db -1 $R1 -2 $R2 --prefix s1,s2,s3 -o out -t 2 > log.txt 2>&1 ); echo "rc=$?"
cat $A/fakeq/marker 2>/dev/null || echo "fake qcmsa did not run"
echo "== (4) forced baseline through the launcher: same outputs as via AVX2?"
rm -rf $A/inst_run3; mkdir -p $A/inst_run3
( cd $A/inst_run3 && PROTAL_NO_AVX2=1 taskset -c 0,1 $P/bin/protal --db $C/db -1 $R1 -2 $R2 --prefix s1,s2,s3 -o out -t 2 > log.txt 2>&1 ); echo "rc=$?"
h() { ( cd $1 && find . -type f | sort | while read f; do case "$f" in *.zst) echo "$(zstd -dc "$f" | md5sum | cut -c1-32)  ${f%.zst}";; *) echo "$(md5sum < "$f" | cut -c1-32)  $f";; esac; done ); }
h $A/inst_run/out > $A/inst_run.md5; h $A/inst_run3/out > $A/inst_run3.md5
echo "files: $(wc -l < $A/inst_run.md5); differing: $(diff $A/inst_run.md5 $A/inst_run3.md5 | grep -c '^[<>]')"; diff $A/inst_run.md5 $A/inst_run3.md5 | head
echo DONE
