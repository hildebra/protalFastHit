#!/bin/bash
# Run time, old against split: callgrind instruction counts of one run each (one thread, a map of a paired-end, a
# PacBio and an ONT sample from the e2e tests, profiled, no strains) on the mini database; the outputs compared too.
R=~/split-rt
rm -rf $R
mkdir -p $R/reads
E=~/split-cmp/tmp_base
cp $E/protal_e2e_reads__*/sa_R1.fq $E/protal_e2e_reads__*/sa_R2.fq $E/protal_e2e_read_types_*/la.fq $E/protal_e2e_read_types_*/oa.fq $R/reads/
for v in base split; do
    bin=~/split-base/build/protal; [ $v = split ] && bin=~/split-build/build/protal
    {
        printf '#OUTPUT_DIR\t%s\n' $R/out_$v
        printf '#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\n'
        printf 'pe\tpe\t%s\t%s\t-\n' $R/reads/sa_R1.fq $R/reads/sa_R2.fq
        printf 'pb\tpb\t%s\t-\tpb\n' $R/reads/la.fq
        printf 'ont\tont\t%s\t-\tont\n' $R/reads/oa.fq
    } > $R/map_$v.tsv
    /usr/bin/time -f "$v: %e s wall under callgrind" taskset -c 0-3 nice -n 5 valgrind --tool=callgrind --callgrind-out-file=$R/cg_$v.out \
        $bin --db $(ls -d ~/split-cmp/tmp_base/protal_e2e_read_type_db_*/db) --map $R/map_$v.tsv -t 1 --no_strains > $R/run_$v.log 2> $R/valgrind_$v.log
    tail -1 $R/valgrind_$v.log
done
for v in base split; do
    total=$(callgrind_annotate $R/cg_$v.out 2>/dev/null | grep -m1 'PROGRAM TOTALS' | awk '{print $1}')
    callgrind_annotate --inclusive=yes $R/cg_$v.out 2>/dev/null > $R/incl_$v.txt
    pick() { grep -m1 -E "$1" $R/incl_$v.txt | awk '{print $1}'; }
    echo "$v: total $total; RunWrapper $(pick 'protal::RunWrapper<protal::NoBenchmark>'); ProfileWrapper $(pick 'protal::ProfileWrapper\(protal::ProfilingContext'); LoadIndex $(pick 'protal::LoadIndex\('); RunPairedEnd $(pick 'classify::RunPairedEnd<'); RunLongReads $(pick 'classify::RunLongReads<')"
done
echo "== outputs (profiles; SAM records as sets)"
for s in pe pb ont; do
    for f in $s.profile $s.profile.log $s.profile.genes.log $s.profile.composition; do
        cmp -s $R/out_base/profiles/$f $R/out_split/profiles/$f 2>/dev/null && echo "$f identical" || { [ -e $R/out_base/$f ] && cmp -s $R/out_base/$f $R/out_split/$f && echo "$f identical" || echo "$f DIFFERS or missing"; }
    done
done
find $R/out_base -name '*.profile' | head -3
