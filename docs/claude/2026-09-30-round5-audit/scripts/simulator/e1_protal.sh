#!/bin/bash
# E1b: profile the E1 simulation with protal (mini DB), then compare with the manifest truth.
set -u
W=~/audit6/simulator
cd $W/e1
rm -rf prot
/usr/bin/time -v ~/strain-build/bin/protal --db ~/strain-build/mini_db/protal_db/database.protal --map sim/protal.meta -t 2 > protal.log 2> protal.err
echo "exit $?"
grep -E "Maximum resident|Elapsed" protal.err
grep -E "TP|Sample" protal.log | head
ls -R prot | head -40
for s in s_1 s_2; do echo "== $s profile"; cat prot/profiles/$s.profile; echo "== truth (manifest, species sums)"; awk -F'\t' -v s=$s '$1==s{a[$3]+=$8; v[$3]+=$7} END{for(k in a) print k"\t"a[k]"\tvcov="v[k]}' sim/manifest.tsv; done
head -3 prot/profiles/s_1.profile.log | cut -c1-400
head -3 prot/profiles/s_1.profile.truth_annotated | cut -c1-300
