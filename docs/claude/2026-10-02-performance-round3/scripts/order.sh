#!/usr/bin/env bash
W=$HOME/mt-work/perf3
cmp -s <(zstdcat $W/pe5.ref/s.sam.zst) <(zstdcat $W/pe5.ref2/s.sam.zst) && echo "reference twice: SAM in the same order" || echo "reference twice: SAM records in another order"
cmp -s <(zstdcat $W/pe5.ref/s.sam.zst) <(zstdcat $W/pe5.lr/s.sam.zst) && echo "reference vs prototype: SAM in the same order" || echo "reference vs prototype: SAM records in another order"
# profile the prototype's SAM with the reference binary: same gene log as the prototype's own run?
rm -rf $W/pe5.prof; nice $W/ref/build/protal --db $HOME/bench071/V071/protal_db --profile_only $W/pe5.lr/s.sam.zst --prefix s -o $W/pe5.prof -t 6 --no_qcmsa > /dev/null 2>&1
cmp -s $W/pe5.prof/s.profile.gene.log $W/pe5.lr/s.profile.gene.log && echo "prototype's SAM profiled by the reference: same gene log as the prototype" || echo "prototype's SAM profiled by the reference: gene log differs"
