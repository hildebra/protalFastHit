#!/usr/bin/env bash
# oneb_nocpu.sh's gdb run without clearing the feature bits: on this CPU the .arch_x86_64_v3 clones should run.
W=$HOME/mt-work/onebcpu
grep -v '^set var' $W/gdb.cmd | sed 's|out\.nocpu|out.v3|g; s|log\.nocpu|log.v3|g' > $W/gdb_v3.cmd
grep '^run' $W/gdb_v3.cmd | cut -c1-200
gdb -batch -x $W/gdb_v3.cmd ${CL:-$HOME/mt-work/oneb/src/build/protal} > $W/gdb_v3.log 2>&1
grep -E 'before:|exited normally|already hit|^[0-9]+ +breakpoint' $W/gdb_v3.log | sed -E 's/\(.*\[clone/ [clone/' | cut -c1-150
cmp -s <(zstdcat $W/out.base/s.sam.zst) <(zstdcat $W/out.v3/s.sam.zst) && echo "SAM text same as plain x86-64" || echo "SAM DIFFER"
diff -r -q -x '*_runtime.tsv' -x '*.sam.zst' $W/out.base $W/out.v3 && echo "outputs same"
