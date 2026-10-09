#!/usr/bin/env bash
# The protal-run timings of a build_gtdb_database.py share archive, as used in the report (README.md):
# per protal run, the alignment, profiling and start-up stages; per read type, the alignment time; per design point,
# samples, reads, alignment seconds and reads per second.
#
#   run_timings.sh LOGS_DIR        (the archive's logs/, e.g. local/v19/protal0.7.9_r226_v19/logs)
#
# Only awk and grep: it reads protal_runs_training.log and protal_runs_test.log (protal.log of every run, concatenated
# with "==> runN <==" lines by build_gtdb_database.py's report()).
set -euo pipefail
logs=${1:?usage: run_timings.sh LOGS_DIR}

# protal's Benchmark durations ("1h 2m 3s 456ms") in seconds.
secs='function secs(s,   t,n,i,a){ t=0; n=split(s,a," "); for(i=1;i<=n;i++){ if(a[i]~/h$/)t+=a[i]*3600; else if(a[i]~/ms$/)t+=a[i]/1000; else if(a[i]~/m$/)t+=a[i]*60; else if(a[i]~/s$/)t+=a[i]+0;} return t}'

for f in protal_runs_training.log protal_runs_test.log; do
    echo "== $f: per protal run (seconds)"
    awk "$secs"'
    function flush() { if (run != "") printf "%-6s samples %4d  align %6.0f  processing %6.0f  profiling %5.0f  index %4.0f  preload %4.0f  tables %4.0f  run %6.0f\n", run, na, al, proc, prof, ld, pre, tab, rp }
    /^==> run/ { flush(); run=$2; al=na=proc=prof=ld=pre=tab=rp=0 }
    /^Aligning reads took/ { sub(/^Aligning reads took /,""); al+=secs($0); na++ }
    /^Processing all samples took/ { sub(/^Processing all samples took /,""); proc=secs($0) }
    /^Profiling took/ { sub(/^Profiling took /,""); prof=secs($0) }
    /^Load Index took/ { sub(/^Load Index took /,""); ld=secs($0) }
    /^Preload genomes took/ { sub(/^Preload genomes took /,""); pre=secs($0) }
    /^Loading the taxonomy, models and tables took/ { sub(/^Loading the taxonomy, models and tables took /,""); tab=secs($0) }
    /^Run protal took/ { sub(/^Run protal took /,""); rp=secs($0) }
    END { flush() }' "$logs/$f"
done

echo "== alignment by read type, both collections"
cat "$logs/protal_runs_training.log" "$logs/protal_runs_test.log" | awk "$secs"'
/^Align the / { ty=$3 }
/^Aligning reads took/ { sub(/^Aligning reads took /,""); t=secs($0); al[ty]+=t; tot+=t; n[ty]++ }
END { for (t in al) printf "%-11s %4d samples %6.0f s (%.1f%%)\n", t, n[t], al[t], 100*al[t]/tot; printf "total %.0f s\n", tot }'

echo "== per design point (both collections): samples, reads (pairs for paired-end), alignment s, reads/s"
cat "$logs/protal_runs_training.log" "$logs/protal_runs_test.log" | awk "$secs"'
/^Align the / { name=$0; sub(/.*of sample /,"",name); sub(/ \(-a.*/,"",name); dp=name; sub(/_s_[0-9]+(_se)?$/,"",dp); if (name ~ /_se$/) dp=dp "_se" }
/^Sample .*: [0-9]+ read/ { match($0,/: [0-9]+ read/); reads[dp]+=substr($0,RSTART+2,RLENGTH-7)+0 }
/^Aligning reads took/ { sub(/^Aligning reads took /,""); al[dp]+=secs($0); cnt[dp]++ }
END { for (d in al) printf "%-34s %4d %11d %8.1f %9.0f\n", d, cnt[d], reads[d], al[d], reads[d]/(al[d]+1e-9) }' | sort -k4 -g -r

echo "== host against community samples (the same read type, depth and length)"
for s in sc_host_pe_p10000000_s_1 sc_moderate_pe_p10000000_s_1 sc_host_ont_b3000000000_s_1 sc_moderate_ont_b3000000000_s_1 sc_host_pb_b3000000000_s_1; do
    grep -h -A14 "of sample $s (" "$logs/protal_runs_training.log" | grep -E "^Sample $s: |^Aligning reads took" | head -2
done
