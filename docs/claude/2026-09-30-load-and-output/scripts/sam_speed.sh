#!/bin/bash
# Alignment with SAM output, 8 threads, --no_profile, map runs, alternated: OLD's .sam.gz (pigz
# afterwards) against NEW's .sam.gz, .sam.zst and .sam; then --profile_only (1 thread) on each file.
#   sam_speed.sh OLD_BIN NEW_BIN DB NAME R1 R2 [NAME R1 R2 ...]
source "$(dirname "$0")/env.sh"
old=$1; new=$2; db=$3; shift 3
O=$OUT/sam_speed; rm -rf $O; mkdir -p $O
secs() { awk -v s="$1" 'BEGIN { n = split(s, a, " "); t = 0; for (i = 1; i <= n; i++) { v = a[i]; if (v ~ /ms$/) t += substr(v, 1, length(v) - 2) / 1000; else if (v ~ /min$/) t += 60 * substr(v, 1, length(v) - 3); else if (v ~ /s$/) t += substr(v, 1, length(v) - 1) } printf "%.2f", t }'; }
names=()
while [ $# -ge 3 ]; do
  ds=$1; r1=$2; r2=$3; shift 3; names+=($ds)
  for rep in 1 2; do
    for v in old_gz new_gz new_zst new_sam; do
      b=$new; [ $v = old_gz ] && b=$old
      ext=${v#*_}; name=$ds.sam; [ $ext != sam ] && name=$ds.sam.$ext
      d=$O/${ds}_$v; rm -rf $d
      printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tSAM\n%s\t%s\t%s\t%s\t%s\n' $d $ds $ds $r1 $r2 $name > $O/$v.map
      /usr/bin/time -f "%e %U" $b --db $db --map $O/$v.map -t 8 --no_qcmsa --no_profile > $O/$v.log 2> $O/$v.time
      f=$(find $d -name "$name")
      align=$(secs "$(grep 'Aligning reads took' $O/$v.log | sed 's/Aligning reads took //')")
      echo -e "$ds rep$rep $v\twall $(cut -d' ' -f1 $O/$v.time) s\tuser $(cut -d' ' -f2 $O/$v.time) s\taligning $align s\t$(stat -c %s $f) bytes\tload $(cut -d' ' -f1 /proc/loadavg)"
    done
  done
done
echo "== --profile_only, 1 thread"
for ds in "${names[@]}"; do
  for rep in 1 2; do
    for v in old_gz new_gz new_zst new_sam; do
      ext=${v#*_}; name=$ds.sam; [ $ext != sam ] && name=$ds.sam.$ext
      f=$(find $O/${ds}_$v -name "$name")
      rm -rf $O/po; mkdir -p $O/po
      $new --db $db --profile_only $f --prefix $ds -o $O/po -t 1 --no_qcmsa > $O/po.log 2>&1
      echo -e "$ds rep$rep $v\t$(grep 'Profiling took' $O/po.log)"
    done
  done
done
