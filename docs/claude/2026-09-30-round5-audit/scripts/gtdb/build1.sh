#!/bin/bash
# End-to-end build_gtdb_database.py --inputs on the downloaded synthetic world, 2 threads.
# usage: build1.sh OUTNAME [extra build_gtdb_database.py args]
OUT=~/audit6/gtdb/${1:-b1}; shift
PY=~/protal-train/bin/python
S=~/audit6/gtdb/src/scripts
BIN=~/audit6/gtdb/bin
mkdir -p $BIN
# protal wrapper: logs each call and its peak memory
cat > $BIN/protal <<EOF
#!/bin/bash
echo "\$(date +%T) start \$\$ \$*" >> $BIN/protal_calls.log
/usr/bin/time -f "%e s %M kB" -o $BIN/time.\$\$ ~/strain-build/bin/protal "\$@"
rc=\$?
echo "\$(date +%T) end \$\$ rc=\$rc \$(cat $BIN/time.\$\$) \$*" >> $BIN/protal_calls.log
exit \$rc
EOF
chmod +x $BIN/protal
echo "=== $(date +%T) build into $OUT $*" >> $BIN/protal_calls.log
cd ~/audit6/gtdb
/usr/bin/time -v $PY $S/build_gtdb_database.py --inputs ~/audit6/gtdb/dl/inputs --outdir $OUT \
    --protal $BIN/protal --simulator ~/strain-build/bin/simulate_metagenomes -t 2 \
    --samples 2 --read-pairs 1000,20000 --read-setups 100:HS20:300:40,150:HS25:350:50 \
    --species-per-sample 8-12 --archaea 1 --holdout-max-share 0.1 \
    --holdout-clades phylum:1,class:1,family:1,genus:2 "$@"
echo "exit $?"
