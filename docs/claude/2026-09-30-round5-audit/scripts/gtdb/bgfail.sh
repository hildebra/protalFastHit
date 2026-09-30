#!/bin/bash
# E-c: the background build of the finished database fails early (as an out-of-memory kill would);
# when does build_gtdb_database.py notice?
OUT=~/audit6/gtdb/b4
rm -rf $OUT
PY=~/protal-train/bin/python
S=~/audit6/gtdb/src/scripts
BIN=~/audit6/gtdb/binfail
mkdir -p $BIN
cat > $BIN/protal <<EOF
#!/bin/bash
echo "\$(date +%T) start \$*" >> $BIN/calls.log
case " \$* " in
  *" --build "*"--db $OUT/protal_db "*) sleep 3; echo "\$(date +%T) fake failure of the final build" >> $BIN/calls.log; exit 137;;
esac
~/strain-build/bin/protal "\$@"
rc=\$?
echo "\$(date +%T) end rc=\$rc \$*" >> $BIN/calls.log
exit \$rc
EOF
chmod +x $BIN/protal
: > $BIN/calls.log
START=$(date +%s)
$PY $S/build_gtdb_database.py --inputs ~/audit6/gtdb/dl/inputs --outdir $OUT \
    --protal $BIN/protal --simulator ~/strain-build/bin/simulate_metagenomes -t 2 \
    --samples 2 --read-pairs 1000,20000 --read-setups 100:HS20:300:40,150:HS25:350:50 \
    --species-per-sample 8-12 --archaea 1 --holdout-max-share 0.1 \
    --holdout-clades phylum:1,class:1,family:1,genus:2 --evaluation basic
echo "exit $? after $(( $(date +%s) - START )) s"
cat $BIN/calls.log | cut -c1-160
ls $OUT
