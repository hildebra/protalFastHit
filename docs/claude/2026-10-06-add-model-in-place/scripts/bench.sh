#!/bin/bash
# Times protal --add_model on a ~20 GB database.protal: HEAD (15006c2, full rewrite) against the in-place
# replacement, and checks that HEAD reads a database the new binary changed in place.
set -e
W=~/protal-addmodel/bench
NEW=~/protal-addmodel/build/protal
OLD=~/protal-addmodel/head/build/protal
MINI=/home/falk/testaudit/src/data/mini_db/protal_db/database.protal
SCR=$(cd "$(dirname "$0")" && pwd)  # this folder (make_bigdb.cpp, make_journal.cpp)
PAD_GB=${PAD_GB:-20}
rm -rf $W && mkdir -p $W/models $W/logs && cd $W
echo "start $(date); $(uptime)"
for f in trained_model trained_model_se trained_model_pb trained_model_ont; do
  cp /mnt/c/Users/hildebra/Documents/locDev/protal/local/v5_adjacency/$f.xml models/
done
ls -l models
M=models/trained_model.xml,models/trained_model_se.xml,models/trained_model_pb.xml,models/trained_model_ont.xml
T=pe,se,pb,ont

g++ -O2 -std=c++20 -I ~/protal-addmodel/src/src -o make_bigdb $SCR/make_bigdb.cpp -lzstd -lpthread
/usr/bin/time -f "make_bigdb: %e s wall, %M KB" ./make_bigdb $MINI big.protal $PAD_GB 6 > logs/make_bigdb.log 2>&1
cat logs/make_bigdb.log
ls -l big.protal

run() {  # label binary args...
  local label=$1 bin=$2; shift 2
  /usr/bin/time -f "$label: %e s wall, %U s user, %S s sys, %M KB max RSS" -o logs/$label.time \
    $bin "$@" > logs/$label.log 2>&1 || { echo "$label failed"; tail -20 logs/$label.log; exit 1; }
  cat logs/$label.time
  grep -E "^Replace|^Rewrite|database.protal: |took|Stored" logs/$label.log | grep -v "^Stored" | head -5
  ls -l big.protal | awk '{print "  size " $5}'
}
run old_rewrite  $OLD --add_model $M --read_type $T --db big.protal -t 6
run new_first    $NEW --add_model $M --read_type $T --db big.protal -t 6
run new_in_place $NEW --add_model $M --read_type $T --db big.protal -t 6
run new_in_place_again $NEW --add_model $M --read_type $T --db big.protal -t 6
run new_one_model $NEW --add_model models/trained_model_se.xml --read_type se --db big.protal -t 6
run new_one_model_1t $NEW --add_model models/trained_model_se.xml --read_type se --db big.protal -t 1
run new_in_place_1t $NEW --add_model $M --read_type $T --db big.protal -t 1
rm -f big.protal

# HEAD reads what the new binary wrote in place (the mini database).
cp $MINI mini.protal
$NEW --add_model $M --read_type $T --db mini.protal -t 6 > logs/mini_first.log 2>&1
$NEW --add_model models/trained_model_pb.xml,models/trained_model_ont.xml --read_type ont,pb --db mini.protal -t 6 > logs/mini_in_place.log 2>&1
grep -E "^Replace|^Rewrite" logs/mini_first.log logs/mini_in_place.log
$OLD --unpack_db --db mini.protal --unpack_dir unpacked_old > logs/unpack_old.log 2>&1
$NEW --unpack_db --db mini.protal --unpack_dir unpacked_new > logs/unpack_new.log 2>&1
diff -r unpacked_old unpacked_new && echo "HEAD and the new binary unpack the same files"
cmp unpacked_old/model_PB.xml models/trained_model_ont.xml && cmp unpacked_old/model_ONT.xml models/trained_model_pb.xml && \
  cmp unpacked_old/model_se.xml models/trained_model_se.xml && cmp unpacked_old/model_pe.xml models/trained_model.xml && echo "models as given"
$OLD --add_model models/trained_model.xml --read_type pe --db mini.protal -t 6 > logs/mini_old_rewrite.log 2>&1 && echo "HEAD rewrites it"
echo "end $(date); $(uptime)"
