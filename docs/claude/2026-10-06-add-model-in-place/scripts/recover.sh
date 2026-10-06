#!/bin/bash
# A stopped in-place --add_model through the command line: a run names the journal, --add_model writes the old
# bytes back and then replaces the models. Then the whole unit suite but the 13 GB PackedIndex test.
set -e
W=~/protal-addmodel/recover
NEW=~/protal-addmodel/build/protal
SCR=$(cd "$(dirname "$0")" && pwd)  # this folder (make_bigdb.cpp, make_journal.cpp)
M=~/protal-addmodel/bench/models
rm -rf $W && mkdir -p $W && cd $W
g++ -O2 -std=c++20 -I ~/protal-addmodel/src/src -o make_journal $SCR/make_journal.cpp -lzstd -lpthread
cp /home/falk/testaudit/src/data/mini_db/protal_db/database.protal db.protal
$NEW --add_model $M/trained_model.xml,$M/trained_model_se.xml,$M/trained_model_pb.xml,$M/trained_model_ont.xml \
     --read_type pe,se,pb,ont --db db.protal -t 4 > first.log 2>&1
grep -E "^Rewrite|^Replace" first.log
cp db.protal old.protal
$NEW --add_model $M/trained_model_ont.xml --read_type se --db db.protal -t 4 > in_place.log 2>&1
grep -E "^Rewrite|^Replace" in_place.log
./make_journal old.protal db.protal
truncate -s -20 db.protal   # stopped before the seek table was complete
echo "--- a run on the broken database:"
$NEW --unpack_db --db db.protal --unpack_dir u1 > run.log 2>&1 && echo "unexpectedly ran" || grep -iE "journal|cannot read" run.log | head -3
echo "--- --add_model:"
$NEW --add_model $M/trained_model_pb.xml --read_type ont --db db.protal -t 4 > recover.log 2>&1 || { tail -20 recover.log; exit 1; }
grep -E "^Wrote the old|^Replace|^Rewrite|Stored" recover.log
ls db.protal.journal 2>/dev/null && echo "journal left" || echo "no journal left"
$NEW --unpack_db --db db.protal --unpack_dir u2 > unpack.log 2>&1
cmp u2/model_se.xml $M/trained_model_se.xml && cmp u2/model_ONT.xml $M/trained_model_pb.xml && cmp u2/model_pe.xml $M/trained_model.xml \
  && echo "models: se as before the stopped run, ont the new one"
echo "--- a journal that is not of the file:"
cp old.protal other.protal
./make_journal old.protal other.protal
truncate -s -20 other.protal
d=$(od -An -t u8 -j 24 -N 8 other.protal.journal | tr -d ' ')   # the directory's size: the check bytes follow it
printf 'ZZZZZZZZ' | dd of=other.protal.journal bs=1 seek=$((48 + d)) conv=notrunc 2>/dev/null
$NEW --add_model $M/trained_model_pb.xml --read_type ont --db other.protal -t 4 > other.log 2>&1 && echo "unexpectedly ran" || grep -iE "journal" other.log | head -3

cd ~/protal-addmodel
./build/tests/protal_tests --gtest_filter='-PackedIndex.*' 2>&1 | tail -5
