#!/bin/bash
# Experiment 2: reads around and above the 65 kb chunk size, slid over a reference genome so that
# genes fall on chunk edges and core boundaries; both strands. HiFi and ONT errors.
set -u
S=$(dirname "$0")
W=~/audit6/longreads
P=${P:-~/strain-build/bin/protal}
G=~/strain-build/mini_db/gtdb_r226
D=$W/mini_db
ACC=GCF_999002001.1
cd $W
mkdir -p e2 && cd e2
gen() {  # name platform length offsets
  [ -f $1.fq ] || python3 $S/gen_reads.py sliding $1 --gtdb $G --db $D --acc $ACC --platform $2 --length $3 --offsets $4 --seed $RANDOM --name "$1_{}"
}
gen h65000 hifi 65000 0:204979:2503
gen h65001 hifi 65001 0:204978:2501
gen h65535 hifi 65535 0:204444:2497
gen h65536 hifi 65536 0:204443:2493
gen h70000 hifi 70000 0:199979:2503
gen h130000 hifi 130000 0:139979:1511
gen h200000 hifi 200000 0:69979:1999
gen o130000 ont 130000 0:139979:1511
gen o200000 ont 200000 0:69979:1999
HIFI="h65000 h65001 h65535 h65536 h70000 h130000 h200000"
ONT="o130000 o200000"
files() { local out=""; for s in "$@"; do out="$out${out:+,}$s.fq"; done; echo $out; }
prefixes() { local out=""; for s in "$@"; do out="$out${out:+,}$s"; done; echo $out; }
$P --db $D -1 $(files $HIFI) --prefix $(prefixes $HIFI) --read_type pb --model_pb $D/model_pe.xml -o out -t 2 --sam_format sam --no_qcmsa --no_strains > hifi.log 2>&1; echo "hifi rc=$?"
$P --db $D -1 $(files $ONT) --prefix $(prefixes $ONT) --read_type ont --model_ont $D/model_pe.xml -o out -t 2 --sam_format sam --no_qcmsa --no_strains > ont.log 2>&1; echo "ont rc=$?"
grep -E "Align the|settled|chunks" hifi.log ont.log
for s in $HIFI $ONT; do
  echo "######## $s"
  ex=0.98; case $s in o*) ex=0.95;; esac
  python3 $S/eval_sam.py out/$s.sam $s.fq $s.truth.pkl --gtdb $G --db $D --exact $ex --show 4 | grep -vE "^mapq"
done
