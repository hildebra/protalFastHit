#!/bin/bash
# download_gtdb.py against a localhost mirror of the synthetic world (no internet).
exec > ~/audit6/gtdb/dl.out 2>&1
set -u
PY=~/protal-train/bin/python
S=~/audit6/gtdb/src/scripts
W=~/audit6/gtdb/w
T=$(dirname "$0")
D=~/audit6/gtdb/dl
rm -rf $D; mkdir -p $D
echo '{}' > $D/faults.json
$PY $T/rangeserver.py $W/mirror 18226 $D/faults.json 2> $D/server.log &
SERVER=$!
sleep 1
M=http://127.0.0.1:18226
export FAKE_TABLE=$W/gtdb/simulation/genomes.tsv
# one strain NCBI "no longer has": pick the first non-representative accession
SUP=$(awk -F'\t' 'NR>1 && $5=="f" {print $1; exit}' $FAKE_TABLE | sed 's/^GB_//;s/^RS_//')
export FAKE_SUPPRESSED=$SUP
echo "suppressed: $SUP"

cat > $D/datasets <<EOF
#!/bin/bash
echo "\$*" >> $D/datasets_calls.log
exec $W/datasets "\$@"
EOF
chmod +x $D/datasets
echo "=== D1: full download"
$PY $S/download_gtdb.py -o $D/inputs --mirror $M --datasets $D/datasets --species 30 --per_species 2 \
    --rep_only_species 20 --batch 8 -t 2; echo "exit $?"
echo "datasets calls:"; grep -c '^download' $D/datasets_calls.log
ls $D/inputs $D/inputs/release | head -30
ls $D/inputs/genomes | wc -l
wc -l $D/inputs/simulation_species.txt $D/inputs/genomes.tsv $D/inputs/missing.txt
cat $D/inputs/missing.txt
$PY -c "import json; d=json.load(open('$D/inputs/download.json')); print(json.dumps(d['genomes'])); print(d['release'])"
echo "=== D1b: rerun"
$PY $S/download_gtdb.py -o $D/inputs --mirror $M --datasets $D/datasets --species 30 --per_species 2 \
    --rep_only_species 20 --batch 8 -t 2; echo "exit $?"

echo "=== D2: complete .part left by a kill before the checksum (416 on resume)"
mkdir -p $D/i2/release
cp $W/mirror/release226/226.0/bac120_taxonomy_r226.tsv.gz $D/i2/release/bac120_taxonomy_r226.tsv.gz.part
$PY $S/download_gtdb.py -o $D/i2 --mirror $M --no_genomes; echo "exit $?"
ls -la $D/i2/release | head

echo "=== D3: connection closes mid-transfer (Content-Length not reached)"
echo '{"truncate": {"bac120_marker_genes_all_r226.tar.gz": 5000}}' > $D/faults.json
$PY $S/download_gtdb.py -o $D/i3 --mirror $M --no_genomes 2>&1 | tail -8; echo "exit ${PIPESTATUS[0]}"
ls -la $D/i3/release/genomic_files_all/
echo '{}' > $D/faults.json
echo "--- rerun without the fault"
$PY $S/download_gtdb.py -o $D/i3 --mirror $M --no_genomes 2>&1 | tail -12; echo "exit ${PIPESTATUS[0]}"
grep 'bac120_marker_genes_all' $D/server.log | tail -3

echo "=== D4: server ignores Range (200) on a partial .part"
mkdir -p $D/i4/release
head -c 100 $W/mirror/release226/226.0/bac120_taxonomy_r226.tsv.gz > $D/i4/release/bac120_taxonomy_r226.tsv.gz.part
echo '{"no_range": true}' > $D/faults.json
$PY $S/download_gtdb.py -o $D/i4 --mirror $M --no_genomes 2>&1 | tail -4; echo "exit ${PIPESTATUS[0]}"
echo '{}' > $D/faults.json

echo "=== D5: tar members outside the release folder"
rm -rf $D/mirror2; cp -r $W/mirror $D/mirror2
B2=$D/mirror2/release226/226.0
mkdir -p $D/evil/ar53_marker_genes_reps_r226
tar -xzf $B2/genomic_files_reps/ar53_marker_genes_reps_r226.tar.gz -C $D/evil
cd $D/evil
echo pwned > evil.txt
ln -s $D/outside_target linkdir
mkdir -p $D/outside_target
$PY - <<EOF
import tarfile, io
with tarfile.open("$B2/genomic_files_reps/ar53_marker_genes_reps_r226.tar.gz", "w:gz") as tar:
    tar.add("ar53_marker_genes_reps_r226", arcname="ar53_marker_genes_reps_r226")
    tar.add("evil.txt", arcname="../../../evil_dotdot.txt")
    tar.add("evil.txt", arcname="$D/evil_absolute.txt")
    info = tarfile.TarInfo("ar53_marker_genes_reps_r226/linkdir"); info.type = tarfile.SYMTYPE; info.linkname = "$D/outside_target"
    tar.addfile(info)
    tar.add("evil.txt", arcname="ar53_marker_genes_reps_r226/linkdir/through_symlink.txt")
EOF
cd $B2 && md5sum ./VERSION.txt ./*.gz ./genomic_files_reps/*.gz ./genomic_files_all/*.gz > MD5SUM.txt; cd -
kill $SERVER; sleep 0.5
$PY $T/rangeserver.py $D/mirror2 18227 $D/faults.json 2> $D/server2.log &
SERVER=$!
sleep 1
$PY $S/download_gtdb.py -o $D/i5 --mirror http://127.0.0.1:18227 --no_genomes 2>&1 | tail -6; echo "exit ${PIPESTATUS[0]}"
ls -la $D/i5 $D/i5/release/genomic_files_reps/ar53_marker_genes_reps_r226 | head -20
find $D -name 'evil_*' -o -name through_symlink.txt
ls -la $D/outside_target
tar --version | head -1
kill $SERVER
