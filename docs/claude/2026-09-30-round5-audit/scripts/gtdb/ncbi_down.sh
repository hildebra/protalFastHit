#!/bin/bash
# NCBI unreachable: a datasets that always fails. How many requests, and what does the run report?
exec > ~/audit6/gtdb/ncbi_down.out 2>&1
PY=~/protal-train/bin/python
S=~/audit6/gtdb/src/scripts
D=~/audit6/gtdb/dl
rm -rf $D/i6; mkdir -p $D/i6
cp -a $D/inputs/release $D/inputs/download.json $D/i6/
cat > $D/datasets_down <<EOF
#!/bin/bash
[ "\$1" = download ] && echo "\$(grep -c . \$(echo "\$@" | sed 's/.*--inputfile \([^ ]*\).*/\1/'))" >> $D/down_calls.log
echo "Error: connection timed out" >&2; exit 1
EOF
chmod +x $D/datasets_down; : > $D/down_calls.log
# the release must resolve: serve the mirror again
$PY $(dirname "$0")/rangeserver.py ~/audit6/gtdb/w/mirror 18228 /dev/null 2>/dev/null &
SERVER=$!; sleep 1
$PY $S/download_gtdb.py -o $D/i6 --mirror http://127.0.0.1:18228 --datasets $D/datasets_down --species 30 --per_species 2 \
    --rep_only_species 20 --batch 16 -t 2 | tail -4
echo "exit ${PIPESTATUS[0]}"
echo "datasets download calls: $(wc -l < $D/down_calls.log); accessions per call: $(sort -n $D/down_calls.log | uniq -c | tr '\n' ' ')"
wc -l < $D/i6/missing.txt
kill $SERVER
