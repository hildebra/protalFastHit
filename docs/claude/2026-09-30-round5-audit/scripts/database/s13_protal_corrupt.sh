#!/usr/bin/env bash
# Protal-level corruption of a bundle + the partial-unpack shadow scenario. At most 2 threads.
set -u
P=~/strain-build/bin/protal
W=~/audit6/database
S=$(dirname "$0")
R1=$W/reads/sa_R1.fq; R2=$W/reads/sa_R2.fq
prof(){ /usr/bin/timeout 120 "$P" --db "$1" -1 $R1 -2 $R2 --prefix sa -o "$2" -t 2 --no_qcmsa --no_strains >"$3" 2>&1; echo $?; }

echo "=== A. flip a byte in the reference.fna frame of the bundle (offset 436618, len 86757)"
cp $W/bundle_only/database.protal $W/refflip.protal
python3 - <<PY
p="$W/refflip.protal"; d=bytearray(open(p,'rb').read())
d[436618+40000]^=0x5A; open(p,'wb').write(d)
print("flipped one byte inside reference.fna's frame")
PY
rc=$(prof $W/refflip.protal $W/out_refflip $W/refflip.log)
echo "  exit=$rc"; grep -iE 'error|corrupt|checksum|reference|cannot|invalid' $W/refflip.log | head -4 | sed 's/^/    /'

echo "=== B. truncate the bundle inside the reference frame (drop last 90000 bytes -> cuts model+ref tail+seek table)"
head -c $((530700)) $W/bundle_only/database.protal > $W/trunc.protal   # cut at model_pe start: removes model frame + seek table
rc=$(prof $W/trunc.protal $W/out_trunc $W/trunc.log)
echo "  exit=$rc"; grep -iE 'error|seek table|missing|corrupt|truncat|neither' $W/trunc.log | head -4 | sed 's/^/    /'

echo "=== C. shadow: a folder with index.prx.zst + database.protal but NO reference.fna"
rm -rf $W/shadow; mkdir -p $W/shadow
cp $W/folder/index.prx.zst $W/shadow/
cp $W/bundle_only/database.protal $W/shadow/database.protal
# (intentionally omit reference.fna / reference.map / taxonomy / model / unique_kmers)
rc=$(prof $W/shadow $W/out_shadow $W/shadow.log)
echo "  exit=$rc  (index.prx.zst present -> Locate picks separate-file mode, bundle ignored)"
grep -iE 'does not exist|Sequence file|neither|database:|Index file' $W/shadow.log | head -6 | sed 's/^/    /'

echo "=== D. shadow, complete: folder with ALL separate files + a bundle; both present -> which loads?"
rm -rf $W/shadow2; mkdir -p $W/shadow2
cp $W/folder/* $W/shadow2/
cp $W/bundle_only/database.protal $W/shadow2/database.protal
rc=$(prof $W/shadow2 $W/out_shadow2 $W/shadow2.log)
echo "  exit=$rc"; grep -iE 'Load index|database:|not /' $W/shadow2.log | head -4 | sed 's/^/    /'
