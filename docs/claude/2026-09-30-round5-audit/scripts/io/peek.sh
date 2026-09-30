#!/bin/bash
cd ~/audit6/io/t6
grep -a -A22 'ERROR: AddressSanitizer' fablse.log | cut -c1-200
echo ====
grep -a -B2 -A12 'IndexCodec.h:203' pe_zst.log | head -30 | cut -c1-200
sed -n 195,210p ~/audit6/io/src/src/Hash/IndexCodec.h
