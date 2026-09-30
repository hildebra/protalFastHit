#!/bin/bash
set -u
cd ~/audit6/simulator/e1
grep -v "^\s*$" protal.err | grep -iv "resident\|page\|context\|time\|size\|swaps\|socket\|signals\|file system\|percent\|exit status\|command being" | tail -30
grep -i "fail\|error\|warn" protal.log | head -20
tail -15 protal.log
ls prot/profiles prot/strains | head
