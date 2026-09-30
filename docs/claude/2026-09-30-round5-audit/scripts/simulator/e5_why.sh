#!/bin/bash
cd ~/audit6/simulator/e5
grep -v "^\[=*>" err_base.txt | head -30
grep -iv "^\[=*>" log_base.txt | head -60 | tail -40
ls -la r1.fq r2.fq; head -4 r1.fq
