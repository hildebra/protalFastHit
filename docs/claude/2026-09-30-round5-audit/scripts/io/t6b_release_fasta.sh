#!/bin/bash
# T6b: the release binary on the FASTA inputs that trip ASan in StripString (blank lines, CRLF):
# does it read every record?
set -u
source $(dirname "$0")/lib.sh
T=$W/t6; cd $T; I=$T/in
head -n 800 $R/sa_R1.fq | awk 'NR%4==1{print ">" substr($0,2)} NR%4==2{print}' > $I/plain_R1.fa
head -n 800 $R/sa_R2.fq | awk 'NR%4==1{print ">" substr($0,2)} NR%4==2{print}' > $I/plain_R2.fa
for v in plain blank crlf; do
  prun rel_$v.log --db $DB -1 $I/${v}_R1.fa -2 $I/${v}_R2.fa --prefix rel_$v -o rel -t 2 --no_qcmsa --no_strains --sam_format sam > rel_$v.rc
  echo "$v $(cat rel_$v.rc) records=$(grep -vc '^@' rel/rel_$v.sam 2>/dev/null) qnames=$(grep -v '^@' rel/rel_$v.sam 2>/dev/null | cut -f1 | sort -u | wc -l)"
done
diff <(grep -v '^@' rel/rel_plain.sam | cut -f1-9 | sort) <(grep -v '^@' rel/rel_blank.sam | cut -f1-9 | sort) > /dev/null && echo "blank == plain (cols 1-9)"
diff <(grep -v '^@' rel/rel_plain.sam | cut -f1-9 | sort) <(grep -v '^@' rel/rel_crlf.sam | cut -f1-9 | sort) > /dev/null && echo "crlf == plain (cols 1-9)"
