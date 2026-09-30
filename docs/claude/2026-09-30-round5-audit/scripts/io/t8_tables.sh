#!/bin/bash
# T8: malformed reference.map / unique_kmers.tsv, parsed in parallel chunks (files > 1 MB): the
# first problem in the file is reported with its physical line number, for 1 and 2 threads.
set -u
source $(dirname "$0")/lib.sh
T=$W/t8; rm -rf $T; mkdir -p $T; cd $T
SAM=$W/t1/sam_t1/sa.sam
fna=$(stat -c %s $DB/reference.fna)
# a big map: the real 348 lines, then 70000 fake genes (taxid 100000+i) inside the reference
{ cat $DB/reference.map; for i in $(seq 1 70000); do printf '%d\t1\t%d\t%d\n' $((100000+i)) $((i % 1000)) $((i % 1000 + 50)); done; } > big.map
echo "big.map: $(wc -l < big.map) lines, $(stat -c %s big.map) bytes"
mkmapdb() { # mkmapdb DIR MAPFILE
  mkdir -p $1; for f in $DB/*; do [ "$(basename $f)" = reference.map ] || ln -sf $f $1/; done; cp $2 $1/reference.map
}
case_map() { # case_map NAME (map in $NAME.map)
  local name=$1
  mkmapdb db_$name $name.map
  for t in 1 2; do
    timeout 120 $P --db db_$name -1 $R/sa_R1.fq -2 $R/sa_R2.fq --prefix sa -o out_$name -t $t --no_qcmsa > $name.t$t.log 2>&1
    echo "$name -t $t rc=$? | $(grep -a -E 'Invalid|overlaps|different reference|rror' $name.t$t.log | head -1 | cut -c1-200)"
  done
}
L=45678
awk -v L=$L 'NR==L{print $1 "\t" $2 "\t" $3; next}{print}' big.map > cols3.map; case_map cols3
awk -v L=$L 'NR==L{print $1 "\tx" $2 "\t" $3 "\t" $4; next}{print}' big.map > nonnum.map; case_map nonnum
awk -v L=$L 'NR==L{print $1 "\t" $2 "\t" $4 "\t" $3; next}{print}' big.map > endstart.map; case_map endstart
awk -v L=$L -v F=$fna 'NR==L{print $1 "\t" $2 "\t" $3 "\t" F+1; next}{print}' big.map > pastend.map; case_map pastend
awk -v L=$L 'NR==L{print "100005\t1\t10\t20"; next}{print}' big.map > dup.map; case_map dup
awk -v L=$L 'NR==L{print $1 "\t" $2 "\t" $3; next}{print}' big.map | sed 's/$/\r/' > crlf.map; case_map crlf
awk -v L=$L 'NR==200{print ""} NR==30000{print ""; print ""} NR==L{print $1 "\t" $2 "\t" $3; next}{print}' big.map > blanks.map; case_map blanks   # 3 blank lines before: expect line L+3
awk 'NR==10000 || NR==60000{print $1 "\t" $2 "\t" $3; next}{print}' big.map > two.map; case_map two     # expect 10000
{ cat big.map; printf '5\t1\t0'; } > lastnonl.map; case_map lastnonl            # expect line 70349
cp big.map ok.map; case_map ok                                                  # parses; stops later (overlaps/fingerprint)
: > empty.map; case_map empty
printf '1\t1\t5\t410' > onenonl.map; case_map onenonl
echo "--- unique_kmers.tsv (--profile_only, no index)"
U=$DB/unique_kmers.tsv
for k in $(seq 1 100); do cat $U; done > bigu.tsv
echo "bigu.tsv: $(wc -l < bigu.tsv) lines, $(stat -c %s bigu.tsv) bytes"
mkudb() { mkdir -p $1; for f in $DB/*; do [ "$(basename $f)" = unique_kmers.tsv ] || ln -sf $f $1/; done; cp $2 $1/unique_kmers.tsv; }
case_u() {
  local name=$1
  mkudb udb_$name $name.tsv
  for t in 1 2; do
    timeout 120 $P --db udb_$name --profile_only $SAM --prefix u -o uout_${name}_$t -t $t --no_qcmsa --no_strains > u_$name.t$t.log 2>&1
    echo "u_$name -t $t rc=$? | $(grep -a -E 'Invalid|rror' u_$name.t$t.log | head -1 | cut -c1-200) | profile $(cmp -s uout_${name}_$t/u.profile $W/t1/sam_t1/sa.profile && echo SAME || echo diff/none)"
  done
}
L=23456
cp bigu.tsv uok.tsv; case_u uok
awk -v L=$L 'BEGIN{OFS="\t"} NR==L{$3="-1"; print; next}{print}' bigu.tsv > uneg.tsv; case_u uneg
awk -v L=$L 'BEGIN{OFS="\t"} NR==L{NF=8; print; next}{print}' bigu.tsv > ucols.tsv; case_u ucols
awk -v L=$L 'BEGIN{OFS="\t"} NR==L{$1="999"; print; next}{print}' bigu.tsv > unomap.tsv; case_u unomap
sed 's/$/\r/' uneg.tsv > ucrlf.tsv; case_u ucrlf
awk -v L=$L 'BEGIN{OFS="\t"} NR==L{$2="99999999999999999999999"; print; next}{print}' bigu.tsv > ubig.tsv; case_u ubig
: > uempty.tsv; case_u uempty
