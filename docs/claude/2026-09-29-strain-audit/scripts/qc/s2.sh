cd ~/audit5/smoke/protal
cat misc/s__Calidella_fervens.statistics.tsv
echo ---
head -3 strains/s__Calidella_fervens.meta.tsv
echo ---
cut -f1 strains/s__Calidella_fervens.meta.tsv | sort | uniq -c
echo --- snp_stats
cut -f1-5 strains/s__Calidella_fervens.snp_stats.tsv
echo --- partition
head -3 strains/s__Calidella_fervens.raw.partition.txt; wc -l strains/*.raw.partition.txt
echo --- misc hcov
cut -f1-5 misc/s__Calidella_fervens.hcov.tsv
awk -F'\t' '{print NF}' misc/s__Calidella_fervens.hcov.tsv | sort | uniq -c
grep -n "Calidella\|qcmsa\|No gene\|No good\|left out" ../protal.log | head -30
