#!/bin/bash
# E6: protal_profile_utils and protal_map_utils on real outputs and on crafted read folders.
set -u
W=~/audit6/simulator
U=~/strain-build/src/scripts
PU=$U/protal_profile_utils
MU=$U/protal_map_utils
cmp $PU /mnt/c/Users/hildebra/Documents/locDev/protal/.claude/worktrees/strain-fixes/scripts/protal_profile_utils && echo "profile_utils == worktree"
cmp $MU /mnt/c/Users/hildebra/Documents/locDev/protal/.claude/worktrees/strain-fixes/scripts/protal_map_utils && echo "map_utils == worktree"
rm -rf $W/e6; mkdir -p $W/e6; cd $W/e6
P=~/audit5/accuracy/prot_A_P/profiles
echo "=== profile_utils merge --input '$P/*' (glob as a user would give it)"
python3 $PU merge --input "$P/*" > m_all.tsv 2> m_all.err; echo "exit $?"; cat m_all.err | cut -c1-300
echo "=== merge --input '$P/*.profile'"
python3 $PU merge --input "$P/*.profile" > m_prof.tsv 2> m_prof.err; echo "exit $?"; head -c 600 m_prof.tsv; echo; wc -l m_prof.tsv; cat m_prof.err
echo "--- rows vs distinct lineages in the inputs"
cat $P/*.profile | cut -f2 | sort -u | wc -l
echo "--- column sums (should be 1 per sample)"
awk -F'\t' 'NR>1{for(i=2;i<=NF;i++) s[i]+=$i} NR==1{n=NF} END{for(i=2;i<=n;i++) printf "%.6f ", s[i]; print ""}' m_prof.tsv | cut -c1-200
echo "=== merge of .profile.gz"
cp $P/A00.profile a.profile; gzip -c a.profile > b.profile.gz
python3 $PU merge --input a.profile b.profile.gz > m_gz.tsv 2>&1; echo "exit $?"; head -3 m_gz.tsv
echo "=== merge of an empty profile (no taxon passed) and one with taxa"
: > empty.profile
python3 $PU merge --input empty.profile a.profile > m_empty.tsv 2>&1; echo "exit $?"; head -3 m_empty.tsv
echo "=== same sample name in two folders"
mkdir -p d1 d2; cp $P/A00.profile d1/; cp $P/A01.profile d2/A00.profile
python3 $PU merge --input d1/A00.profile d2/A00.profile > /dev/null 2> dup.err; echo "exit $?"; cat dup.err
python3 $PU merge --resolve-samples --input d1/A00.profile d2/A00.profile 2>&1 | head -2
G=~/genes_study/cong/out/profiles
echo "=== genes_study profiles"
python3 $PU merge --input "$G/*.profile" > m_g.tsv 2> m_g.err; echo "exit $?"; head -2 m_g.tsv | cut -c1-300; cat m_g.err

echo "=== map_utils generate on a mixed read folder"
mkdir -p reads/runA reads/runB
for f in s1_R1.fq.gz s1_R2.fq.gz s2_1.fastq.gz s2_2.fastq.gz s3_S1_L001_R1_001.fastq.gz s3_S1_L001_R2_001.fastq.gz \
         se_only.fq.gz ont_reads.fastq.gz s4.R1.fq s4.R2.fq; do : > reads/runA/$f; done
for f in s1_R1.fq.gz s1_R2.fq.gz; do : > reads/runB/$f; done
python3 $MU generate --input reads/runA > gen_a.map 2> gen_a.err; echo "exit $?"; cat gen_a.map; cat gen_a.err
echo "--- two runs with the same file names"
python3 $MU generate --input reads/runA reads/runB > gen_ab.map 2> gen_ab.err; echo "exit $?"; cat gen_ab.map; cat gen_ab.err
python3 $MU generate --id-from-folder --input reads/runA reads/runB > gen_abf.map 2> gen_abf.err; echo "exit $?"; cat gen_abf.map; cat gen_abf.err

echo "=== validate / flatten / merge with single-end rows ('-' in SECOND, as protal documents)"
cat > se.map <<EOF
#INPUT_DIR	$W/e6/reads/runA
#OUTPUT_DIR	$W/e6/out
#SAMPLEID	FIRST	SECOND	PREFIX	READ_TYPE
se1	se_only.fq.gz	-	se1	se
ont1	ont_reads.fastq.gz	-	ont1	ont
pe1	s1_R1.fq.gz	s1_R2.fq.gz	pe1	pe
EOF
python3 $MU validate --map se.map; echo "validate exit $?"
python3 $MU flatten --map se.map > se_flat.map; echo "flatten exit $?"; cat se_flat.map
cat > se_nosecond.map <<EOF
#INPUT_DIR	$W/e6/reads/runA
#SAMPLEID	FIRST	PREFIX
se1	se_only.fq.gz	se1
EOF
python3 $MU validate --map se_nosecond.map; echo "validate (no SECOND column) exit $?"
cat > pe2.map <<EOF
#INPUT_DIR	$W/e6/reads/runB
#OUTPUT_DIR	$W/e6/out2
#SAMPLEID	FIRST	SECOND	PREFIX	READ_TYPE
pe2	s1_R1.fq.gz	s1_R2.fq.gz	pe2	pe
EOF
python3 $MU merge --map se.map pe2.map > merged.map 2> merged.err; echo "merge exit $?"; cat merged.map; cat merged.err
echo "=== the simulator's protal.meta"
python3 $MU validate --map $W/e1/sim/protal.meta; echo "validate exit $?"
