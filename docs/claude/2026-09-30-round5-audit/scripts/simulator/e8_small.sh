#!/bin/bash
# E8: small re-checks: GeneCov0 column, profile_utils with a folder glob, dupname merge.
set -u
W=~/audit6/simulator
P=$W/e1/prot/profiles
awk -F'\t' 'FNR==1{delete c; for(i=1;i<=NF;i++) if($i ~ /^GeneCov[0-9]+$/) c[i]=$i; next}
            {for(i in c) if($i!=0) nz[c[i]]=1; for(i in c) all[c[i]]=1}
            END{s=""; for(k in all) if(!(k in nz)) s=s" "k; print "GeneCov columns 0 in every row of s_1 and s_2:"s}' $P/s_1.profile.log $P/s_2.profile.log
echo "=== profile_utils merge --input 'profiles/*'"
cd $W/e6
python3 ~/strain-build/src/scripts/protal_profile_utils merge --input "$HOME/audit5/accuracy/prot_A_P/profiles/*" > /dev/null 2> glob.err; echo "exit $?"; cut -c1-260 glob.err
python3 ~/strain-build/src/scripts/protal_profile_utils merge --input "$P/s_1.profile*" > /dev/null 2> glob2.err; echo "exit $?"; cut -c1-260 glob2.err
echo "=== profile_utils on the dupname profile (two taxa, one lineage)"
python3 ~/strain-build/src/scripts/protal_profile_utils merge --input $W/e5/out_dupname/r.profile
