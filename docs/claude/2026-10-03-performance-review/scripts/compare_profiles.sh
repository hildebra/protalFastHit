#!/usr/bin/env bash
# compare_profiles.sh A B [prefixA prefixB]: two protal output folders' profiles: taxa called in one and not the other
# (<prefix>.profile), and over all taxa of <prefix>.profile.log (joined by TaxID): rows whose Predicted differs, the
# largest |difference| of Probability, and whether anything but Predicted/Probability differs.
A=$1; B=$2; pa=${3:-s}; pb=${4:-$pa}
fa=$(find $A -name "$pa.profile" | head -1); fb=$(find $B -name "$pb.profile" | head -1)
la=$(find $A -name "$pa.profile.log" | head -1); lb=$(find $B -name "$pb.profile.log" | head -1)
ca=$(cut -f1 $fa | sort); cb=$(cut -f1 $fb | sort)
echo "   calls: $(echo "$ca" | wc -l) vs $(echo "$cb" | wc -l); only in A: $(comm -23 <(echo "$ca") <(echo "$cb") | wc -l), only in B: $(comm -13 <(echo "$ca") <(echo "$cb") | wc -l)"
awk -F'\t' 'NR==FNR { if (FNR>1) { pred[$10]=$1; prob[$10]=$2; $1=""; $2=""; rest[$10]=$0 } next }
  FNR>1 { n++; if ($10 in pred) { if (pred[$10]!=$1) dp++; d=prob[$10]-$2; if (d<0) d=-d; if (d>maxd) maxd=d; p1=$1; p2=$2; $1=""; $2=""; if (rest[$10]!=$0) other++ } else missing++ }
  END { printf "   profile.log: %d taxa; Predicted differs for %d; max |dProbability| %.6f; other columns differ for %d; taxa not in A: %d\n", n, dp, maxd, other, missing }' $la $lb
