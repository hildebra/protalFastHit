# hunks.awk -v keep="1,2": a unified diff of one file with only the listed hunks (1-based), header kept.
BEGIN { n = split(keep, k, ","); for (i = 1; i <= n; i++) want[k[i]] = 1; h = 0 }
/^diff --git/ || /^index / || /^--- / || /^\+\+\+ / { print; next }
/^@@/ { h++ }
{ if (h == 0 || (h in want)) print }
