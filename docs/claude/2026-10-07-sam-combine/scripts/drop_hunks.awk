# awk -v drop="414 416 546 2823" -f drop_hunks.awk patch: a zero-context patch without the hunks whose old start line
# is listed in `drop` (another session's hunks in a file this session also changed).
BEGIN { n = split(drop, d, " "); for (i = 1; i <= n; i++) skip[d[i]] = 1 }
/^@@ / {
    split($2, a, ",")
    start = substr(a[1], 2)
    keep = !(start in skip)
}
/^(diff|index|---|\+\+\+) / { print; next }
/^@@ / || keep { if (keep) print }
