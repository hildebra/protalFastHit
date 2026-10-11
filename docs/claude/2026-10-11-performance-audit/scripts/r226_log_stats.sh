#!/usr/bin/env bash
# GTDB-scale numbers from the r226 builds' protal logs: start-up parts, SAM finishing, profiling steps.
# Usage: log_stats.sh LOG...
for f in "$@"; do
    echo "=== $f"
    echo "--- start-up per run (gene tables, preload, tables bracket, index load), seconds"
    awk '
        function secs(s,   t, m) { t = 0
            if (match(s, /([0-9]+)m /)) { t += 60 * substr(s, RSTART, RLENGTH - 2) }
            if (match(s, /([0-9]+)s /)) { t += substr(s, RSTART, RLENGTH - 2) }
            if (match(s, /([0-9]+)ms/)) { t += substr(s, RSTART, RLENGTH - 2) / 1000 }
            return t }
        /^Loading the gene tables took/ { gt = secs($0 " ") }
        /^Preload genomes took/ { pre = secs($0 " ") }
        /^Loading the taxonomy, models and tables took/ { tab = secs($0 " ") }
        /^Load Index took/ { idx = secs($0 " "); printf "gene_tables %.1f preload %.1f tables_bracket %.1f index %.1f\n", gt, pre, tab, idx }
    ' "$f"
    echo "--- SAM finishing: copied vs placed"
    grep -E "^SAM header:" "$f" | awk '
        /copied behind/ { n++; match($0, /records \(([0-9]+) bytes\)/); b += substr($0, RSTART + 9, RLENGTH - 16);
                          match($0, /behind it in [0-9.]+ s/); s += substr($0, RSTART + 13, RLENGTH - 15) }
        /placed in its room|written before the records/ { p++ }
        END { printf "copied %d (%.2f GB, %.1f s, %.0f MB/s), placed %d\n", n, b / 1e9, s, b / 1e6 / (s > 0 ? s : 1), p }'
    echo "--- aligning, SAM writing, profiling: total seconds over samples"
    awk '
        function secs(s,   t) { t = 0
            if (match(s, /([0-9]+)m /)) { t += 60 * substr(s, RSTART, RLENGTH - 2) }
            if (match(s, /([0-9]+)s /)) { t += substr(s, RSTART, RLENGTH - 2) }
            if (match(s, /([0-9]+)ms/)) { t += substr(s, RSTART, RLENGTH - 2) / 1000 }
            return t }
        /^Aligning reads took/ { al += secs($0 " "); na++ }
        /^Writing the SAM header and file took/ { sw += secs($0 " ") }
        /^Profiling sample .* took/ { np++
            match($0, /took [0-9.]+s/); pr += substr($0, RSTART + 5, RLENGTH - 6)
            rest = $0; sub(/^[^:]*: /, "", rest)
            n = split(rest, parts, /, /)
            for (i = 1; i <= n; i++) { if (match(parts[i], / [0-9.]+s$/)) { k = substr(parts[i], 1, RSTART - 1); step[k] += substr(parts[i], RSTART + 1, RLENGTH - 2) } } }
        END { printf "aligning %.0f s over %d samples; SAM writing %.0f s; profiling %.0f s over %d samples\n", al, na, sw, pr, np
              for (k in step) printf "  profiling step %-40s %.1f s\n", k, step[k] }
    ' "$f"
    echo "--- run totals"
    grep -E "^Run protal took|^Memory after (loading the index|aligning|profiling)" "$f" | head -20
done
