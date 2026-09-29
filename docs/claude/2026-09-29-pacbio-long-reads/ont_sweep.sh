#!/usr/bin/env bash
# Error-rate sweep on the close-relatives database: native ONT alignment at -a 0.9 (default) and 0.85,
# and the reads' own error rate as the alignments show it (X, I, D per aligned base).
E=$(cd "$(dirname "$0")" && pwd)
B=${B:-$HOME/protal-lr-build}  # holds src/ (a checkout), build/protal, mini_db/, mini_db_close/
for err in 0.016 0.03 0.05; do
    echo "== error $err"
    for a in 0.9 0.85; do
        ERR=$err TAG=a$a EXTRA="-a $a" MINI=$B/mini_db_close bash $E/run_ont.sh 2>&1 | grep -E '^native|settled'
    done
    W=$B/ont_eval_mini_db_close_5_err$err
    awk -F'\t' '!/^@/ && !and($2, 256) { n = split($6, a, /[A-Z=]/); m = 0; e = 0;
        while (match($6, /[0-9]+[MXID=]/)) { v = substr($6, RSTART, RLENGTH - 1); op = substr($6, RSTART + RLENGTH - 1, 1);
            if (op == "M" || op == "=") m += v; else e += v; $6 = substr($6, RSTART + RLENGTH) }
        M += m; E += e } END { printf "differences per aligned base (errors + divergence): %.4f\n", E / (M + E) }' "$W/out/native_a0.9.sam"
done
