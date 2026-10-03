# Summary (score.py)

### F1 / Bray-Curtis, full database (mean over samples)

| reads | scenario | 0.7.3 | 0.7.5 |
|---|---|---|---|
| pb | pb_b3000000 | 0.667 / 0.265 | 0.667 / 0.264 |
| pb | pb_b90000000 | 0.986 / 0.082 | 0.985 / 0.082 |
| ont | ont_b3000000 | 0.753 / 0.238 | 0.754 / 0.239 |
| ont | ont_b90000000 | 0.991 / 0.082 | 0.991 / 0.082 |

### F1 / Bray-Curtis, missing database (mean over samples)

| reads | scenario | 0.7.3 | 0.7.5 |
|---|---|---|---|
| pb | pb_b3000000 | 0.621 / 0.273 | 0.643 / 0.253 |
| pb | pb_b90000000 | 0.979 / 0.085 | 0.979 / 0.082 |
| ont | ont_b3000000 | 0.742 / 0.234 | 0.740 / 0.238 |
| ont | ont_b90000000 | 0.974 / 0.094 | 0.979 / 0.087 |

### False positives per sample / precision / share of the reported abundance / congeners of a species in the sample, full database

| reads | scenario | 0.7.3 | 0.7.5 |
|---|---|---|---|
| pb | pb_b3000000 | 1.2 / 0.980 / 0.0052 / 5 of 5 | 0.0 / 1.000 / 0.0000 / 0 of 0 |
| pb | pb_b90000000 | 0.0 / 1.000 / 0.0000 / 0 of 0 | 0.2 / 0.999 / 0.0002 / 1 of 1 |
| ont | ont_b3000000 | 1.0 / 0.979 / 0.0025 / 4 of 4 | 0.2 / 0.997 / 0.0006 / 1 of 1 |
| ont | ont_b90000000 | 0.2 / 0.999 / 0.0000 / 1 of 1 | 0.5 / 0.998 / 0.0001 / 2 of 2 |

### False positives per sample / precision / share of the reported abundance / congeners of a species in the sample, missing database

| reads | scenario | 0.7.3 | 0.7.5 |
|---|---|---|---|
| pb | pb_b3000000 | 1.5 / 0.949 / 0.0104 / 6 of 6 | 0.0 / 1.000 / 0.0000 / 0 of 0 |
| pb | pb_b90000000 | 0.2 / 0.998 / 0.0003 / 1 of 1 | 0.2 / 0.998 / 0.0002 / 1 of 1 |
| ont | ont_b3000000 | 2.0 / 0.955 / 0.0037 / 7 of 8 | 0.8 / 0.967 / 0.0039 / 3 of 3 |
| ont | ont_b90000000 | 0.2 / 0.998 / 0.0149 / 1 of 1 | 0.8 / 0.994 / 0.0154 / 3 of 3 |

### Archaea: found of present (recall), archaeal false positives, full database

| reads | scenario | 0.7.3 | 0.7.5 |
|---|---|---|---|
| pb | pb_b3000000 | 36/81 = 0.44, FP 1 | 36/81 = 0.44, FP 0 |
| pb | pb_b90000000 | 109/115 = 0.95, FP 0 | 108/115 = 0.94, FP 0 |
| ont | ont_b3000000 | 42/81 = 0.52, FP 1 | 42/81 = 0.52, FP 0 |
| ont | ont_b90000000 | 112/115 = 0.97, FP 0 | 113/115 = 0.98, FP 0 |

### Archaea: found of present (recall), archaeal false positives, missing database

| reads | scenario | 0.7.3 | 0.7.5 |
|---|---|---|---|
| pb | pb_b3000000 | 18/42 = 0.43, FP 0 | 16/37 = 0.43, FP 0 |
| pb | pb_b90000000 | 49/54 = 0.91, FP 0 | 41/45 = 0.91, FP 0 |
| ont | ont_b3000000 | 21/42 = 0.50, FP 0 | 19/37 = 0.51, FP 0 |
| ont | ont_b90000000 | 48/54 = 0.89, FP 0 | 42/45 = 0.93, FP 0 |

### By true relative abundance: recall / median abs log2(reported / true) of those found, full database

| reads | scenario | abundance | 0.7.3 | 0.7.5 |
|---|---|---|---|---|
| pb | pb_b3000000 | 0.1-1% | 0.35 / 1.26 (n=445) | 0.34 / 1.26 (n=445) |
| pb | pb_b3000000 | >= 1% | 0.94 / 0.39 (n=71) | 0.94 / 0.39 (n=71) |
| pb | pb_b90000000 | 0.01-0.1% | 0.93 / 0.54 (n=370) | 0.93 / 0.54 (n=370) |
| pb | pb_b90000000 | 0.1-1% | 1.00 / 0.48 (n=321) | 1.00 / 0.48 (n=321) |
| pb | pb_b90000000 | >= 1% | 1.00 / 0.15 (n=69) | 1.00 / 0.15 (n=69) |
| ont | ont_b3000000 | 0.1-1% | 0.48 / 1.01 (n=445) | 0.47 / 1.01 (n=445) |
| ont | ont_b3000000 | >= 1% | 0.99 / 0.45 (n=71) | 0.97 / 0.44 (n=71) |
| ont | ont_b90000000 | 0.01-0.1% | 0.96 / 0.61 (n=370) | 0.96 / 0.60 (n=370) |
| ont | ont_b90000000 | 0.1-1% | 0.99 / 0.49 (n=321) | 1.00 / 0.49 (n=321) |
| ont | ont_b90000000 | >= 1% | 1.00 / 0.17 (n=69) | 1.00 / 0.17 (n=69) |

### By true relative abundance: recall / median abs log2(reported / true) of those found, missing database

| reads | scenario | abundance | 0.7.3 | 0.7.5 |
|---|---|---|---|---|
| pb | pb_b3000000 | 0.1-1% | 0.28 / 1.53 (n=254) | 0.27 / 1.66 (n=211) |
| pb | pb_b3000000 | >= 1% | 0.84 / 0.45 (n=64) | 0.85 / 0.43 (n=61) |
| pb | pb_b90000000 | 0.01-0.1% | 0.89 / 0.63 (n=159) | 0.87 / 0.64 (n=119) |
| pb | pb_b90000000 | 0.1-1% | 0.98 / 0.58 (n=259) | 0.98 / 0.57 (n=235) |
| pb | pb_b90000000 | >= 1% | 1.00 / 0.16 (n=51) | 1.00 / 0.17 (n=53) |
| ont | ont_b3000000 | 0.1-1% | 0.43 / 1.07 (n=254) | 0.41 / 1.04 (n=211) |
| ont | ont_b3000000 | >= 1% | 0.98 / 0.59 (n=64) | 0.97 / 0.72 (n=61) |
| ont | ont_b90000000 | 0.01-0.1% | 0.91 / 0.69 (n=159) | 0.93 / 0.74 (n=119) |
| ont | ont_b90000000 | 0.1-1% | 0.95 / 0.50 (n=259) | 0.96 / 0.49 (n=235) |
| ont | ont_b90000000 | >= 1% | 1.00 / 0.20 (n=51) | 1.00 / 0.20 (n=53) |

### Wall time / CPU time / peak RSS, full database (mean over samples)

| reads | scenario | 0.7.3 | 0.7.5 |
|---|---|---|---|
| pb | pb_b3000000 | 1.1 s / 4 s / 3.44 GB | 1.1 s / 4 s / 3.41 GB |
| pb | pb_b90000000 | 6.4 s / 35 s / 3.48 GB | 5.2 s / 28 s / 3.44 GB |
| ont | ont_b3000000 | 0.9 s / 3 s / 3.44 GB | 1.0 s / 4 s / 3.42 GB |
| ont | ont_b90000000 | 9.8 s / 52 s / 3.52 GB | 7.6 s / 41 s / 3.46 GB |

Load average before the runs: median 5.61, max 6.00.
