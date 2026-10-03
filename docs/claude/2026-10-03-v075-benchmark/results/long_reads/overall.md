# Means over all samples (score.py)

Per run as score.py scores it; archaea recall over all present archaea. Times: /usr/bin/time -v.

### full database

| reads | version | samples | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall | wall (s) | CPU (s) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pb | 0.7.3 | 8 | 0.8262 | 0.9899 | 0.7488 | 0.62 | 0.1735 | 0.5821 | 0.740 | 3.8 | 19.3 | 3.46 |
| pb | 0.7.5 | 8 | 0.8261 | 0.9994 | 0.7467 | 0.12 | 0.1728 | 0.5766 | 0.735 | 3.1 | 15.7 | 3.42 |
| ont | 0.7.3 | 8 | 0.8720 | 0.9889 | 0.8065 | 0.62 | 0.1600 | 0.6476 | 0.786 | 5.4 | 27.4 | 3.48 |
| ont | 0.7.5 | 8 | 0.8727 | 0.9973 | 0.8034 | 0.38 | 0.1607 | 0.6448 | 0.791 | 4.3 | 22.3 | 3.44 |

### missing database

| reads | version | samples | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall | wall (s) | CPU (s) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pb | 0.7.3 | 8 | 0.8001 | 0.9735 | 0.7230 | 0.88 | 0.1786 | 0.6240 | 0.698 | 3.2 | 16.4 | 3.33 |
| pb | 0.7.5 | 8 | 0.8114 | 0.9989 | 0.7263 | 0.12 | 0.1674 | 0.5774 | 0.695 | 2.6 | 13.0 | 3.27 |
| ont | 0.7.3 | 8 | 0.8581 | 0.9765 | 0.7903 | 1.12 | 0.1641 | 0.6522 | 0.719 | 4.6 | 24.4 | 3.37 |
| ont | 0.7.5 | 8 | 0.8596 | 0.9803 | 0.7922 | 0.75 | 0.1626 | 0.6465 | 0.744 | 3.9 | 20.2 | 3.33 |

