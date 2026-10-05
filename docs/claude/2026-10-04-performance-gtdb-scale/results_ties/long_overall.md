# Means over all samples (score.py)

Per run as score.py scores it; archaea recall over all present archaea. Times: /usr/bin/time -v.

### full database

| reads | version | samples | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall | wall (s) | CPU (s) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pb | 0.7.3 | 8 | 0.8261 | 0.9994 | 0.7467 | 0.12 | 0.1728 | 0.5766 | 0.735 | 2.4 | 12.7 | 3.41 |
| pb | 0.7.5 | 8 | 0.8261 | 0.9994 | 0.7467 | 0.12 | 0.1728 | 0.5766 | 0.735 | 2.3 | 12.5 | 3.41 |
| ont | 0.7.3 | 8 | 0.8724 | 0.9973 | 0.8029 | 0.38 | 0.1610 | 0.6436 | 0.786 | 5.5 | 29.8 | 3.45 |
| ont | 0.7.5 | 8 | 0.8724 | 0.9973 | 0.8029 | 0.38 | 0.1610 | 0.6436 | 0.786 | 5.0 | 27.0 | 3.45 |

### missing database

| reads | version | samples | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall | wall (s) | CPU (s) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pb | 0.7.3 | 8 | 0.8114 | 0.9989 | 0.7263 | 0.12 | 0.1674 | 0.5774 | 0.695 | 1.9 | 10.5 | 3.27 |
| pb | 0.7.5 | 8 | 0.8114 | 0.9989 | 0.7263 | 0.12 | 0.1674 | 0.5774 | 0.695 | 2.0 | 10.8 | 3.26 |
| ont | 0.7.3 | 8 | 0.8623 | 0.9869 | 0.7922 | 0.62 | 0.1620 | 0.6450 | 0.744 | 4.8 | 25.9 | 3.33 |
| ont | 0.7.5 | 8 | 0.8623 | 0.9869 | 0.7922 | 0.62 | 0.1620 | 0.6450 | 0.744 | 4.6 | 24.6 | 3.33 |

