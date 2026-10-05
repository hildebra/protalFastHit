# Means over all samples (score.py)

Per run as score.py scores it; archaea recall over all present archaea. Times: /usr/bin/time -v.

### full database

| reads | version | samples | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall | wall (s) | CPU (s) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pe | 0.7.3 | 24 | 0.9080 | 0.9883 | 0.8659 | 1.17 | 0.0893 | 0.2777 | 0.720 | 2.8 | 12.4 | 3.30 |
| pe | 0.7.5 | 24 | 0.9080 | 0.9883 | 0.8659 | 1.17 | 0.0893 | 0.2777 | 0.720 | 2.6 | 11.6 | 3.31 |

### missing database

| reads | version | samples | F1 | precision | recall | FP per sample | Bray-Curtis | median abs log2 | archaea recall | wall (s) | CPU (s) | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pe | 0.7.3 | 24 | 0.8974 | 0.9696 | 0.8609 | 1.25 | 0.1006 | 0.3075 | 0.679 | 2.2 | 9.7 | 3.21 |
| pe | 0.7.5 | 24 | 0.8974 | 0.9696 | 0.8609 | 1.25 | 0.1006 | 0.3075 | 0.679 | 2.1 | 9.5 | 3.20 |

