# Own read cluster features (own_cluster.py)

Seeds [1, 2, 3]; CV = species held out (5 folds) at 0.5; test = the independent test set; delta = test F1 against the base with a paired bootstrap (mean, 95% interval).

## pe

| features | CV F1 | test F1 | test FP | test FN | FN with a missing congener in the sample | delta (interval) |
|---|---|---|---|---|---|---|
| base | 0.9793 | 0.9690 | 18.7 | 74.0 | 39.3 | +0.0000 (+0.0000, +0.0000) |
| base + own cluster within 0.01 | 0.9803 | 0.9672 | 18.7 | 79.3 | 40.0 | -0.0019 (-0.0047, +0.0007) |
| base + own cluster within 0.02 | 0.9801 | 0.9694 | 16.7 | 74.7 | 38.3 | +0.0004 (-0.0020, +0.0031) |
| base + own cluster within 0.04 | 0.9796 | 0.9669 | 20.0 | 78.7 | 41.0 | -0.0021 (-0.0047, +0.0005) |
| base + identity shape | 0.9797 | 0.9695 | 16.7 | 74.3 | 38.7 | +0.0005 (-0.0017, +0.0027) |
| base + own within 0.01 + shape | 0.9805 | 0.9686 | 19.7 | 74.3 | 37.7 | -0.0005 (-0.0031, +0.0019) |

Borrowed reads (simulated from another species) among the present taxa's records, test set: median 0.010; with a missing congener in the sample 0.062 of all records and 0.004 of the own cluster's (0.01); absent taxa: 1.000.

## se

| features | CV F1 | test F1 | test FP | test FN | FN with a missing congener in the sample | delta (interval) |
|---|---|---|---|---|---|---|
| base | 0.9757 | 0.9646 | 16.0 | 86.7 | 40.3 | +0.0000 (+0.0000, +0.0000) |
| base + own cluster within 0.01 | 0.9752 | 0.9641 | 17.3 | 87.0 | 40.3 | -0.0005 (-0.0037, +0.0021) |
| base + own cluster within 0.02 | 0.9765 | 0.9647 | 18.0 | 84.7 | 41.7 | +0.0001 (-0.0027, +0.0028) |
| base + own cluster within 0.04 | 0.9751 | 0.9638 | 18.7 | 86.7 | 42.0 | -0.0009 (-0.0040, +0.0025) |
| base + identity shape | 0.9753 | 0.9639 | 16.3 | 88.3 | 42.0 | -0.0007 (-0.0039, +0.0021) |
| base + own within 0.01 + shape | 0.9742 | 0.9637 | 18.3 | 87.0 | 42.0 | -0.0009 (-0.0037, +0.0013) |

Borrowed reads (simulated from another species) among the present taxa's records, test set: median 0.018; with a missing congener in the sample 0.101 of all records and 0.014 of the own cluster's (0.01); absent taxa: 1.000.

