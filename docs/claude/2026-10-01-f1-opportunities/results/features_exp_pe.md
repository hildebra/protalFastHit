# Feature, model and threshold experiments (features_exp.py)

0.7.1's pipeline tables (/home/falk/bench071/V071); seeds [1, 2, 3]; CV = species held out (5 folds), at 0.5; test = the independent test set; delta = test F1 against the base, paired bootstrap (mean, 95% interval).

## pe (training 7783 taxa, test 3795)

| variant | CV F1 | test F1 | test FP | test FN | delta (interval) |
|---|---|---|---|---|---|
| base (0.7.1's 28 features) | 0.9793 | 0.9690 | 18.7 | 74.0 | +0.0000 (+0.0000, +0.0000) |
| + sample | 0.9842 | 0.9520 | 7.7 | 132.3 | -0.0167 (-0.0393, +0.0026) |
| + congener | 0.9818 | 0.9675 | 16.7 | 80.3 | -0.0014 (-0.0054, +0.0023) |
| + conserv | 0.9815 | 0.9707 | 16.7 | 71.0 | +0.0017 (-0.0007, +0.0044) |
| + nbr | 0.9796 | 0.9681 | 16.7 | 78.3 | -0.0009 (-0.0031, +0.0012) |
| base, CV-optimal threshold | 0.9793 | 0.9677 | 16.3 | 80.0 | -0.0013 (-0.0033, +0.0003) |
| base, threshold per depth bin | 0.9793 | 0.9698 | 17.3 | 73.0 | +0.0009 (-0.0040, +0.0057) |
| base, 256 trees | 0.9797 | 0.9690 | 18.0 | 74.7 | -0.0000 (-0.0017, +0.0017) |
| base, 512 leaves | 0.9790 | 0.9690 | 19.3 | 73.3 | +0.0000 (-0.0006, +0.0007) |
| base, gradient boosting | 0.9816 | 0.9679 | 20.0 | 76.0 | -0.0011 (-0.0047, +0.0022) |
| base, 25% of the training samples | 0.9736 | 0.9648 | 23.3 | 81.7 | -0.0042 (-0.0087, -0.0003) |
| base, 50% of the training samples | 0.9782 | 0.9647 | 17.0 | 88.0 | -0.0043 (-0.0100, +0.0001) |
| base, 75% of the training samples | 0.9776 | 0.9666 | 19.3 | 80.3 | -0.0024 (-0.0055, +0.0003) |

