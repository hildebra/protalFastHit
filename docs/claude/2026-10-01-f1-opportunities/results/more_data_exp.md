# Twice the training data (more_data_exp.py)

Test F1 on the pipeline's independent test set, mean of seeds 1-3; delta against A with a paired bootstrap over test samples (mean, 95% interval).

## pe (A 7783 taxa, B 7699, test 3795)

| training | test F1 | FP | FN | delta against A (interval) |
|---|---|---|---|---|
| A (the pipeline's) | 0.9690 | 18.7 | 74.0 | +0.0000 (+0.0000, +0.0000) |
| B (seed 101) | 0.9667 | 13.7 | 85.3 | -0.0023 (-0.0075, +0.0021) |
| A + B | 0.9709 | 12.3 | 74.3 | +0.0021 (-0.0006, +0.0049) |

## se (A 7138 taxa, B 7095, test 3564)

| training | test F1 | FP | FN | delta against A (interval) |
|---|---|---|---|---|
| A (the pipeline's) | 0.9646 | 16.0 | 86.7 | +0.0000 (+0.0000, +0.0000) |
| B (seed 101) | 0.9616 | 13.0 | 98.0 | -0.0031 (-0.0088, +0.0010) |
| A + B | 0.9642 | 13.0 | 90.7 | -0.0004 (-0.0041, +0.0025) |

## pb (A 1085 taxa, B 1102, test 803)

| training | test F1 | FP | FN | delta against A (interval) |
|---|---|---|---|---|
| A (the pipeline's) | 0.9666 | 5.7 | 36.3 | +0.0000 (+0.0000, +0.0000) |
| B (seed 101) | 0.9654 | 7.3 | 36.3 | -0.0016 (-0.0088, +0.0055) |
| A + B | 0.9667 | 7.3 | 34.7 | -0.0000 (-0.0042, +0.0036) |

## ont (A 1424 taxa, B 1421, test 968)

| training | test F1 | FP | FN | delta against A (interval) |
|---|---|---|---|---|
| A (the pipeline's) | 0.9492 | 10.0 | 56.3 | +0.0000 (+0.0000, +0.0000) |
| B (seed 101) | 0.9545 | 6.3 | 53.0 | +0.0055 (-0.0038, +0.0171) |
| A + B | 0.9565 | 7.7 | 49.3 | +0.0074 (-0.0008, +0.0189) |

