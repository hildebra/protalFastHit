# Stacking the quality-calibrated divergence and the conservation pattern (combo.py)

## pe

| features | CV F1 | test F1 | test FP | test FN | delta (interval) |
|---|---|---|---|---|---|
| base | 0.9793 | 0.9690 | 18.7 | 74.0 | +0.0000 (+0.0000, +0.0000) |
| + excess | 0.9800 | 0.9710 | 15.0 | 71.7 | +0.0021 (-0.0008, +0.0052) |
| + conservation pattern | 0.9815 | 0.9707 | 16.7 | 71.0 | +0.0017 (-0.0007, +0.0044) |
| + both | 0.9824 | 0.9717 | 14.3 | 70.3 | +0.0028 (-0.0001, +0.0058) |

## ont

| features | CV F1 | test F1 | test FP | test FN | delta (interval) |
|---|---|---|---|---|---|
| base | 0.9737 | 0.9492 | 10.0 | 56.3 | +0.0000 (+0.0000, +0.0000) |
| + excess | 0.9806 | 0.9534 | 5.7 | 55.0 | +0.0045 (-0.0028, +0.0136) |
| + conservation pattern | 0.9735 | 0.9513 | 9.0 | 54.7 | +0.0021 (-0.0019, +0.0073) |
| + both | 0.9791 | 0.9557 | 7.3 | 50.7 | +0.0068 (-0.0016, +0.0195) |

