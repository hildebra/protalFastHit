# The conservation pattern by class of taxon (feature_by_class.py)

0.7.1 pipeline, paired-end training samples. Only taxa with 5 or more hit genes in the second part.

| taxa | rows | ratio median (all) | ratio quartiles (all) | rows with >= 5 genes | ratio median (>= 5 genes) | conserved share median (>= 5 genes) |
|---|---|---|---|---|---|---|
| absent, congener of a held-out species in the sample | 1618 | +0.000 | +0.000 / +0.607 | 840 | +0.327 | 0.438 |
| absent, other | 3794 | +0.000 | +0.000 / +0.000 | 1030 | +0.098 | 0.833 |
| present, a congener held out in the sample | 536 | +0.000 | -0.235 / +0.202 | 407 | -0.023 | 0.463 |
| present, no congener held out | 1835 | +0.000 | -0.176 / +0.224 | 1495 | +0.008 | 0.471 |

genes with factor < 1: 84 of 168
