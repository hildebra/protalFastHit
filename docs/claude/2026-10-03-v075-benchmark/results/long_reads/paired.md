# Paired differences per sample (score.py)

Mean over the samples both versions ran, 95% bootstrap interval. wall: seconds per run.

### full database

| reads | comparison | samples | F1 | FP per sample | Bray-Curtis | wall (s) |
|---|---|---|---|---|---|---|
| pb | 0.7.5 less 0.7.3 | 8 | -0.0001 (-0.0010, +0.0009) | -0.50 (-1.25, +0.12) | -0.0007 (-0.0015, -0.0001) | -0.7 (-1.1, -0.2) |
| ont | 0.7.5 less 0.7.3 | 8 | +0.0007 (-0.0060, +0.0093) | -0.25 (-0.88, +0.25) | +0.0007 (-0.0007, +0.0022) | -1.1 (-1.9, -0.2) |

### missing database

| reads | comparison | samples | F1 | FP per sample | Bray-Curtis | wall (s) |
|---|---|---|---|---|---|---|
| pb | 0.7.5 less 0.7.3 | 8 | +0.0113 (+0.0004, +0.0246) | -0.75 (-1.38, -0.25) | -0.0112 (-0.0187, -0.0048) | -0.6 (-1.0, -0.1) |
| ont | 0.7.5 less 0.7.3 | 8 | +0.0015 (-0.0103, +0.0147) | -0.38 (-1.50, +0.50) | -0.0015 (-0.0073, +0.0061) | -0.7 (-1.3, -0.1) |

