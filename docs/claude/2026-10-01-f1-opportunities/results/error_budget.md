# Error budget (error_budget.py)

Summed over the samples of each point; FN and FP by class; F1 at the knob, at the F1-optimal threshold of these samples (an oracle), and the ceiling if every present species with reads were called and nothing else.

| version | db | reads | point | TP | FN, no reads | FN, seen | FP, relative of a missing species | FP, congener of a present species | FP, other | F1 at 0.5 | best threshold: F1 | ceiling |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| v060 | full | pe | rl100_p1000 | 24 | 297 | 204 | 0 | 0 | 0 | 0.087 | 0.05: 0.562 | 0.606 |
| v060 | full | pe | rl100_p10000 | 167 | 14 | 186 | 0 | 0 | 0 | 0.625 | 0.09: 0.932 | 0.981 |
| v060 | full | pe | rl100_p500000 | 738 | 0 | 55 | 13 | 0 | 0 | 0.956 | 0.41: 0.989 | 1.000 |
| v060 | full | pe | rl150_p1000 | 56 | 249 | 211 | 0 | 0 | 0 | 0.196 | 0.05: 0.647 | 0.682 |
| v060 | full | pe | rl150_p10000 | 288 | 74 | 398 | 1 | 0 | 0 | 0.549 | 0.09: 0.909 | 0.949 |
| v060 | full | pe | rl150_p500000 | 323 | 0 | 17 | 9 | 0 | 0 | 0.961 | 0.39: 0.983 | 1.000 |
| v060 | full | pe | rl150_p5000000 | 353 | 0 | 7 | 3 | 0 | 0 | 0.986 | 0.44: 0.994 | 1.000 |
| v060 | missing | pe | rl100_p1000 | 18 | 186 | 135 | 0 | 0 | 0 | 0.101 | 0.05: 0.564 | 0.622 |
| v060 | missing | pe | rl100_p10000 | 98 | 5 | 119 | 0 | 0 | 0 | 0.613 | 0.13: 0.918 | 0.989 |
| v060 | missing | pe | rl100_p500000 | 467 | 0 | 32 | 18 | 0 | 1 | 0.948 | 0.39: 0.978 | 1.000 |
| v060 | missing | pe | rl150_p1000 | 33 | 155 | 130 | 0 | 0 | 0 | 0.188 | 0.05: 0.617 | 0.678 |
| v060 | missing | pe | rl150_p10000 | 179 | 44 | 246 | 1 | 0 | 0 | 0.552 | 0.11: 0.892 | 0.951 |
| v060 | missing | pe | rl150_p500000 | 189 | 0 | 9 | 14 | 0 | 0 | 0.943 | 0.41: 0.963 | 1.000 |
| v060 | missing | pe | rl150_p5000000 | 219 | 0 | 7 | 5 | 0 | 0 | 0.973 | 0.42: 0.987 | 1.000 |
| v070 | full | pe | rl100_p1000 | 235 | 246 | 44 | 2 | 0 | 0 | 0.617 | 0.05: 0.656 | 0.694 |
| v070 | full | pe | rl100_p10000 | 348 | 8 | 11 | 6 | 1 | 2 | 0.961 | 0.70: 0.964 | 0.989 |
| v070 | full | pe | rl100_p500000 | 789 | 0 | 4 | 10 | 9 | 9 | 0.980 | 0.83: 0.995 | 1.000 |
| v070 | full | pe | rl150_p1000 | 278 | 218 | 20 | 1 | 0 | 0 | 0.699 | 0.09: 0.707 | 0.732 |
| v070 | full | pe | rl150_p10000 | 694 | 42 | 24 | 11 | 0 | 2 | 0.946 | 0.36: 0.950 | 0.972 |
| v070 | full | pe | rl150_p500000 | 338 | 0 | 2 | 6 | 1 | 2 | 0.984 | 0.83: 0.996 | 1.000 |
| v070 | full | pe | rl150_p5000000 | 360 | 0 | 0 | 4 | 4 | 7 | 0.980 | 0.92: 0.994 | 1.000 |
| v070 | full | se | rl100_p1000 | 203 | 277 | 45 | 1 | 1 | 1 | 0.555 | 0.05: 0.611 | 0.642 |
| v070 | full | se | rl100_p10000 | 332 | 12 | 23 | 4 | 0 | 2 | 0.942 | 0.35: 0.952 | 0.983 |
| v070 | full | se | rl100_p500000 | 790 | 0 | 3 | 17 | 12 | 16 | 0.971 | 0.90: 0.992 | 1.000 |
| v070 | full | se | rl150_p1000 | 244 | 240 | 32 | 3 | 0 | 0 | 0.640 | 0.10: 0.674 | 0.697 |
| v070 | full | se | rl150_p10000 | 666 | 56 | 38 | 6 | 0 | 1 | 0.930 | 0.14: 0.936 | 0.962 |
| v070 | full | se | rl150_p500000 | 338 | 0 | 2 | 3 | 2 | 0 | 0.990 | 0.87: 0.996 | 1.000 |
| v070 | full | pb | pb_b3000000 | 455 | 39 | 22 | 5 | 0 | 0 | 0.932 | 0.29: 0.944 | 0.961 |
| v070 | full | pb | pb_b90000000 | 752 | 0 | 8 | 7 | 1 | 0 | 0.989 | 0.55: 0.992 | 1.000 |
| v070 | full | ont | ont_b3000000 | 479 | 11 | 26 | 9 | 1 | 1 | 0.952 | 0.19: 0.961 | 0.989 |
| v070 | full | ont | ont_b90000000 | 755 | 0 | 5 | 7 | 0 | 1 | 0.991 | 0.58: 0.995 | 1.000 |
| v070 | missing | pe | rl100_p1000 | 150 | 149 | 40 | 4 | 0 | 0 | 0.609 | 0.09: 0.653 | 0.718 |
| v070 | missing | pe | rl100_p10000 | 208 | 3 | 11 | 4 | 2 | 3 | 0.948 | 0.45: 0.950 | 0.993 |
| v070 | missing | pe | rl100_p500000 | 490 | 0 | 9 | 5 | 0 | 1 | 0.985 | 0.43: 0.985 | 1.000 |
| v070 | missing | pe | rl150_p1000 | 169 | 127 | 22 | 6 | 0 | 1 | 0.684 | 0.56: 0.691 | 0.750 |
| v070 | missing | pe | rl150_p10000 | 421 | 20 | 28 | 12 | 0 | 4 | 0.929 | 0.27: 0.938 | 0.978 |
| v070 | missing | pe | rl150_p500000 | 195 | 0 | 3 | 7 | 0 | 2 | 0.970 | 0.68: 0.985 | 1.000 |
| v070 | missing | pe | rl150_p5000000 | 223 | 0 | 3 | 1 | 0 | 0 | 0.991 | 0.52: 0.993 | 1.000 |
| v070 | missing | se | rl100_p1000 | 135 | 174 | 30 | 4 | 0 | 0 | 0.565 | 0.06: 0.619 | 0.655 |
| v070 | missing | se | rl100_p10000 | 199 | 5 | 18 | 2 | 1 | 3 | 0.932 | 0.25: 0.945 | 0.989 |
| v070 | missing | se | rl100_p500000 | 488 | 0 | 11 | 8 | 2 | 4 | 0.975 | 0.60: 0.980 | 1.000 |
| v070 | missing | se | rl150_p1000 | 142 | 141 | 35 | 6 | 0 | 0 | 0.609 | 0.10: 0.638 | 0.715 |
| v070 | missing | se | rl150_p10000 | 406 | 28 | 35 | 7 | 1 | 3 | 0.916 | 0.29: 0.925 | 0.969 |
| v070 | missing | se | rl150_p500000 | 195 | 0 | 3 | 4 | 3 | 1 | 0.973 | 0.87: 0.982 | 1.000 |
| v070 | missing | pb | pb_b3000000 | 282 | 23 | 13 | 5 | 0 | 1 | 0.931 | 0.32: 0.936 | 0.962 |
| v070 | missing | pb | pb_b90000000 | 456 | 0 | 13 | 12 | 0 | 2 | 0.971 | 0.52: 0.974 | 1.000 |
| v070 | missing | ont | ont_b3000000 | 284 | 9 | 25 | 8 | 0 | 3 | 0.927 | 0.26: 0.934 | 0.986 |
| v070 | missing | ont | ont_b90000000 | 454 | 0 | 15 | 1 | 1 | 1 | 0.981 | 0.43: 0.983 | 1.000 |
| v071 | full | pe | rl100_p1000 | 242 | 249 | 34 | 3 | 1 | 0 | 0.628 | 0.08: 0.667 | 0.689 |
| v071 | full | pe | rl100_p10000 | 348 | 8 | 11 | 3 | 1 | 1 | 0.967 | 0.48: 0.967 | 0.989 |
| v071 | full | pe | rl100_p500000 | 790 | 0 | 3 | 11 | 4 | 5 | 0.986 | 0.74: 0.996 | 1.000 |
| v071 | full | pe | rl150_p1000 | 284 | 218 | 14 | 2 | 0 | 0 | 0.708 | 0.13: 0.712 | 0.732 |
| v071 | full | pe | rl150_p10000 | 700 | 44 | 16 | 9 | 0 | 0 | 0.953 | 0.35: 0.955 | 0.970 |
| v071 | full | pe | rl150_p500000 | 339 | 0 | 1 | 5 | 0 | 4 | 0.985 | 0.94: 0.996 | 1.000 |
| v071 | full | pe | rl150_p5000000 | 360 | 0 | 0 | 2 | 1 | 3 | 0.992 | 0.71: 0.999 | 1.000 |
| v071 | full | se | rl100_p1000 | 224 | 277 | 24 | 0 | 1 | 1 | 0.597 | 0.08: 0.624 | 0.642 |
| v071 | full | se | rl100_p10000 | 339 | 12 | 16 | 3 | 0 | 1 | 0.955 | 0.44: 0.955 | 0.983 |
| v071 | full | se | rl100_p500000 | 790 | 0 | 3 | 10 | 2 | 6 | 0.987 | 0.72: 0.996 | 1.000 |
| v071 | full | se | rl150_p1000 | 258 | 240 | 18 | 3 | 0 | 0 | 0.664 | 0.07: 0.676 | 0.697 |
| v071 | full | se | rl150_p10000 | 675 | 56 | 29 | 5 | 0 | 2 | 0.936 | 0.31: 0.941 | 0.962 |
| v071 | full | se | rl150_p500000 | 338 | 0 | 2 | 4 | 2 | 6 | 0.980 | 0.82: 0.996 | 1.000 |
| v071 | full | pb | pb_b3000000 | 458 | 39 | 19 | 4 | 0 | 0 | 0.937 | 0.26: 0.945 | 0.961 |
| v071 | full | pb | pb_b90000000 | 753 | 1 | 6 | 6 | 0 | 0 | 0.991 | 0.52: 0.993 | 0.999 |
| v071 | full | ont | ont_b3000000 | 483 | 11 | 22 | 6 | 2 | 1 | 0.958 | 0.15: 0.964 | 0.989 |
| v071 | full | ont | ont_b90000000 | 759 | 0 | 1 | 5 | 1 | 0 | 0.995 | 0.63: 0.997 | 1.000 |
| v071 | missing | pe | rl100_p1000 | 158 | 152 | 29 | 4 | 0 | 1 | 0.629 | 0.21: 0.668 | 0.711 |
| v071 | missing | pe | rl100_p10000 | 208 | 3 | 11 | 3 | 1 | 3 | 0.952 | 0.50: 0.952 | 0.993 |
| v071 | missing | pe | rl100_p500000 | 489 | 0 | 10 | 2 | 0 | 0 | 0.988 | 0.43: 0.992 | 1.000 |
| v071 | missing | pe | rl150_p1000 | 171 | 127 | 20 | 6 | 0 | 1 | 0.690 | 0.11: 0.693 | 0.750 |
| v071 | missing | pe | rl150_p10000 | 424 | 22 | 23 | 8 | 1 | 3 | 0.937 | 0.30: 0.942 | 0.976 |
| v071 | missing | pe | rl150_p500000 | 195 | 0 | 3 | 5 | 1 | 0 | 0.977 | 0.69: 0.985 | 1.000 |
| v071 | missing | pe | rl150_p5000000 | 224 | 0 | 2 | 1 | 0 | 0 | 0.993 | 0.41: 0.996 | 1.000 |
| v071 | missing | se | rl100_p1000 | 147 | 174 | 18 | 3 | 0 | 0 | 0.601 | 0.09: 0.635 | 0.655 |
| v071 | missing | se | rl100_p10000 | 201 | 5 | 16 | 1 | 1 | 3 | 0.939 | 0.46: 0.942 | 0.989 |
| v071 | missing | se | rl100_p500000 | 491 | 0 | 8 | 4 | 1 | 1 | 0.986 | 0.41: 0.988 | 1.000 |
| v071 | missing | se | rl150_p1000 | 149 | 141 | 28 | 7 | 0 | 0 | 0.629 | 0.10: 0.652 | 0.715 |
| v071 | missing | se | rl150_p10000 | 407 | 28 | 34 | 7 | 1 | 3 | 0.918 | 0.28: 0.932 | 0.969 |
| v071 | missing | se | rl150_p500000 | 195 | 0 | 3 | 1 | 0 | 1 | 0.987 | 0.72: 0.990 | 1.000 |
| v071 | missing | pb | pb_b3000000 | 282 | 23 | 13 | 3 | 0 | 1 | 0.934 | 0.30: 0.938 | 0.962 |
| v071 | missing | pb | pb_b90000000 | 454 | 1 | 14 | 7 | 0 | 1 | 0.975 | 0.44: 0.976 | 0.999 |
| v071 | missing | ont | ont_b3000000 | 285 | 9 | 24 | 6 | 0 | 2 | 0.933 | 0.35: 0.939 | 0.986 |
| v071 | missing | ont | ont_b90000000 | 457 | 0 | 12 | 4 | 0 | 1 | 0.982 | 0.54: 0.986 | 1.000 |

Present species seen but not called, by what they were simulated from:

- v060 full pe: read pairs simulated: 1 135, 11-100 271, 2-3 164, 4-10 413, >100 95; genome: reference 185, strain 2-4% 555, strain <2% 338; Archaea 271, Bacteria 807
- v060 missing pe: read pairs simulated: 1 93, 11-100 162, 2-3 104, 4-10 262, >100 57; genome: reference 108, strain 2-4% 353, strain <2% 217; Archaea 137, Bacteria 541
- v070 full ont: read pairs simulated: 1 17, 2-3 10, 4-10 4; genome: reference 4, strain 2-4% 21, strain <2% 6; Archaea 9, Bacteria 22
- v070 full pb: read pairs simulated: 1 17, 2-3 7, 4-10 6; genome: reference 2, strain 2-4% 22, strain <2% 6; Archaea 4, Bacteria 26
- v070 full pe: read pairs simulated: 1 42, 11-100 7, 2-3 22, 4-10 28, >100 6; genome: reference 10, strain 2-4% 78, strain <2% 17; Archaea 19, Bacteria 86
- v070 full se: read pairs simulated: 1 43, 11-100 14, 2-3 26, 4-10 55, >100 5; genome: reference 10, strain 2-4% 101, strain <2% 32; Archaea 25, Bacteria 118
- v070 missing ont: read pairs simulated: 1 19, 11-100 1, 2-3 8, 4-10 12; genome: reference 6, strain 2-4% 24, strain <2% 10; Archaea 11, Bacteria 29
- v070 missing pb: read pairs simulated: 1 11, 2-3 6, 4-10 9; genome: reference 2, strain 2-4% 19, strain <2% 5; Archaea 3, Bacteria 23
- v070 missing pe: read pairs simulated: 1 42, 11-100 5, 2-3 21, 4-10 33, >100 15; genome: reference 16, strain 2-4% 76, strain <2% 24; Archaea 18, Bacteria 98
- v070 missing se: read pairs simulated: 1 42, 11-100 17, 2-3 20, 4-10 39, >100 14; genome: reference 15, strain 2-4% 81, strain <2% 36; Archaea 20, Bacteria 112
- v071 full ont: read pairs simulated: 1 15, 2-3 8; genome: reference 3, strain 2-4% 20; Archaea 7, Bacteria 16
- v071 full pb: read pairs simulated: 1 13, 11-100 1, 2-3 7, 4-10 4; genome: reference 2, strain 2-4% 19, strain <2% 4; Archaea 3, Bacteria 22
- v071 full pe: read pairs simulated: 1 32, 11-100 8, 2-3 10, 4-10 25, >100 4; genome: reference 6, strain 2-4% 59, strain <2% 14; Archaea 11, Bacteria 68
- v071 full se: read pairs simulated: 1 22, 11-100 13, 2-3 15, 4-10 37, >100 5; genome: reference 6, strain 2-4% 72, strain <2% 14; Archaea 18, Bacteria 74
- v071 missing ont: read pairs simulated: 1 18, 2-3 10, 4-10 8; genome: reference 5, strain 2-4% 25, strain <2% 6; Archaea 10, Bacteria 26
- v071 missing pb: read pairs simulated: 1 12, 11-100 1, 2-3 5, 4-10 9; genome: reference 2, strain 2-4% 18, strain <2% 7; Archaea 2, Bacteria 25
- v071 missing pe: read pairs simulated: 1 32, 11-100 6, 2-3 14, 4-10 31, >100 15; genome: reference 14, strain 2-4% 63, strain <2% 21; Archaea 16, Bacteria 82
- v071 missing se: read pairs simulated: 1 28, 11-100 14, 2-3 17, 4-10 37, >100 11; genome: reference 12, strain 2-4% 62, strain <2% 33; Archaea 20, Bacteria 87

One threshold per read type (pooled over every depth), F1 at 0.5 and at the best:

- v060 full pe: 0.6916 at 0.5, 0.8517 at 0.14
- v060 missing pe: 0.6849 at 0.5, 0.8285 at 0.15
- v070 full ont: 0.9759 at 0.5, 0.9778 at 0.54
- v070 full pb: 0.9671 at 0.5, 0.9674 at 0.41
- v070 full pe: 0.8973 at 0.5, 0.8981 at 0.68
- v070 full se: 0.8659 at 0.5, 0.8665 at 0.47
- v070 missing ont: 0.9591 at 0.5, 0.9625 at 0.43
- v070 missing pb: 0.9553 at 0.5, 0.9558 at 0.52
- v070 missing pe: 0.8883 at 0.5, 0.8908 at 0.32
- v070 missing se: 0.8554 at 0.5, 0.8601 at 0.33
- v071 full ont: 0.9807 at 0.5, 0.9809 at 0.60
- v071 full pb: 0.9700 at 0.5, 0.9730 at 0.27
- v071 full pe: 0.9037 at 0.5, 0.9039 at 0.52
- v071 full se: 0.8789 at 0.5, 0.8798 at 0.44
- v071 missing ont: 0.9624 at 0.5, 0.9661 at 0.60
- v071 missing pb: 0.9590 at 0.5, 0.9597 at 0.47
- v071 missing pe: 0.8943 at 0.5, 0.8955 at 0.43
- v071 missing se: 0.8667 at 0.5, 0.8696 at 0.40
