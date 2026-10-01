# Where a missing species' reads land, by the gene's rate (trace_relatives.py)

36 paired-end training samples of the 0.7.1 pipeline; source coverage summed: own 2065.4, relative 473.8. Medians over the genes of each class.

| gene rate | genes | own records / coverage | relative records / coverage | R (relative / own) | on a congener | MAPQ < 4 | MAPQ < 10 | own MAPQ < 4 | R of the records kept (MAPQ >= 4) | kept records on the most-hit taxon | relative identity | own identity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rate < 0.7 | 58 | 2.275 | 2.470 | 1.012 | 0.770 | 0.696 | 0.834 | 0.106 | 0.376 | 0.891 | 0.958 | 0.983 |
| 0.7-1 | 33 | 4.602 | 3.341 | 0.772 | 0.924 | 0.462 | 0.654 | 0.025 | 0.426 | 0.826 | 0.936 | 0.977 |
| 1-1.4 | 41 | 4.512 | 2.322 | 0.539 | 0.976 | 0.348 | 0.529 | 0.011 | 0.364 | 0.810 | 0.928 | 0.972 |
| rate >= 1.4 | 36 | 2.485 | 0.538 | 0.192 | 0.998 | 0.233 | 0.371 | 0.005 | 0.141 | 0.835 | 0.919 | 0.964 |

Spearman correlation of the gene's rate with R over 168 genes: -0.975
