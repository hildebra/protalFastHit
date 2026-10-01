# Errors on the independent test set (test_errors.py)

## pe: 76 false negatives, 19 false positives, F1 0.9682

- missed by fragments: 1 36 of 206 (17.5%), 11-100 6 of 296 (2.0%), 2-3 10 of 141 (7.1%), 4-10 10 of 181 (5.5%), >100 14 of 698 (2.0%)
- missed by simulated from the representative: 0.0 69 of 1245 (5.5%), 1.0 7 of 277 (2.5%)
- missed by domain: Archaea 10 of 144 (6.9%), Bacteria 66 of 1378 (4.8%)
- present taxa with no congener of theirs missing from the database: 1035, missed 36 (3.5%), of them with more than 10 fragments 2
- present taxa with a missing congener in the sample: 487, missed 40 (8.2%), of them with more than 10 fragments 18
- false positives by meta_relative_rank: class 3, family 4, genus 10, order 2
- false positives by meta_novel_congener: 0 10, 1 9
- false positives' fragments: 1 7, 11-100 4, 2-3 5, 4-10 2, >100 1

| design point | F1 at 0.5 | best threshold: F1 |
|---|---|---|
| rl100_p10000 | 0.9466 | 0.19: 0.9582 |
| rl100_p500 | 0.9231 | 0.08: 0.9497 |
| rl100_p500000 | 0.9812 | 0.26: 0.9860 |
| rl150_p10000 | 0.9333 | 0.83: 0.9859 |
| rl150_p500 | 0.9422 | 0.09: 0.9748 |
| rl150_p500000 | 0.9861 | 0.43: 0.9908 |
| rl250_p10000 | 0.9690 | 0.16: 0.9718 |
| rl250_p500 | 0.9425 | 0.08: 0.9670 |
| rl250_p500000 | 0.9928 | 0.41: 0.9928 |

## se: 88 false negatives, 16 false positives, F1 0.9642

- missed by fragments: 1 33 of 198 (16.7%), 11-100 7 of 311 (2.3%), 2-3 19 of 143 (13.3%), 4-10 12 of 186 (6.5%), >100 17 of 649 (2.6%)
- missed by simulated from the representative: 0.0 82 of 1224 (6.7%), 1.0 6 of 263 (2.3%)
- missed by domain: Archaea 12 of 136 (8.8%), Bacteria 76 of 1351 (5.6%)
- present taxa with no congener of theirs missing from the database: 1014, missed 46 (4.5%), of them with more than 10 fragments 4
- present taxa with a missing congener in the sample: 473, missed 42 (8.9%), of them with more than 10 fragments 20
- false positives by meta_relative_rank: class 1, family 1, genus 12, order 2
- false positives by meta_novel_congener: 0 6, 1 10
- false positives' fragments: 1 7, 11-100 5, 2-3 3, 4-10 1

| design point | F1 at 0.5 | best threshold: F1 |
|---|---|---|
| rl100_p10000_se | 0.9235 | 0.22: 0.9634 |
| rl100_p500000_se | 0.9780 | 0.33: 0.9860 |
| rl100_p500_se | 0.9150 | 0.06: 0.9518 |
| rl150_p10000_se | 0.9333 | 0.94: 1.0000 |
| rl150_p500000_se | 0.9815 | 0.36: 0.9885 |
| rl150_p500_se | 0.9442 | 0.05: 0.9662 |
| rl250_p10000_se | 0.9651 | 0.46: 0.9690 |
| rl250_p500000_se | 0.9928 | 0.22: 0.9952 |
| rl250_p500_se | 0.9630 | 0.15: 0.9767 |

## pb: 38 false negatives, 5 false positives, F1 0.9658

- missed by fragments: 1 25 of 105 (23.8%), 11-100 3 of 244 (1.2%), 2-3 6 of 117 (5.1%), 4-10 2 of 111 (1.8%), >100 2 of 68 (2.9%)
- missed by simulated from the representative: 0.0 34 of 517 (6.6%), 1.0 4 of 128 (3.1%)
- missed by domain: Archaea 8 of 62 (12.9%), Bacteria 30 of 583 (5.1%)
- present taxa with no congener of theirs missing from the database: 418, missed 24 (5.7%), of them with more than 10 fragments 1
- present taxa with a missing congener in the sample: 227, missed 14 (6.2%), of them with more than 10 fragments 4
- false positives by meta_relative_rank: family 1, genus 4
- false positives by meta_novel_congener: 0 1, 1 4
- false positives' fragments: 1 3, 11-100 2

| design point | F1 at 0.5 | best threshold: F1 |
|---|---|---|
| pb_b150000 | 0.9256 | 0.05: 0.9764 |
| pb_b3000000 | 0.9544 | 0.22: 0.9742 |
| pb_b90000000 | 0.9875 | 0.41: 0.9891 |

## ont: 54 false negatives, 11 false positives, F1 0.9503

- missed by fragments: 1 32 of 90 (35.6%), 11-100 5 of 263 (1.9%), 2-3 8 of 132 (6.1%), 4-10 6 of 99 (6.1%), >100 3 of 92 (3.3%)
- missed by simulated from the representative: 0.0 47 of 541 (8.7%), 1.0 7 of 135 (5.2%)
- missed by domain: Archaea 7 of 67 (10.4%), Bacteria 47 of 609 (7.7%)
- present taxa with no congener of theirs missing from the database: 433, missed 21 (4.8%), of them with more than 10 fragments 0
- present taxa with a missing congener in the sample: 243, missed 33 (13.6%), of them with more than 10 fragments 8
- false positives by meta_relative_rank: class 1, family 1, genus 8, order 1
- false positives by meta_novel_congener: 0 3, 1 8
- false positives' fragments: 1 6, 11-100 2, 2-3 2, >100 1

| design point | F1 at 0.5 | best threshold: F1 |
|---|---|---|
| ont_b150000 | 0.8550 | 0.05: 0.9470 |
| ont_b3000000 | 0.9673 | 0.33: 0.9674 |
| ont_b90000000 | 0.9798 | 0.68: 0.9843 |

