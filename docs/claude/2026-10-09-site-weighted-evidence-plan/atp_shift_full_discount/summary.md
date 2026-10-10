# Ancestry true-positive test: /home/falk/atp_shift

Worlds: test (seed 11), train (seed 23); 36 test genera. See docs/claude/2026-10-09-ancestry-true-positive-test/README.md for the design and the predictions.

## Verdicts

| hypothesis | prediction | outcome | verdict |
|---|---|---|---|
| H1a table | stored alleles are the allele genomes' differences (>= 0.95 exact) | min share exact 1.000 | PASS |
| H3 sites | consensus sites = stem + representative's lineage (no recombination, >= 3 congeners): >= 0.98 of them, recall >= 0.95 | share 0.864, recall 0.651 | FAIL |
| H3 fallback | with 2 congeners the nearest congener's own derived states enter (share in history clearly below 1) | mean share 0.838 | FAIL |
| H7 novel reads | the allele scores draw some of a novel species' reads onto the target (it carries the alleles' ancestral base at the representative's own sites) | 5789 of 116100 reads (0.050) | INFO |
| H7 settling | the strains' reads the allele scores move go to their species (>= 0.9) | 0.505 of 6284 moved | FAIL |
| H1b/H2b oracle | protal's ancestry, allele and polymorphic features equal the oracle's (error-free reads, every row within 1e-6) | worst share within 1e-6: 1.000 | PASS |
| H1 leak | an allele genome's own reads: allele_explained_share >= 0.95, polymorphic_novel_share <= 0.02 | 0.952, 0.000 | PASS |
| H2 calibration | planted A_a: ancestry_agreement = a (+-0.03) | worst |diff| 0.018 | PASS |
| H1 calibration | planted P_f: allele_explained_share rises with f (a step near 0.5: a read is explained once it shares more than half of the allele's edits in its span) | f 0: 0.000, f 0.25: 0.072, f 0.5: 0.263, f 0.75: 0.531, f 1: 0.709 | PASS |
| H4 ancestry, clean strata | the same with >= 3 congeners and no recombination (the consensus applies) | b 0: 0.978, b 0.5: 0.821 | INFO |
| H5 fixed sites | deep strains with their representative-lineage sites covered (>= 0.8) against late novel species (b 0.8, 0.95), genera with alleles: ancestry_fixed_gain adds to ancestry_agreement (CI above 0) | AUC 0.616 -> 0.814, gain 0.197 [0.112, 0.286]; ancestry_fixed_agreement alone 0.818, the true fixed-site agreement alone 0.818 (159 Q, 285 N rows) | PASS |
| H4 ancestry | ancestry_agreement separates deep strains from novel species that left the stem early (b <= 0.5): AUC >= 0.9 | b 0: 0.935, b 0.5: 0.836 | FAIL |
| H6 confound | allele_explained_share favours strains with a stored relative (Q_near) over late twins, but not lone strains | Q_near 0.883, Q_lone 0.344 | PASS |
| H9 decoy | the decoy (a tip of the deep strains' clade labelled a species) against Q_deep: no feature separates (AUC 0.5 +- 0.1) | largest |AUC - 0.5|: uniqueness 0.542 (306 Q, 195 decoy rows) | PASS |
| H9b crown lineage | N_1, a novel lineage leaving at S's crown, against Q_deep: the ancestry features cannot tell them apart; the allele features can where Q_deep shares its side of the crown with allele genomes | ancestry_agreement 0.529, allele_explained_share 0.621 | INFO |
| H4/H5 gain | the three groups add separation over identity and depth (pe, CI above 0) | 0.119 [0.083, 0.152] | PASS |
| H5 ceiling | the true fixed-site agreement alone (the oracle's) adds separation; compare with '+all three' for the gap | 0.001 [-0.005, 0.006] | INFO |
| H5 ceiling, genera with alleles | the true fixed-site agreement (the oracle's) over identity and depth, where the fixed sites exist; compare with '+all three' | 0.043 [0.033, 0.054] | INFO |
| H8 calls pe target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.065 [0.037, 0.096] | PASS |
| H8 calls pe hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.119 [0.058, 0.177] | PASS |
| H8 calls hifi target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.060 [0.033, 0.086] | PASS |
| H8 calls hifi hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.161 [0.085, 0.224] | PASS |
| H9 shuffled | shuffled alleles cost the default model separation (the alleles carry information) | AUC 0.914 -> 0.834 | PASS |

## L0: the database's tables against the trees

Per test genus: copies of the target's genes with stored alleles, alleles per copy, and the share of stored alleles whose edits are exactly an allele genome's differences in the allele's range.
| genus | congeners | twin | allele_genomes | regime | copies | copies_with_alleles | alleles_per_copy | share_exact |
|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 0 | 0 | none | 120 | 0 | 0.000 | - |
| 2 | 2 | 0 | 0 | low | 120 | 0 | 0.000 | - |
| 3 | 2 | 0 | 0 | high | 120 | 0 | 0.000 | - |
| 4 | 2 | 0 | 2 | none | 120 | 120 | 1.975 | 1.000 |
| 5 | 2 | 0 | 2 | low | 120 | 119 | 1.917 | 1.000 |
| 6 | 2 | 0 | 2 | high | 120 | 118 | 1.892 | 1.000 |
| 7 | 2 | 0 | 8 | none | 120 | 118 | 3.100 | 1.000 |
| 8 | 2 | 0 | 8 | low | 120 | 119 | 3.900 | 1.000 |
| 9 | 2 | 0 | 8 | high | 120 | 116 | 3.758 | 1.000 |
| 10 | 2 | 1 | 0 | none | 120 | 0 | 0.000 | - |
| 11 | 2 | 1 | 0 | low | 120 | 0 | 0.000 | - |
| 12 | 2 | 1 | 0 | high | 120 | 0 | 0.000 | - |
| 13 | 2 | 1 | 2 | none | 120 | 118 | 1.917 | 1.000 |
| 14 | 2 | 1 | 2 | low | 120 | 119 | 1.883 | 1.000 |
| 15 | 2 | 1 | 2 | high | 120 | 116 | 1.775 | 1.000 |
| 16 | 2 | 1 | 8 | none | 120 | 120 | 3.758 | 1.000 |
| 17 | 2 | 1 | 8 | low | 120 | 117 | 3.825 | 1.000 |
| 18 | 2 | 1 | 8 | high | 120 | 117 | 3.742 | 1.000 |
| 19 | 6 | 0 | 0 | none | 120 | 0 | 0.000 | - |
| 20 | 6 | 0 | 0 | low | 120 | 0 | 0.000 | - |
| 21 | 6 | 0 | 0 | high | 120 | 0 | 0.000 | - |
| 22 | 6 | 0 | 2 | none | 120 | 120 | 1.750 | 1.000 |
| 23 | 6 | 0 | 2 | low | 120 | 118 | 1.925 | 1.000 |
| 24 | 6 | 0 | 2 | high | 120 | 116 | 1.808 | 1.000 |
| 25 | 6 | 0 | 8 | none | 120 | 120 | 3.992 | 1.000 |
| 26 | 6 | 0 | 8 | low | 120 | 120 | 3.933 | 1.000 |
| 27 | 6 | 0 | 8 | high | 120 | 117 | 3.708 | 1.000 |
| 28 | 6 | 1 | 0 | none | 120 | 0 | 0.000 | - |
| 29 | 6 | 1 | 0 | low | 120 | 0 | 0.000 | - |
| 30 | 6 | 1 | 0 | high | 120 | 0 | 0.000 | - |
| 31 | 6 | 1 | 2 | none | 120 | 118 | 1.900 | 1.000 |
| 32 | 6 | 1 | 2 | low | 120 | 94 | 1.142 | 1.000 |
| 33 | 6 | 1 | 2 | high | 120 | 109 | 1.558 | 1.000 |
| 34 | 6 | 1 | 8 | none | 120 | 84 | 1.367 | 1.000 |
| 35 | 6 | 1 | 8 | low | 120 | 115 | 3.575 | 1.000 |
| 36 | 6 | 1 | 8 | high | 120 | 114 | 2.717 | 1.000 |

Per strain with an allele genome in its species: the share of genes in which its nearest allele genome is among the stored alleles (K = 4, chosen for the sample's coverage since 2026-10-09 evening, farthest first before; the build rejects alleles as far as the nearest congener's copy).
| genome | genus | role | allele_genomes | twin | mrca_nearest_allele | share_genes_nearest_stored |
|---|---|---|---|---|---|---|
| GCA_998004004.1 | 4 | Q_deep | 2 | 0 | 0.000919 | 1.000 |
| GCA_998004005.1 | 4 | Q_deep | 2 | 0 | 0.000919 | 1.000 |
| GCA_998004017.1 | 4 | Q_ils | 2 | 0 | 0.000919 | 1.000 |
| GCA_998004013.1 | 4 | Q_imp0.2 | 2 | 0 | 0.000465 | 1.000 |
| GCA_998004015.1 | 4 | Q_imp0.4 | 2 | 0 | 0.000465 | 1.000 |
| GCA_998004003.1 | 4 | Q_lone | 2 | 0 | 0.003716 | 0.975 |
| GCA_998004001.1 | 4 | Q_near | 2 | 0 | 0.000465 | 1.000 |
| GCA_998004002.1 | 4 | Q_near | 2 | 0 | 0.000678 | 0.975 |
| GCA_998004006.1 | 4 | Q_rep | 2 | 0 | 0.003716 | 0.975 |
| GCA_998005004.1 | 5 | Q_deep | 2 | 0 | 0.004292 | 0.950 |
| GCA_998005005.1 | 5 | Q_deep | 2 | 0 | 0.004292 | 0.950 |
| GCA_998005017.1 | 5 | Q_ils | 2 | 0 | 0.004292 | 0.950 |
| GCA_998005013.1 | 5 | Q_imp0.2 | 2 | 0 | 0.000542 | 0.950 |
| GCA_998005015.1 | 5 | Q_imp0.4 | 2 | 0 | 0.000542 | 0.950 |
| GCA_998005003.1 | 5 | Q_lone | 2 | 0 | 0.007043 | 0.967 |
| GCA_998005001.1 | 5 | Q_near | 2 | 0 | 0.000542 | 0.950 |
| GCA_998005002.1 | 5 | Q_near | 2 | 0 | 9.2e-05 | 0.967 |
| GCA_998005006.1 | 5 | Q_rep | 2 | 0 | 0.000373 | 0.967 |
| GCA_998006004.1 | 6 | Q_deep | 2 | 0 | 0.008763 | 0.933 |
| GCA_998006005.1 | 6 | Q_deep | 2 | 0 | 0.008763 | 0.933 |
| GCA_998006017.1 | 6 | Q_ils | 2 | 0 | 0.008763 | 0.933 |
| GCA_998006013.1 | 6 | Q_imp0.2 | 2 | 0 | 0.001551 | 0.958 |
| GCA_998006015.1 | 6 | Q_imp0.4 | 2 | 0 | 0.001551 | 0.958 |
| GCA_998006003.1 | 6 | Q_lone | 2 | 0 | 0.004133 | 0.958 |
| GCA_998006001.1 | 6 | Q_near | 2 | 0 | 0.001551 | 0.958 |
| GCA_998006002.1 | 6 | Q_near | 2 | 0 | 4.3e-05 | 0.933 |
| GCA_998006006.1 | 6 | Q_rep | 2 | 0 | 0.004133 | 0.958 |
| GCA_998007004.1 | 7 | Q_deep | 8 | 0 | 0.011217 | 0.658 |
| GCA_998007005.1 | 7 | Q_deep | 8 | 0 | 0.011217 | 0.658 |
| GCA_998007017.1 | 7 | Q_ils | 8 | 0 | 0.011217 | 0.658 |
| GCA_998007013.1 | 7 | Q_imp0.2 | 8 | 0 | 0.000103 | 0.658 |
| GCA_998007015.1 | 7 | Q_imp0.4 | 8 | 0 | 0.000103 | 0.658 |
| GCA_998007001.1 | 7 | Q_near | 8 | 0 | 0.000103 | 0.658 |
| GCA_998007002.1 | 7 | Q_near | 8 | 0 | 2.3e-05 | 0.567 |
| GCA_998007003.1 | 7 | Q_near | 8 | 0 | 6.1e-05 | 0.442 |
| GCA_998007006.1 | 7 | Q_rep | 8 | 0 | 0.000784 | 0.442 |
| GCA_998008004.1 | 8 | Q_deep | 8 | 0 | 0.008927 | 0.567 |
| GCA_998008005.1 | 8 | Q_deep | 8 | 0 | 0.008927 | 0.567 |
| GCA_998008017.1 | 8 | Q_ils | 8 | 0 | 0.008927 | 0.567 |
| GCA_998008013.1 | 8 | Q_imp0.2 | 8 | 0 | 8e-06 | 0.567 |
(171 more rows in l0_nearest_allele_stored.tsv)

The build's line: Strain alleles: 7538 alleles (102032 edits) of 2762 gene copies of 24 species, up to 4 each (of 38880 full-reference copies of the database's species and genes: 0 of genomes outside --allele_genome_share 1, 26074 identical to the representative's, 2353 repeated, 19 covering less than 0.5 of it or more than 0.1 apart, 386 as far as the nearest congener's copy or farther; 10048 alleles offered): /home/falk/atp_shift/test/db/strain_alleles.tsv

The ancestry sites (the oracle's port of protal's rule, from the database's references; the oracle equals protal, L2a) against the trees: the share on S's stem or the representative's lineage, the recall of each, and why history sites are missed or other sites taken (by the congeners compared at the position).
| genus | congeners | twin | regime | sites | share_in_history | stem_recall | rep_lineage_recall | missed: < 3 congeners, the nearest agrees | missed: a congener carries another base | missed: not compared | other: consensus | other: the nearest's own (< 3 congeners) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 0 | none | 2529 | 0.875 | 0.697 | 0.673 | 925 |  | 83 |  | 317 |
| 2 | 2 | 0 | low | 4022 | 0.864 | 0.756 | 0.774 | 902 |  | 94 |  | 548 |
| 3 | 2 | 0 | high | 4221 | 0.712 | 0.768 | 0.690 | 867 |  | 82 |  | 1216 |
| 4 | 2 | 0 | none | 3301 | 0.875 | 0.697 | 0.678 | 1165 |  | 94 |  | 411 |
| 5 | 2 | 0 | low | 2772 | 0.875 | 0.695 | 0.677 | 1010 |  | 77 |  | 346 |
| 6 | 2 | 0 | high | 3790 | 0.857 | 0.720 | 0.743 | 1040 |  | 89 |  | 542 |
| 7 | 2 | 0 | none | 2832 | 0.883 | 0.708 | 0.704 | 943 |  | 65 |  | 330 |
| 8 | 2 | 0 | low | 4911 | 0.730 | 0.817 | 0.796 | 693 |  | 106 |  | 1326 |
| 9 | 2 | 0 | high | 5075 | 0.724 | 0.737 | 0.780 | 994 |  | 96 |  | 1400 |
| 10 | 2 | 1 | none | 1496 | 0.862 | 0.807 | 0.779 | 294 |  | 27 |  | 207 |
| 11 | 2 | 1 | low | 1904 | 0.685 | 0.839 | 0.861 | 195 |  | 33 |  | 599 |
| 12 | 2 | 1 | high | 1840 | 0.699 | 0.765 | 0.770 | 346 |  | 23 |  | 554 |
| 13 | 2 | 1 | none | 1324 | 0.834 | 0.795 | 0.810 | 244 |  | 20 |  | 220 |
| 14 | 2 | 1 | low | 1676 | 0.835 | 0.760 | 0.776 | 382 |  | 38 |  | 276 |
| 15 | 2 | 1 | high | 1551 | 0.799 | 0.718 | 0.692 | 464 |  | 30 |  | 312 |
| 16 | 2 | 1 | none | 1950 | 0.698 | 0.824 | 0.855 | 219 |  | 29 |  | 588 |
| 17 | 2 | 1 | low | 1742 | 0.691 | 0.836 | 0.854 | 178 |  | 28 |  | 539 |
| 18 | 2 | 1 | high | 1582 | 0.802 | 0.705 | 0.647 | 562 |  | 30 |  | 314 |
| 19 | 6 | 0 | none | 3761 | 0.910 | 0.756 | 0.744 | 12 | 1000 | 81 | 328 | 12 |
| 20 | 6 | 0 | low | 4189 | 0.893 | 0.751 | 0.745 | 25 | 1042 | 119 | 445 | 3 |
| 21 | 6 | 0 | high | 3309 | 0.926 | 0.659 | 0.598 | 17 | 1540 | 68 | 245 | 1 |
| 22 | 6 | 0 | none | 3038 | 0.867 | 0.795 | 0.792 | 11 | 600 | 52 | 403 | 1 |
| 23 | 6 | 0 | low | 3790 | 0.902 | 0.708 | 0.705 | 24 | 1210 | 97 | 369 | 4 |
| 24 | 6 | 0 | high | 3125 | 0.890 | 0.633 | 0.574 | 20 | 1546 | 83 | 341 | 3 |
| 25 | 6 | 0 | none | 4390 | 0.893 | 0.804 | 0.800 | 20 | 802 | 93 | 469 | 2 |
| 26 | 6 | 0 | low | 4410 | 0.871 | 0.784 | 0.781 | 24 | 857 | 134 | 562 | 5 |
| 27 | 6 | 0 | high | 3724 | 0.885 | 0.678 | 0.648 | 24 | 1419 | 76 | 428 |  |
| 28 | 6 | 1 | none | 1475 | 0.864 | 0.731 | 0.749 | 7 | 398 | 34 | 196 | 4 |
| 29 | 6 | 1 | low | 1593 | 0.842 | 0.808 | 0.826 | 8 | 255 | 25 | 250 | 2 |
| 30 | 6 | 1 | high | 1417 | 0.793 | 0.673 | 0.453 | 11 | 764 | 23 | 285 | 9 |
| 31 | 6 | 1 | none | 1144 | 0.889 | 0.651 | 0.708 | 6 | 426 | 31 | 127 |  |
| 32 | 6 | 1 | low | 1558 | 0.820 | 0.799 | 0.808 | 3 | 282 | 21 | 280 | 1 |
| 33 | 6 | 1 | high | 1806 | 0.698 | 0.680 | 0.515 | 15 | 724 | 27 | 544 | 2 |
| 34 | 6 | 1 | none | 1360 | 0.871 | 0.712 | 0.662 | 6 | 464 | 32 | 169 | 6 |
| 35 | 6 | 1 | low | 1678 | 0.779 | 0.767 | 0.695 | 9 | 405 | 21 | 362 | 8 |
| 36 | 6 | 1 | high | 1224 | 0.809 | 0.656 | 0.544 | 7 | 584 | 29 | 233 | 1 |

## L1: the focus genomes' reads, with and without the allele scores

The paired-end focus reads (best record of each mate): the share on the target without and with the allele scores, and the reads the scores moved (Q: the target is right; N: no species is).
| class | twin | reads | on_target_noallele | on_target_default | moved | moved_to_target | share_moves_to_target |
|---|---|---|---|---|---|---|---|
| Q | 0 | 66509 | 0.939 | 0.934 | 1322 | 470 | 0.356 |
| Q | 1 | 65698 | 0.889 | 0.896 | 4962 | 2704 | 0.545 |
| N | 0 | 57611 | 0.830 | 0.811 | 5254 | 1878 | 0.357 |
| N | 1 | 58489 | 0.757 | 0.734 | 9323 | 3911 | 0.420 |
| class | twin | from | to | reads |
|---|---|---|---|---|
| D | 0 | congener | congener | 7 |
| D | 0 | congener | other | 1 |
| D | 0 | congener | sister | 2 |
| D | 0 | congener | target | 18 |
| D | 0 | other | other | 4 |
| D | 0 | sister | congener | 6 |
| D | 0 | sister | target | 63 |
| D | 0 | target | congener | 21 |
| D | 0 | target | other | 3 |
| D | 0 | target | sister | 136 |
| D | 1 | congener | congener | 2 |
| D | 1 | congener | other | 1 |
| D | 1 | congener | target | 11 |
| D | 1 | other | other | 1 |
| D | 1 | sister | congener | 4 |
| D | 1 | sister | other | 2 |
| D | 1 | sister | target | 353 |
| D | 1 | target | congener | 3 |
| D | 1 | target | other | 2 |
| D | 1 | target | sister | 448 |
| N | 0 | congener | congener | 112 |
| N | 0 | congener | other | 17 |
| N | 0 | congener | sister | 58 |
| N | 0 | congener | target | 133 |
| N | 0 | other | congener | 7 |
| N | 0 | other | other | 43 |
| N | 0 | other | sister | 6 |
| N | 0 | other | target | 4 |
| N | 0 | sister | congener | 140 |
| N | 0 | sister | other | 28 |
| N | 0 | sister | target | 1741 |
| N | 0 | target | congener | 416 |
| N | 0 | target | other | 92 |
| N | 0 | target | sister | 2457 |
| N | 1 | congener | congener | 29 |
| N | 1 | congener | sister | 34 |
| N | 1 | congener | target | 59 |
| N | 1 | other | other | 7 |
| N | 1 | other | sister | 1 |
| N | 1 | other | target | 1 |
(29 more rows in l1_moves.tsv)
The protal log's 'strain alleles:' lines (every sample of the default run): 343414 unsure reads took the shifts, 57840 of them to another species.

## L2a: protal's features against the oracle's recomputation

Every target row: |protal - oracle|. The oracle ports AncestrySites.h, StrainAlleles.h and the profiler's sums; a feature off here computes something else than its definition.
| run | set | feature | rows | max_abs_diff | share_within_1e-6 |
|---|---|---|---|---|---|
| default | exact | ancestry_sites_per_record | 864 | 0.000 | 1.000 |
| default | exact | ancestry_agreement | 864 | 0.000 | 1.000 |
| default | exact | ancestry_congener_share | 864 | 0.000 | 1.000 |
| default | exact | ancestry_indel_sites_per_record | 864 | 0.000 | 1.000 |
| default | exact | ancestry_indel_congener_share | 864 | 0.000 | 1.000 |
| default | exact | allele_explained_share | 864 | 0.000 | 1.000 |
| default | exact | allele_identity_gain | 864 | 0.000 | 1.000 |
| default | exact | polymorphic_known_share | 864 | 0.000 | 1.000 |
| default | exact | polymorphic_novel_share | 864 | 0.000 | 1.000 |
| default | exact | ancestry_fixed_gain | 864 | 0.000 | 1.000 |
| default | exact | ancestry_fixed_agreement | 864 | 0.000 | 1.000 |
| default | exact | conserved_mismatch_ratio | 864 | 0.000 | 1.000 |
| default | exact | conserved_mismatch_rate | 864 | 0.000 | 1.000 |
| default | exact | ancestry_agreement_weighted | 864 | 0.000 | 1.000 |
| default | exact | nonsynonymous_share | 864 | 0.000 | 1.000 |
| default | exact | nonsynonymous_conserved_rate | 864 | 0.000 | 1.000 |
| default | exact | allele_copy_share | 864 | 0.000 | 1.000 |
| default | exact | allele_sites_per_kb | 864 | 0.000 | 1.000 |
| default | exact | ancestry_fixed_share | 864 | 0.000 | 1.000 |
| default | exact | column_weight_coverage | 864 | 0.000 | 1.000 |
| default | pe | ancestry_sites_per_record | 3284 | 0.000 | 1.000 |
| default | pe | ancestry_agreement | 3284 | 0.000 | 1.000 |
| default | pe | ancestry_congener_share | 3284 | 0.000 | 1.000 |
| default | pe | ancestry_indel_sites_per_record | 3284 | 0.000 | 1.000 |
| default | pe | ancestry_indel_congener_share | 3284 | 0.000 | 1.000 |
| default | pe | allele_explained_share | 3284 | 0.000 | 1.000 |
| default | pe | allele_identity_gain | 3284 | 0.000 | 1.000 |
| default | pe | polymorphic_known_share | 3284 | 0.000 | 1.000 |
| default | pe | polymorphic_novel_share | 3284 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_gain | 3284 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_agreement | 3284 | 0.000 | 1.000 |
| default | pe | conserved_mismatch_ratio | 3284 | 0.000 | 1.000 |
| default | pe | conserved_mismatch_rate | 3284 | 0.000 | 1.000 |
| default | pe | ancestry_agreement_weighted | 3284 | 0.000 | 1.000 |
| default | pe | nonsynonymous_share | 3284 | 0.000 | 1.000 |
| default | pe | nonsynonymous_conserved_rate | 3284 | 0.000 | 1.000 |
| default | pe | allele_copy_share | 3284 | 0.000 | 1.000 |
| default | pe | allele_sites_per_kb | 3284 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_share | 3284 | 0.000 | 1.000 |
| default | pe | column_weight_coverage | 3284 | 0.000 | 1.000 |
(40 more rows in l2_oracle.tsv)

## L2b: the planted genomes and the leaky control (error-free reads)

Mean features of the planted genomes (P_f: a share f of the farthest allele genome's differences; A_a: the congeners' base at a share 1 - a of the consensus sites) and of Q_leak (an allele genome itself). truth_agreement: the genome's agreement at the consensus sites by position.
| role | rows | truth_agreement | ancestry_agreement | allele_explained_share | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | identity |
|---|---|---|---|---|---|---|---|---|
| A_a0 | 7 | 0.002 | 0.003 | 0.127 | 0.283 | 0.053 | 0.000 | 0.978 |
| A_a0.25 | 7 | 0.253 | 0.268 | 0.079 | 0.210 | 0.034 | 0.008 | 0.980 |
| A_a0.5 | 6 | 0.502 | 0.502 | 0.042 | 0.141 | 0.023 | -0.009 | 0.980 |
| A_a0.75 | 7 | 0.752 | 0.749 | 0.015 | 0.071 | 0.012 | -0.005 | 0.980 |
| A_a1 | 6 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0 | 7 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0.25 | 6 | 0.873 | 0.921 | 0.072 | 0.204 | 0.000 | 0.075 | 0.980 |
| P_f0.5 | 7 | 0.925 | 0.840 | 0.263 | 0.409 | 0.000 | 0.153 | 0.980 |
| P_f0.75 | 6 | 0.817 | 0.791 | 0.531 | 0.582 | 0.000 | 0.199 | 0.980 |
| P_f1 | 6 | 0.513 | 0.732 | 0.709 | 0.798 | 0.000 | 0.257 | 0.978 |
| Q_leak | 110 | 0.855 | 0.755 | 0.952 | 0.498 | 0.000 | 0.195 | 0.989 |

## L2c: separation, pe (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1689 | 1400 | 0.774 | 0.685 | 0.795 | 0.505 | 0.500 | 0.823 | 0.188 | 0.625 | 0.564 | 0.549 | 0.394 | 0.505 | 0.651 | 0.669 | 0.500 | 0.824 | 0.503 | 0.500 | 0.833 | 0.650 |
| Q_deep vs N_0 | 306 | 172 | 0.825 | 0.811 | 0.878 | 0.383 | 0.555 | 0.935 | 0.069 | 0.697 | 0.650 | 0.660 | 0.431 | 0.690 | 0.764 | 0.596 | 0.500 | 0.946 | 0.527 | 0.500 | 0.946 | 0.761 |
| Q_deep vs N_0.5 | 306 | 192 | 0.692 | 0.679 | 0.773 | 0.529 | 0.535 | 0.836 | 0.187 | 0.690 | 0.654 | 0.669 | 0.448 | 0.634 | 0.712 | 0.580 | 0.500 | 0.843 | 0.486 | 0.500 | 0.872 | 0.712 |
| Q_deep vs N_0.8 | 306 | 210 | 0.614 | 0.617 | 0.671 | 0.519 | 0.520 | 0.664 | 0.348 | 0.632 | 0.604 | 0.610 | 0.441 | 0.552 | 0.654 | 0.531 | 0.500 | 0.669 | 0.527 | 0.500 | 0.704 | 0.653 |
| Q_deep vs N_0.95 | 306 | 214 | 0.550 | 0.544 | 0.560 | 0.531 | 0.499 | 0.578 | 0.432 | 0.642 | 0.627 | 0.638 | 0.463 | 0.560 | 0.620 | 0.557 | 0.500 | 0.577 | 0.497 | 0.500 | 0.627 | 0.620 |
| Q_deep vs N_1 | 306 | 221 | 0.521 | 0.537 | 0.530 | 0.525 | 0.509 | 0.529 | 0.481 | 0.621 | 0.606 | 0.621 | 0.467 | 0.539 | 0.577 | 0.503 | 0.500 | 0.532 | 0.484 | 0.500 | 0.565 | 0.576 |
| Q_near vs N_0.95 | 329 | 214 | 0.757 | 0.651 | 0.744 | 0.569 | 0.482 | 0.774 | 0.231 | 0.883 | 0.777 | 0.739 | 0.414 | 0.610 | 0.899 | 0.684 | 0.500 | 0.766 | 0.493 | 0.500 | 0.907 | 0.899 |
| Q_lone vs N_0.95 | 271 | 214 | 0.760 | 0.668 | 0.784 | 0.476 | 0.460 | 0.785 | 0.216 | 0.344 | 0.323 | 0.309 | 0.332 | 0.303 | 0.353 | 0.701 | 0.500 | 0.784 | 0.503 | 0.500 | 0.652 | 0.351 |
| Q_rep vs N_0 | 199 | 172 | 0.996 | 0.922 | 0.986 | 0.382 | 0.531 | 0.999 | 0.006 | 0.485 | 0.397 | 0.336 | 0.378 | 0.418 | 0.695 | 0.893 | 0.500 | 0.998 | 0.488 | 0.500 | 0.989 | 0.690 |
| Q_imp vs N (b <= 0.95) | 393 | 788 | 0.800 | 0.728 | 0.829 | 0.504 | 0.505 | 0.855 | 0.160 | 0.678 | 0.608 | 0.580 | 0.375 | 0.516 | 0.661 | 0.687 | 0.500 | 0.853 | 0.529 | 0.500 | 0.858 | 0.660 |
| Q (no import) vs N_imp | 1296 | 391 | 0.800 | 0.650 | 0.817 | 0.511 | 0.487 | 0.870 | 0.145 | 0.609 | 0.544 | 0.545 | 0.388 | 0.512 | 0.668 | 0.681 | 0.500 | 0.872 | 0.503 | 0.500 | 0.886 | 0.667 |
| Q_ils vs N_0.95 | 191 | 214 | 0.505 | 0.513 | 0.513 | 0.527 | 0.500 | 0.467 | 0.533 | 0.607 | 0.600 | 0.593 | 0.456 | 0.530 | 0.507 | 0.520 | 0.500 | 0.469 | 0.460 | 0.500 | 0.507 | 0.507 |
| decoy: Q_deep vs decoy | 306 | 195 | 0.497 | 0.485 | 0.503 | 0.542 | 0.509 | 0.500 | 0.503 | 0.515 | 0.516 | 0.520 | 0.498 | 0.526 | 0.539 | 0.481 | 0.500 | 0.501 | 0.493 | 0.500 | 0.541 | 0.538 |
| allele genomes 0 | 559 | 472 | 0.764 | 0.689 | 0.789 | 0.482 | 0.510 | 0.832 | 0.179 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.673 | 0.500 | 0.831 | 0.503 | 0.500 | 0.832 | 0.500 |
| allele genomes 2 | 568 | 464 | 0.778 | 0.678 | 0.800 | 0.514 | 0.493 | 0.824 | 0.192 | 0.781 | 0.666 | 0.635 | 0.322 | 0.505 | 0.815 | 0.688 | 0.500 | 0.824 | 0.529 | 0.500 | 0.864 | 0.815 |
| allele genomes 8 | 562 | 464 | 0.780 | 0.689 | 0.796 | 0.522 | 0.495 | 0.827 | 0.181 | 0.758 | 0.601 | 0.582 | 0.295 | 0.499 | 0.843 | 0.650 | 0.500 | 0.827 | 0.478 | 0.500 | 0.868 | 0.843 |
| congeners 2 | 841 | 690 | 0.783 | 0.680 | 0.803 | 0.516 | 0.507 | 0.818 | 0.198 | 0.633 | 0.576 | 0.542 | 0.388 | 0.503 | 0.641 | 0.664 | 0.500 | 0.819 | 0.503 | 0.500 | 0.822 | 0.641 |
| congeners 6 | 848 | 710 | 0.770 | 0.690 | 0.789 | 0.494 | 0.494 | 0.829 | 0.178 | 0.618 | 0.553 | 0.557 | 0.400 | 0.505 | 0.660 | 0.690 | 0.500 | 0.830 | 0.504 | 0.500 | 0.843 | 0.659 |
| twin 0 | 864 | 745 | 0.876 | 0.767 | 0.866 | 0.569 | 0.495 | 0.870 | 0.137 | 0.638 | 0.556 | 0.519 | 0.333 | 0.498 | 0.686 | 0.665 | 0.500 | 0.874 | 0.498 | 0.500 | 0.863 | 0.685 |
| twin 1 | 825 | 655 | 0.760 | 0.600 | 0.714 | 0.475 | 0.549 | 0.769 | 0.246 | 0.615 | 0.578 | 0.577 | 0.445 | 0.515 | 0.624 | 0.682 | 0.500 | 0.767 | 0.508 | 0.500 | 0.797 | 0.622 |
| regime low | 561 | 467 | 0.753 | 0.685 | 0.788 | 0.490 | 0.512 | 0.801 | 0.209 | 0.645 | 0.574 | 0.573 | 0.392 | 0.525 | 0.653 | 0.673 | 0.500 | 0.804 | 0.499 | 0.500 | 0.798 | 0.653 |
| regime high | 561 | 468 | 0.792 | 0.717 | 0.832 | 0.520 | 0.502 | 0.848 | 0.173 | 0.639 | 0.556 | 0.504 | 0.382 | 0.469 | 0.651 | 0.662 | 0.500 | 0.849 | 0.490 | 0.500 | 0.854 | 0.649 |
| regime none | 567 | 465 | 0.778 | 0.657 | 0.769 | 0.503 | 0.491 | 0.833 | 0.171 | 0.595 | 0.560 | 0.565 | 0.406 | 0.519 | 0.649 | 0.675 | 0.500 | 0.832 | 0.520 | 0.500 | 0.861 | 0.648 |
| depth 3 | 407 | 289 | 0.745 | 0.754 | 0.748 | 0.526 | 0.463 | 0.780 | 0.240 | 0.602 | 0.562 | 0.542 | 0.421 | 0.501 | 0.629 | 0.603 | 0.500 | 0.781 | 0.508 | 0.500 | 0.792 | 0.627 |
| depth 10 | 432 | 354 | 0.765 | 0.736 | 0.790 | 0.475 | 0.491 | 0.828 | 0.186 | 0.628 | 0.577 | 0.557 | 0.413 | 0.529 | 0.651 | 0.684 | 0.500 | 0.828 | 0.478 | 0.500 | 0.836 | 0.650 |
| depth 30 | 424 | 380 | 0.796 | 0.694 | 0.819 | 0.505 | 0.511 | 0.846 | 0.163 | 0.636 | 0.562 | 0.548 | 0.368 | 0.498 | 0.653 | 0.686 | 0.500 | 0.846 | 0.522 | 0.500 | 0.846 | 0.653 |
| depth 100 | 426 | 377 | 0.794 | 0.640 | 0.833 | 0.503 | 0.531 | 0.837 | 0.167 | 0.645 | 0.561 | 0.558 | 0.361 | 0.499 | 0.679 | 0.720 | 0.500 | 0.839 | 0.520 | 0.500 | 0.853 | 0.679 |
| all, within identity bands of 0.005 | 1689 | 1400 | 0.564 | 0.555 | 0.664 | 0.657 | 0.695 | 0.710 | 0.304 | 0.673 | 0.670 | 0.628 | 0.476 | 0.581 | 0.606 | 0.565 | 0.500 | 0.707 | 0.503 | 0.500 | 0.749 | 0.605 |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 156 | 84 | 0.789 | 0.866 | 0.924 | 0.722 | 0.668 | 0.787 | 0.936 | 0.787 |
| Q_deep vs N_0 | congeners 6 | 150 | 88 | 0.857 | 0.888 | 0.947 | 0.672 | 0.653 | 0.742 | 0.958 | 0.735 |
| Q_deep vs N_0 | no twin | 162 | 92 | 0.977 | 0.959 | 0.980 | 0.760 | 0.687 | 0.795 | 0.983 | 0.792 |
| Q_deep vs N_0 | twin | 144 | 80 | 0.919 | 0.843 | 0.895 | 0.636 | 0.630 | 0.738 | 0.904 | 0.734 |
| Q_deep vs N_0 | no allele genomes | 79 | 60 | 0.833 | 0.873 | 0.927 | 0.500 | 0.500 | 0.500 | 0.927 | 0.500 |
| Q_deep vs N_0 | allele genomes | 227 | 112 | 0.824 | 0.883 | 0.941 | 0.808 | 0.742 | 0.941 | 0.953 | 0.941 |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 19 | 0.799 | 0.809 | 0.978 | 0.773 | 0.798 | 0.944 | 0.996 | 0.945 |
| Q_deep vs N_0 | recombination none | 105 | 57 | 0.836 | 0.872 | 0.980 | 0.702 | 0.702 | 0.774 | 0.968 | 0.772 |
| Q_deep vs N_0 | recombination high | 102 | 56 | 0.813 | 0.846 | 0.918 | 0.715 | 0.641 | 0.781 | 0.960 | 0.774 |
| Q_deep vs N_0.5 | congeners 2 | 156 | 94 | 0.684 | 0.759 | 0.805 | 0.728 | 0.684 | 0.740 | 0.868 | 0.740 |
| Q_deep vs N_0.5 | congeners 6 | 150 | 98 | 0.702 | 0.786 | 0.863 | 0.650 | 0.655 | 0.686 | 0.883 | 0.684 |
| Q_deep vs N_0.5 | no twin | 162 | 100 | 0.933 | 0.920 | 0.907 | 0.727 | 0.674 | 0.795 | 0.916 | 0.795 |
| Q_deep vs N_0.5 | twin | 144 | 92 | 0.722 | 0.682 | 0.770 | 0.657 | 0.655 | 0.650 | 0.813 | 0.649 |
| Q_deep vs N_0.5 | no allele genomes | 79 | 69 | 0.710 | 0.755 | 0.817 | 0.500 | 0.500 | 0.500 | 0.817 | 0.500 |
| Q_deep vs N_0.5 | allele genomes | 227 | 123 | 0.687 | 0.783 | 0.854 | 0.795 | 0.744 | 0.837 | 0.910 | 0.837 |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 20 | 0.611 | 0.642 | 0.821 | 0.715 | 0.813 | 0.846 | 0.948 | 0.847 |
| Q_deep vs N_0.5 | recombination none | 105 | 61 | 0.725 | 0.782 | 0.879 | 0.685 | 0.712 | 0.740 | 0.909 | 0.739 |
| Q_deep vs N_0.5 | recombination high | 102 | 67 | 0.707 | 0.797 | 0.852 | 0.707 | 0.645 | 0.717 | 0.936 | 0.718 |
| Q_deep vs N_0.8 | congeners 2 | 156 | 108 | 0.625 | 0.686 | 0.655 | 0.666 | 0.611 | 0.691 | 0.719 | 0.691 |
| Q_deep vs N_0.8 | congeners 6 | 150 | 102 | 0.606 | 0.655 | 0.675 | 0.599 | 0.607 | 0.615 | 0.690 | 0.614 |
| Q_deep vs N_0.8 | no twin | 162 | 111 | 0.712 | 0.732 | 0.681 | 0.674 | 0.643 | 0.713 | 0.756 | 0.713 |
| Q_deep vs N_0.8 | twin | 144 | 99 | 0.680 | 0.626 | 0.648 | 0.594 | 0.572 | 0.616 | 0.635 | 0.615 |
| Q_deep vs N_0.8 | no allele genomes | 79 | 66 | 0.616 | 0.655 | 0.660 | 0.500 | 0.500 | 0.500 | 0.660 | 0.500 |
| Q_deep vs N_0.8 | allele genomes | 227 | 144 | 0.616 | 0.677 | 0.678 | 0.714 | 0.670 | 0.744 | 0.743 | 0.745 |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 24 | 0.559 | 0.615 | 0.568 | 0.653 | 0.756 | 0.633 | 0.690 | 0.639 |
| Q_deep vs N_0.8 | recombination none | 105 | 72 | 0.593 | 0.613 | 0.682 | 0.641 | 0.657 | 0.648 | 0.721 | 0.649 |
| Q_deep vs N_0.8 | recombination high | 102 | 73 | 0.632 | 0.739 | 0.697 | 0.626 | 0.588 | 0.679 | 0.765 | 0.679 |
| Q_deep vs N_0.95 | congeners 2 | 156 | 103 | 0.548 | 0.537 | 0.576 | 0.687 | 0.654 | 0.638 | 0.639 | 0.636 |
| Q_deep vs N_0.95 | congeners 6 | 150 | 111 | 0.554 | 0.582 | 0.585 | 0.598 | 0.627 | 0.596 | 0.610 | 0.597 |
| Q_deep vs N_0.95 | no twin | 162 | 115 | 0.565 | 0.590 | 0.594 | 0.657 | 0.658 | 0.650 | 0.662 | 0.649 |
| Q_deep vs N_0.95 | twin | 144 | 99 | 0.591 | 0.516 | 0.559 | 0.621 | 0.614 | 0.606 | 0.592 | 0.605 |
| Q_deep vs N_0.95 | no allele genomes | 79 | 73 | 0.570 | 0.551 | 0.598 | 0.500 | 0.500 | 0.500 | 0.598 | 0.500 |
| Q_deep vs N_0.95 | allele genomes | 227 | 141 | 0.550 | 0.567 | 0.579 | 0.715 | 0.705 | 0.656 | 0.636 | 0.657 |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 24 | 0.543 | 0.618 | 0.546 | 0.671 | 0.802 | 0.587 | 0.586 | 0.584 |
| Q_deep vs N_0.95 | recombination none | 105 | 68 | 0.571 | 0.553 | 0.599 | 0.668 | 0.696 | 0.641 | 0.644 | 0.640 |
| Q_deep vs N_0.95 | recombination high | 102 | 72 | 0.560 | 0.615 | 0.636 | 0.647 | 0.630 | 0.655 | 0.710 | 0.654 |

## L2c: separation, se (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1611 | 1303 | 0.753 | 0.691 | 0.753 | 0.486 | 0.518 | 0.796 | 0.218 | 0.619 | 0.571 | 0.548 | 0.417 | 0.500 | 0.629 | 0.635 | 0.500 | 0.795 | 0.506 | 0.500 | 0.810 | - |
| Q_deep vs N_0 | 285 | 157 | 0.772 | 0.792 | 0.818 | 0.423 | 0.606 | 0.904 | 0.109 | 0.689 | 0.651 | 0.637 | 0.482 | 0.636 | 0.709 | 0.596 | 0.500 | 0.910 | 0.534 | 0.500 | 0.912 | - |
| Q_deep vs N_0.5 | 285 | 180 | 0.652 | 0.672 | 0.710 | 0.520 | 0.560 | 0.790 | 0.237 | 0.676 | 0.650 | 0.647 | 0.454 | 0.595 | 0.689 | 0.568 | 0.500 | 0.793 | 0.478 | 0.500 | 0.814 | - |
| Q_deep vs N_0.8 | 285 | 192 | 0.585 | 0.595 | 0.629 | 0.539 | 0.510 | 0.640 | 0.372 | 0.633 | 0.616 | 0.607 | 0.458 | 0.559 | 0.626 | 0.561 | 0.500 | 0.641 | 0.512 | 0.500 | 0.693 | - |
| Q_deep vs N_0.95 | 285 | 206 | 0.515 | 0.533 | 0.531 | 0.545 | 0.514 | 0.541 | 0.457 | 0.648 | 0.636 | 0.631 | 0.482 | 0.547 | 0.599 | 0.526 | 0.500 | 0.542 | 0.511 | 0.500 | 0.592 | - |
| Q_deep vs N_1 | 285 | 207 | 0.516 | 0.526 | 0.497 | 0.550 | 0.512 | 0.494 | 0.511 | 0.615 | 0.600 | 0.604 | 0.472 | 0.532 | 0.545 | 0.519 | 0.500 | 0.493 | 0.501 | 0.500 | 0.533 | - |
| Q_near vs N_0.95 | 314 | 206 | 0.716 | 0.653 | 0.701 | 0.540 | 0.504 | 0.742 | 0.270 | 0.867 | 0.788 | 0.743 | 0.437 | 0.598 | 0.860 | 0.632 | 0.500 | 0.735 | 0.505 | 0.500 | 0.866 | - |
| Q_lone vs N_0.95 | 256 | 206 | 0.737 | 0.670 | 0.740 | 0.458 | 0.472 | 0.762 | 0.243 | 0.397 | 0.377 | 0.346 | 0.377 | 0.322 | 0.346 | 0.623 | 0.500 | 0.759 | 0.514 | 0.500 | 0.654 | - |
| Q_rep vs N_0 | 197 | 157 | 0.972 | 0.937 | 0.968 | 0.385 | 0.578 | 0.978 | 0.022 | 0.472 | 0.418 | 0.365 | 0.414 | 0.419 | 0.654 | 0.824 | 0.500 | 0.978 | 0.490 | 0.500 | 0.973 | - |
| Q_imp vs N (b <= 0.95) | 377 | 735 | 0.771 | 0.722 | 0.795 | 0.461 | 0.531 | 0.837 | 0.184 | 0.664 | 0.611 | 0.572 | 0.406 | 0.514 | 0.650 | 0.635 | 0.500 | 0.836 | 0.523 | 0.500 | 0.848 | - |
| Q (no import) vs N_imp | 1234 | 361 | 0.780 | 0.671 | 0.758 | 0.491 | 0.513 | 0.832 | 0.183 | 0.602 | 0.552 | 0.546 | 0.413 | 0.507 | 0.644 | 0.652 | 0.500 | 0.832 | 0.511 | 0.500 | 0.853 | - |
| Q_ils vs N_0.95 | 182 | 206 | 0.520 | 0.512 | 0.501 | 0.529 | 0.488 | 0.454 | 0.549 | 0.610 | 0.596 | 0.581 | 0.462 | 0.515 | 0.514 | 0.483 | 0.500 | 0.455 | 0.491 | 0.500 | 0.477 | - |
| decoy: Q_deep vs decoy | 285 | 186 | 0.487 | 0.496 | 0.508 | 0.562 | 0.513 | 0.488 | 0.506 | 0.507 | 0.510 | 0.504 | 0.520 | 0.505 | 0.515 | 0.538 | 0.500 | 0.486 | 0.508 | 0.500 | 0.507 | - |
| allele genomes 0 | 534 | 432 | 0.752 | 0.688 | 0.760 | 0.473 | 0.514 | 0.809 | 0.201 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.625 | 0.500 | 0.807 | 0.518 | 0.500 | 0.809 | - |
| allele genomes 2 | 545 | 439 | 0.759 | 0.688 | 0.746 | 0.480 | 0.531 | 0.794 | 0.229 | 0.746 | 0.649 | 0.617 | 0.362 | 0.501 | 0.786 | 0.651 | 0.500 | 0.793 | 0.524 | 0.500 | 0.829 | - |
| allele genomes 8 | 532 | 432 | 0.749 | 0.698 | 0.752 | 0.506 | 0.512 | 0.794 | 0.219 | 0.742 | 0.627 | 0.586 | 0.362 | 0.499 | 0.781 | 0.630 | 0.500 | 0.793 | 0.475 | 0.500 | 0.834 | - |
| congeners 2 | 803 | 637 | 0.762 | 0.691 | 0.753 | 0.494 | 0.524 | 0.794 | 0.221 | 0.629 | 0.583 | 0.545 | 0.412 | 0.494 | 0.620 | 0.630 | 0.500 | 0.794 | 0.501 | 0.500 | 0.799 | - |
| congeners 6 | 808 | 666 | 0.748 | 0.692 | 0.753 | 0.479 | 0.514 | 0.798 | 0.215 | 0.609 | 0.560 | 0.551 | 0.422 | 0.505 | 0.639 | 0.652 | 0.500 | 0.797 | 0.511 | 0.500 | 0.820 | - |
| twin 0 | 854 | 720 | 0.848 | 0.753 | 0.824 | 0.506 | 0.524 | 0.842 | 0.168 | 0.633 | 0.566 | 0.524 | 0.364 | 0.494 | 0.669 | 0.615 | 0.500 | 0.845 | 0.511 | 0.500 | 0.848 | - |
| twin 1 | 757 | 583 | 0.714 | 0.615 | 0.678 | 0.497 | 0.586 | 0.742 | 0.276 | 0.610 | 0.588 | 0.575 | 0.470 | 0.517 | 0.596 | 0.661 | 0.500 | 0.737 | 0.501 | 0.500 | 0.763 | - |
| regime low | 540 | 439 | 0.731 | 0.692 | 0.750 | 0.462 | 0.536 | 0.768 | 0.246 | 0.625 | 0.574 | 0.560 | 0.427 | 0.507 | 0.630 | 0.635 | 0.500 | 0.768 | 0.503 | 0.500 | 0.769 | - |
| regime high | 538 | 432 | 0.769 | 0.699 | 0.769 | 0.503 | 0.500 | 0.810 | 0.215 | 0.637 | 0.569 | 0.515 | 0.401 | 0.478 | 0.627 | 0.628 | 0.500 | 0.807 | 0.500 | 0.500 | 0.821 | - |
| regime none | 533 | 432 | 0.763 | 0.685 | 0.742 | 0.498 | 0.523 | 0.820 | 0.188 | 0.598 | 0.571 | 0.564 | 0.423 | 0.515 | 0.630 | 0.645 | 0.500 | 0.820 | 0.517 | 0.500 | 0.847 | - |
| depth 3 | 347 | 228 | 0.742 | 0.758 | 0.711 | 0.487 | 0.494 | 0.725 | 0.297 | 0.590 | 0.570 | 0.529 | 0.430 | 0.497 | 0.604 | 0.539 | 0.500 | 0.726 | 0.477 | 0.500 | 0.744 | - |
| depth 10 | 415 | 322 | 0.735 | 0.756 | 0.730 | 0.471 | 0.491 | 0.789 | 0.229 | 0.616 | 0.585 | 0.557 | 0.443 | 0.524 | 0.622 | 0.643 | 0.500 | 0.786 | 0.497 | 0.500 | 0.811 | - |
| depth 30 | 423 | 376 | 0.768 | 0.732 | 0.786 | 0.491 | 0.529 | 0.830 | 0.182 | 0.644 | 0.575 | 0.553 | 0.398 | 0.494 | 0.632 | 0.662 | 0.500 | 0.826 | 0.543 | 0.500 | 0.827 | - |
| depth 100 | 426 | 377 | 0.769 | 0.642 | 0.795 | 0.488 | 0.547 | 0.822 | 0.186 | 0.650 | 0.568 | 0.566 | 0.400 | 0.496 | 0.661 | 0.705 | 0.500 | 0.824 | 0.504 | 0.500 | 0.849 | - |
| all, within identity bands of 0.005 | 1611 | 1303 | 0.545 | 0.580 | 0.604 | 0.606 | 0.693 | 0.695 | 0.320 | 0.685 | 0.683 | 0.646 | 0.496 | 0.603 | 0.592 | 0.567 | 0.500 | 0.694 | 0.510 | 0.500 | 0.745 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 146 | 77 | 0.723 | 0.757 | 0.880 | 0.704 | 0.654 | 0.724 | 0.883 | - |
| Q_deep vs N_0 | congeners 6 | 139 | 80 | 0.818 | 0.878 | 0.928 | 0.670 | 0.622 | 0.695 | 0.941 | - |
| Q_deep vs N_0 | no twin | 158 | 88 | 0.932 | 0.901 | 0.953 | 0.744 | 0.674 | 0.782 | 0.954 | - |
| Q_deep vs N_0 | twin | 127 | 69 | 0.801 | 0.788 | 0.851 | 0.632 | 0.598 | 0.645 | 0.849 | - |
| Q_deep vs N_0 | no allele genomes | 77 | 52 | 0.804 | 0.825 | 0.905 | 0.500 | 0.500 | 0.500 | 0.905 | - |
| Q_deep vs N_0 | allele genomes | 208 | 105 | 0.763 | 0.818 | 0.905 | 0.791 | 0.706 | 0.848 | 0.920 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 32 | 18 | 0.701 | 0.804 | 0.958 | 0.750 | 0.727 | 0.819 | 0.989 | - |
| Q_deep vs N_0 | recombination none | 97 | 51 | 0.764 | 0.823 | 0.979 | 0.680 | 0.673 | 0.697 | 0.967 | - |
| Q_deep vs N_0 | recombination high | 95 | 50 | 0.782 | 0.799 | 0.875 | 0.720 | 0.640 | 0.749 | 0.898 | - |
| Q_deep vs N_0.5 | congeners 2 | 146 | 86 | 0.641 | 0.685 | 0.738 | 0.690 | 0.657 | 0.710 | 0.766 | - |
| Q_deep vs N_0.5 | congeners 6 | 139 | 94 | 0.667 | 0.735 | 0.839 | 0.662 | 0.639 | 0.668 | 0.864 | - |
| Q_deep vs N_0.5 | no twin | 158 | 97 | 0.838 | 0.855 | 0.843 | 0.720 | 0.663 | 0.764 | 0.872 | - |
| Q_deep vs N_0.5 | twin | 127 | 83 | 0.626 | 0.642 | 0.737 | 0.629 | 0.622 | 0.638 | 0.731 | - |
| Q_deep vs N_0.5 | no allele genomes | 77 | 62 | 0.674 | 0.745 | 0.775 | 0.500 | 0.500 | 0.500 | 0.775 | - |
| Q_deep vs N_0.5 | allele genomes | 208 | 118 | 0.643 | 0.690 | 0.810 | 0.773 | 0.719 | 0.814 | 0.838 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 32 | 20 | 0.593 | 0.555 | 0.828 | 0.718 | 0.783 | 0.809 | 0.898 | - |
| Q_deep vs N_0.5 | recombination none | 97 | 61 | 0.664 | 0.694 | 0.859 | 0.693 | 0.693 | 0.718 | 0.864 | - |
| Q_deep vs N_0.5 | recombination high | 95 | 59 | 0.680 | 0.748 | 0.814 | 0.699 | 0.628 | 0.689 | 0.873 | - |
| Q_deep vs N_0.8 | congeners 2 | 146 | 97 | 0.608 | 0.628 | 0.647 | 0.643 | 0.605 | 0.648 | 0.695 | - |
| Q_deep vs N_0.8 | congeners 6 | 139 | 95 | 0.567 | 0.636 | 0.632 | 0.622 | 0.611 | 0.603 | 0.693 | - |
| Q_deep vs N_0.8 | no twin | 158 | 106 | 0.652 | 0.679 | 0.660 | 0.680 | 0.644 | 0.688 | 0.740 | - |
| Q_deep vs N_0.8 | twin | 127 | 86 | 0.591 | 0.635 | 0.618 | 0.585 | 0.566 | 0.585 | 0.640 | - |
| Q_deep vs N_0.8 | no allele genomes | 77 | 59 | 0.651 | 0.669 | 0.669 | 0.500 | 0.500 | 0.500 | 0.669 | - |
| Q_deep vs N_0.8 | allele genomes | 208 | 133 | 0.560 | 0.613 | 0.637 | 0.722 | 0.678 | 0.712 | 0.715 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 32 | 20 | 0.473 | 0.623 | 0.589 | 0.718 | 0.766 | 0.678 | 0.781 | - |
| Q_deep vs N_0.8 | recombination none | 97 | 64 | 0.576 | 0.604 | 0.697 | 0.652 | 0.655 | 0.628 | 0.762 | - |
| Q_deep vs N_0.8 | recombination high | 95 | 67 | 0.565 | 0.635 | 0.614 | 0.656 | 0.607 | 0.618 | 0.702 | - |
| Q_deep vs N_0.95 | congeners 2 | 146 | 101 | 0.509 | 0.512 | 0.550 | 0.668 | 0.631 | 0.610 | 0.593 | - |
| Q_deep vs N_0.95 | congeners 6 | 139 | 105 | 0.521 | 0.548 | 0.535 | 0.625 | 0.632 | 0.586 | 0.590 | - |
| Q_deep vs N_0.95 | no twin | 158 | 114 | 0.524 | 0.532 | 0.560 | 0.676 | 0.661 | 0.654 | 0.650 | - |
| Q_deep vs N_0.95 | twin | 127 | 92 | 0.520 | 0.538 | 0.512 | 0.613 | 0.599 | 0.556 | 0.523 | - |
| Q_deep vs N_0.95 | no allele genomes | 77 | 71 | 0.533 | 0.570 | 0.568 | 0.500 | 0.500 | 0.500 | 0.568 | - |
| Q_deep vs N_0.95 | allele genomes | 208 | 135 | 0.520 | 0.521 | 0.543 | 0.724 | 0.699 | 0.627 | 0.595 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 32 | 22 | 0.586 | 0.578 | 0.516 | 0.675 | 0.740 | 0.663 | 0.575 | - |
| Q_deep vs N_0.95 | recombination none | 97 | 66 | 0.563 | 0.526 | 0.576 | 0.650 | 0.647 | 0.627 | 0.602 | - |
| Q_deep vs N_0.95 | recombination high | 95 | 70 | 0.503 | 0.563 | 0.556 | 0.668 | 0.637 | 0.623 | 0.668 | - |

## L2c: separation, hifi (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 576 | 513 | 0.817 | 0.784 | 0.880 | 0.593 | 0.504 | 0.883 | 0.121 | 0.671 | 0.613 | 0.546 | 0.308 | 0.507 | 0.683 | 0.445 | 0.500 | 0.888 | 0.481 | 0.500 | 0.857 | - |
| Q_deep vs N_0 | 101 | 75 | 0.815 | 0.837 | 0.947 | 0.596 | 0.500 | 0.998 | 0.004 | 0.745 | 0.730 | 0.729 | 0.493 | 0.709 | 0.765 | 0.439 | 0.500 | 0.999 | 0.527 | 0.500 | 0.988 | - |
| Q_deep vs N_0.5 | 101 | 74 | 0.782 | 0.820 | 0.955 | 0.571 | 0.457 | 0.964 | 0.043 | 0.695 | 0.693 | 0.684 | 0.344 | 0.629 | 0.746 | 0.430 | 0.500 | 0.978 | 0.573 | 0.500 | 0.904 | - |
| Q_deep vs N_0.8 | 101 | 65 | 0.651 | 0.705 | 0.787 | 0.584 | 0.508 | 0.767 | 0.233 | 0.678 | 0.692 | 0.695 | 0.401 | 0.577 | 0.687 | 0.474 | 0.500 | 0.778 | 0.494 | 0.500 | 0.742 | - |
| Q_deep vs N_0.95 | 101 | 73 | 0.588 | 0.666 | 0.625 | 0.527 | 0.471 | 0.610 | 0.387 | 0.649 | 0.672 | 0.684 | 0.391 | 0.537 | 0.623 | 0.465 | 0.500 | 0.615 | 0.476 | 0.500 | 0.645 | - |
| Q_deep vs N_1 | 101 | 72 | 0.532 | 0.606 | 0.568 | 0.522 | 0.488 | 0.548 | 0.448 | 0.627 | 0.665 | 0.669 | 0.377 | 0.489 | 0.549 | 0.487 | 0.500 | 0.548 | 0.458 | 0.500 | 0.537 | - |
| Q_near vs N_0.95 | 100 | 73 | 0.822 | 0.837 | 0.812 | 0.658 | 0.510 | 0.814 | 0.193 | 1.000 | 0.895 | 0.724 | 0.282 | 0.611 | 0.998 | 0.466 | 0.500 | 0.812 | 0.477 | 0.500 | 0.999 | - |
| Q_lone vs N_0.95 | 92 | 73 | 0.750 | 0.825 | 0.768 | 0.542 | 0.473 | 0.774 | 0.226 | 0.345 | 0.310 | 0.305 | 0.202 | 0.286 | 0.386 | 0.493 | 0.500 | 0.776 | 0.386 | 0.500 | 0.615 | - |
| Q_rep vs N_0 | 73 | 75 | 1.000 | 0.978 | 0.999 | 0.677 | 0.519 | 1.000 | 0.000 | 0.687 | 0.518 | 0.378 | 0.410 | 0.479 | 0.783 | 0.399 | 0.500 | 1.000 | 0.596 | 0.500 | 1.000 | - |
| Q_imp vs N (b <= 0.95) | 136 | 287 | 0.833 | 0.834 | 0.915 | 0.597 | 0.547 | 0.912 | 0.095 | 0.666 | 0.623 | 0.541 | 0.302 | 0.492 | 0.656 | 0.399 | 0.500 | 0.914 | 0.436 | 0.500 | 0.879 | - |
| Q (no import) vs N_imp | 440 | 154 | 0.859 | 0.713 | 0.929 | 0.601 | 0.477 | 0.947 | 0.057 | 0.678 | 0.602 | 0.545 | 0.293 | 0.534 | 0.724 | 0.459 | 0.500 | 0.957 | 0.504 | 0.500 | 0.917 | - |
| Q_ils vs N_0.95 | 74 | 73 | 0.550 | 0.592 | 0.539 | 0.454 | 0.489 | 0.477 | 0.531 | 0.635 | 0.645 | 0.646 | 0.380 | 0.511 | 0.505 | 0.473 | 0.500 | 0.480 | 0.457 | 0.500 | 0.507 | - |
| decoy: Q_deep vs decoy | 101 | 63 | 0.481 | 0.465 | 0.477 | 0.512 | 0.467 | 0.467 | 0.525 | 0.533 | 0.538 | 0.544 | 0.512 | 0.510 | 0.511 | 0.493 | 0.500 | 0.472 | 0.507 | 0.500 | 0.491 | - |
| allele genomes 0 | 192 | 172 | 0.808 | 0.759 | 0.871 | 0.546 | 0.501 | 0.875 | 0.120 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.434 | 0.500 | 0.885 | 0.428 | 0.500 | 0.875 | - |
| allele genomes 2 | 192 | 171 | 0.824 | 0.797 | 0.883 | 0.642 | 0.506 | 0.889 | 0.114 | 0.881 | 0.789 | 0.678 | 0.099 | 0.546 | 0.907 | 0.452 | 0.500 | 0.891 | 0.533 | 0.500 | 0.911 | - |
| allele genomes 8 | 192 | 170 | 0.820 | 0.807 | 0.885 | 0.608 | 0.501 | 0.897 | 0.118 | 0.901 | 0.724 | 0.548 | 0.079 | 0.486 | 0.916 | 0.444 | 0.500 | 0.895 | 0.490 | 0.500 | 0.922 | - |
| congeners 2 | 288 | 257 | 0.834 | 0.794 | 0.890 | 0.613 | 0.497 | 0.903 | 0.105 | 0.675 | 0.638 | 0.533 | 0.322 | 0.512 | 0.687 | 0.377 | 0.500 | 0.908 | 0.448 | 0.500 | 0.864 | - |
| congeners 6 | 288 | 256 | 0.801 | 0.779 | 0.871 | 0.583 | 0.511 | 0.868 | 0.133 | 0.669 | 0.592 | 0.563 | 0.294 | 0.502 | 0.683 | 0.418 | 0.500 | 0.872 | 0.513 | 0.500 | 0.853 | - |
| twin 0 | 288 | 259 | 0.921 | 0.870 | 0.924 | 0.746 | 0.505 | 0.902 | 0.100 | 0.684 | 0.617 | 0.507 | 0.285 | 0.500 | 0.696 | 0.414 | 0.500 | 0.906 | 0.516 | 0.500 | 0.879 | - |
| twin 1 | 288 | 254 | 0.876 | 0.768 | 0.828 | 0.533 | 0.522 | 0.863 | 0.137 | 0.664 | 0.612 | 0.583 | 0.312 | 0.515 | 0.673 | 0.453 | 0.500 | 0.873 | 0.444 | 0.500 | 0.831 | - |
| regime low | 192 | 168 | 0.792 | 0.771 | 0.858 | 0.610 | 0.504 | 0.858 | 0.146 | 0.684 | 0.630 | 0.595 | 0.326 | 0.528 | 0.678 | 0.479 | 0.500 | 0.872 | 0.464 | 0.500 | 0.811 | - |
| regime high | 192 | 166 | 0.857 | 0.870 | 0.958 | 0.588 | 0.495 | 0.942 | 0.069 | 0.706 | 0.656 | 0.496 | 0.308 | 0.490 | 0.695 | 0.420 | 0.500 | 0.941 | 0.445 | 0.500 | 0.907 | - |
| regime none | 192 | 179 | 0.810 | 0.728 | 0.843 | 0.571 | 0.508 | 0.871 | 0.129 | 0.632 | 0.561 | 0.544 | 0.290 | 0.516 | 0.685 | 0.434 | 0.500 | 0.872 | 0.538 | 0.500 | 0.872 | - |
| depth 1 | 159 | 135 | 0.809 | 0.783 | 0.884 | 0.581 | 0.505 | 0.888 | 0.117 | 0.684 | 0.621 | 0.559 | 0.319 | 0.511 | 0.698 | 0.465 | 0.500 | 0.892 | 0.472 | 0.500 | 0.851 | - |
| depth 2 | 136 | 122 | 0.814 | 0.783 | 0.879 | 0.591 | 0.499 | 0.880 | 0.119 | 0.638 | 0.568 | 0.505 | 0.259 | 0.452 | 0.649 | 0.409 | 0.500 | 0.890 | 0.490 | 0.500 | 0.841 | - |
| depth 4 | 138 | 127 | 0.827 | 0.818 | 0.893 | 0.559 | 0.489 | 0.892 | 0.113 | 0.679 | 0.636 | 0.563 | 0.324 | 0.528 | 0.693 | 0.424 | 0.500 | 0.896 | 0.516 | 0.500 | 0.872 | - |
| depth 0.5 | 143 | 129 | 0.819 | 0.758 | 0.864 | 0.623 | 0.522 | 0.875 | 0.132 | 0.681 | 0.624 | 0.555 | 0.324 | 0.531 | 0.691 | 0.476 | 0.500 | 0.877 | 0.451 | 0.500 | 0.867 | - |
| all, within identity bands of 0.005 | 576 | 513 | 0.594 | 0.562 | 0.724 | 0.658 | 0.726 | 0.750 | 0.242 | 0.625 | 0.630 | 0.581 | 0.328 | 0.513 | 0.583 | 0.494 | 0.500 | 0.766 | 0.454 | 0.500 | 0.727 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 44 | 34 | 0.780 | 0.912 | 1.000 | 0.805 | 0.772 | 0.805 | 0.992 | - |
| Q_deep vs N_0 | congeners 6 | 57 | 41 | 0.828 | 0.970 | 0.996 | 0.696 | 0.696 | 0.738 | 0.994 | - |
| Q_deep vs N_0 | no twin | 48 | 38 | 1.000 | 1.000 | 1.000 | 0.810 | 0.740 | 0.822 | 1.000 | - |
| Q_deep vs N_0 | twin | 53 | 37 | 0.961 | 0.878 | 0.996 | 0.683 | 0.715 | 0.702 | 0.969 | - |
| Q_deep vs N_0 | no allele genomes | 27 | 31 | 0.881 | 0.947 | 1.000 | 0.500 | 0.500 | 0.500 | 1.000 | - |
| Q_deep vs N_0 | allele genomes | 74 | 44 | 0.781 | 0.946 | 0.998 | 0.892 | 0.864 | 0.946 | 0.993 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 11 | 7 | 0.870 | 0.987 | 1.000 | 0.636 | 0.636 | 1.000 | 1.000 | - |
| Q_deep vs N_0 | recombination none | 35 | 25 | 0.842 | 0.943 | 1.000 | 0.688 | 0.683 | 0.755 | 0.990 | - |
| Q_deep vs N_0 | recombination high | 29 | 23 | 0.846 | 0.964 | 1.000 | 0.787 | 0.739 | 0.787 | 1.000 | - |
| Q_deep vs N_0.5 | congeners 2 | 44 | 41 | 0.763 | 0.945 | 0.966 | 0.781 | 0.720 | 0.802 | 0.935 | - |
| Q_deep vs N_0.5 | congeners 6 | 57 | 33 | 0.792 | 0.963 | 0.961 | 0.628 | 0.647 | 0.713 | 0.902 | - |
| Q_deep vs N_0.5 | no twin | 48 | 38 | 0.996 | 0.995 | 0.965 | 0.776 | 0.680 | 0.807 | 0.941 | - |
| Q_deep vs N_0.5 | twin | 53 | 36 | 0.976 | 0.901 | 0.963 | 0.654 | 0.690 | 0.685 | 0.857 | - |
| Q_deep vs N_0.5 | no allele genomes | 27 | 24 | 0.829 | 0.917 | 0.917 | 0.500 | 0.500 | 0.500 | 0.917 | - |
| Q_deep vs N_0.5 | allele genomes | 74 | 50 | 0.765 | 0.971 | 0.979 | 0.836 | 0.815 | 0.940 | 0.953 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 11 | 6 | 0.773 | 0.939 | 0.970 | 0.742 | 0.727 | 1.000 | 1.000 | - |
| Q_deep vs N_0.5 | recombination none | 35 | 22 | 0.822 | 0.934 | 0.971 | 0.719 | 0.714 | 0.774 | 0.961 | - |
| Q_deep vs N_0.5 | recombination high | 29 | 25 | 0.766 | 1.000 | 0.993 | 0.775 | 0.691 | 0.786 | 0.989 | - |
| Q_deep vs N_0.8 | congeners 2 | 44 | 33 | 0.654 | 0.785 | 0.767 | 0.740 | 0.692 | 0.738 | 0.784 | - |
| Q_deep vs N_0.8 | congeners 6 | 57 | 32 | 0.636 | 0.780 | 0.778 | 0.662 | 0.700 | 0.678 | 0.735 | - |
| Q_deep vs N_0.8 | no twin | 48 | 30 | 0.793 | 0.828 | 0.753 | 0.728 | 0.684 | 0.736 | 0.795 | - |
| Q_deep vs N_0.8 | twin | 53 | 35 | 0.794 | 0.736 | 0.777 | 0.660 | 0.691 | 0.648 | 0.697 | - |
| Q_deep vs N_0.8 | no allele genomes | 27 | 23 | 0.718 | 0.731 | 0.734 | 0.500 | 0.500 | 0.500 | 0.734 | - |
| Q_deep vs N_0.8 | allele genomes | 74 | 42 | 0.639 | 0.812 | 0.782 | 0.785 | 0.821 | 0.803 | 0.812 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 11 | 9 | 0.636 | 0.747 | 0.737 | 0.697 | 0.727 | 0.808 | 0.818 | - |
| Q_deep vs N_0.8 | recombination none | 35 | 23 | 0.624 | 0.673 | 0.729 | 0.640 | 0.645 | 0.655 | 0.713 | - |
| Q_deep vs N_0.8 | recombination high | 29 | 17 | 0.617 | 0.996 | 0.939 | 0.736 | 0.738 | 0.773 | 0.858 | - |
| Q_deep vs N_0.95 | congeners 2 | 44 | 38 | 0.575 | 0.638 | 0.636 | 0.709 | 0.682 | 0.698 | 0.743 | - |
| Q_deep vs N_0.95 | congeners 6 | 57 | 35 | 0.608 | 0.635 | 0.591 | 0.623 | 0.685 | 0.597 | 0.601 | - |
| Q_deep vs N_0.95 | no twin | 48 | 37 | 0.655 | 0.680 | 0.637 | 0.692 | 0.683 | 0.656 | 0.680 | - |
| Q_deep vs N_0.95 | twin | 53 | 36 | 0.605 | 0.563 | 0.597 | 0.634 | 0.683 | 0.618 | 0.631 | - |
| Q_deep vs N_0.95 | no allele genomes | 27 | 23 | 0.625 | 0.663 | 0.646 | 0.500 | 0.500 | 0.500 | 0.646 | - |
| Q_deep vs N_0.95 | allele genomes | 74 | 50 | 0.584 | 0.616 | 0.603 | 0.749 | 0.819 | 0.697 | 0.706 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 11 | 7 | 0.727 | 0.675 | 0.623 | 0.727 | 0.701 | 0.818 | 0.805 | - |
| Q_deep vs N_0.95 | recombination none | 35 | 25 | 0.625 | 0.615 | 0.638 | 0.645 | 0.649 | 0.650 | 0.687 | - |
| Q_deep vs N_0.95 | recombination high | 29 | 26 | 0.573 | 0.802 | 0.756 | 0.739 | 0.721 | 0.705 | 0.767 | - |

## L2d: the gain of each group over identity and depth (logistic models, genera held out)

AUC of grouped-CV logistic models on the target rows; gain over identity, top identity and fragments with a 95% interval over resampled genera.
| set | model | auc_base | auc | gain | gain_lo | gain_hi |
|---|---|---|---|---|---|---|
| pe | identity, depth +ancestry | 0.782 | 0.824 | 0.042 | 0.023 | 0.055 |
| pe | identity, depth +alleles | 0.782 | 0.855 | 0.073 | 0.039 | 0.116 |
| pe | identity, depth +polymorphic | 0.782 | 0.819 | 0.036 | 0.008 | 0.063 |
| pe | identity, depth +weights | 0.782 | 0.825 | 0.042 | 0.028 | 0.057 |
| pe | identity, depth +all three | 0.782 | 0.902 | 0.119 | 0.083 | 0.152 |
| pe | identity, depth +ancestry, polymorphic | 0.782 | 0.877 | 0.094 | 0.066 | 0.119 |
| pe | identity, depth +all four | 0.782 | 0.903 | 0.120 | 0.084 | 0.154 |
| pe | identity, depth +oracle fixed agreement | 0.782 | 0.783 | 0.001 | -0.005 | 0.006 |
| pe, genera with alleles | identity, depth +ancestry | 0.786 | 0.825 | 0.039 | 0.018 | 0.056 |
| pe, genera with alleles | identity, depth +alleles | 0.786 | 0.917 | 0.131 | 0.079 | 0.191 |
| pe, genera with alleles | identity, depth +polymorphic | 0.786 | 0.873 | 0.087 | 0.050 | 0.133 |
| pe, genera with alleles | identity, depth +weights | 0.786 | 0.828 | 0.042 | 0.023 | 0.058 |
| pe, genera with alleles | identity, depth +all three | 0.786 | 0.938 | 0.152 | 0.098 | 0.207 |
| pe, genera with alleles | identity, depth +ancestry, polymorphic | 0.786 | 0.907 | 0.121 | 0.082 | 0.167 |
| pe, genera with alleles | identity, depth +all four | 0.786 | 0.938 | 0.152 | 0.098 | 0.208 |
| pe, genera with alleles | identity, depth +oracle fixed agreement | 0.786 | 0.829 | 0.043 | 0.033 | 0.054 |
| se | identity, depth +ancestry | 0.767 | 0.808 | 0.042 | 0.027 | 0.056 |
| se | identity, depth +alleles | 0.767 | 0.848 | 0.081 | 0.046 | 0.124 |
| se | identity, depth +polymorphic | 0.767 | 0.815 | 0.048 | 0.024 | 0.071 |
| se | identity, depth +weights | 0.767 | 0.804 | 0.038 | 0.025 | 0.051 |
| se | identity, depth +all three | 0.767 | 0.891 | 0.125 | 0.095 | 0.161 |
| se | identity, depth +ancestry, polymorphic | 0.767 | 0.862 | 0.095 | 0.069 | 0.119 |
| se | identity, depth +all four | 0.767 | 0.891 | 0.124 | 0.094 | 0.160 |
| se | identity, depth +oracle fixed agreement | 0.767 | 0.767 | 0.000 | 0.000 | 0.000 |
| se, genera with alleles | identity, depth +ancestry | 0.766 | 0.802 | 0.036 | 0.019 | 0.054 |
| se, genera with alleles | identity, depth +alleles | 0.766 | 0.907 | 0.141 | 0.093 | 0.196 |
| se, genera with alleles | identity, depth +polymorphic | 0.766 | 0.853 | 0.087 | 0.051 | 0.127 |
| se, genera with alleles | identity, depth +weights | 0.766 | 0.799 | 0.033 | 0.019 | 0.047 |
| se, genera with alleles | identity, depth +all three | 0.766 | 0.927 | 0.160 | 0.113 | 0.214 |
| se, genera with alleles | identity, depth +ancestry, polymorphic | 0.766 | 0.890 | 0.124 | 0.090 | 0.160 |
| se, genera with alleles | identity, depth +all four | 0.766 | 0.926 | 0.160 | 0.113 | 0.215 |
| se, genera with alleles | identity, depth +oracle fixed agreement | 0.766 | 0.766 | 0.000 | 0.000 | 0.000 |
| hifi | identity, depth +ancestry | 0.808 | 0.874 | 0.066 | 0.038 | 0.094 |
| hifi | identity, depth +alleles | 0.808 | 0.885 | 0.077 | 0.037 | 0.124 |
| hifi | identity, depth +polymorphic | 0.808 | 0.894 | 0.086 | 0.042 | 0.125 |
| hifi | identity, depth +weights | 0.808 | 0.883 | 0.075 | 0.046 | 0.103 |
| hifi | identity, depth +all three | 0.808 | 0.946 | 0.138 | 0.095 | 0.180 |
| hifi | identity, depth +ancestry, polymorphic | 0.808 | 0.934 | 0.125 | 0.086 | 0.166 |
| hifi | identity, depth +all four | 0.808 | 0.952 | 0.143 | 0.104 | 0.186 |
| hifi | identity, depth +oracle fixed agreement | 0.808 | 0.808 | 0.000 | 0.000 | 0.000 |
(8 more rows in l2_increment.tsv)

## Sensitivity: by divergence of the lineage and by recombination with congeners

Per feature, the threshold 5% of all N rows exceed (the feature oriented by its AUC); rate_above: the share of Q rows above it (sensitivity) or of N rows (false-positive rate) in each bin. 'MRCA with nearest allele genome' is how recently the strain shared an ancestor with a strain the database knows; 'genes with a congener's segment' how much of it came from a congener.
| set | class | by | level | rows | identity | ancestry_agreement | allele_explained_share | polymorphic_known_share | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 271 | 0.288 | 0.314 | 0.402 | 0.303 | 0.450 | 0.000 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.078 | 0.109 | 0.062 | 0.047 | 0.109 | 0.000 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 559 | 0.351 | 0.519 | 0.000 | 0.000 | 0.193 | 0.000 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 795 | 0.397 | 0.511 | 0.701 | 0.429 | 0.642 | 0.000 |
| pe | Q | depth | 10 | 432 | 0.361 | 0.444 | 0.394 | 0.250 | 0.428 | 0.000 |
| pe | Q | depth | 100 | 426 | 0.333 | 0.437 | 0.415 | 0.270 | 0.420 | 0.000 |
| pe | Q | depth | 3 | 407 | 0.376 | 0.533 | 0.371 | 0.246 | 0.526 | 0.000 |
| pe | Q | depth | 30 | 424 | 0.340 | 0.455 | 0.406 | 0.243 | 0.399 | 0.000 |
| pe | Q | distance to rep | 0.005-0.015 | 747 | 0.368 | 0.485 | 0.503 | 0.257 | 0.355 | 0.000 |
| pe | Q | distance to rep | 0.015-0.03 | 465 | 0.037 | 0.237 | 0.488 | 0.355 | 0.359 | 0.000 |
| pe | Q | distance to rep | <0.005 | 319 | 0.950 | 0.987 | 0.053 | 0.069 | 0.887 | 0.000 |
| pe | Q | distance to rep | >=0.03 | 158 | 0.000 | 0.006 | 0.316 | 0.297 | 0.203 | 0.000 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 551 | 0.354 | 0.512 | 0.428 | 0.249 | 0.410 | 0.000 |
| pe | Q | genes with a congener's segment | <=0.1 | 453 | 0.316 | 0.426 | 0.364 | 0.194 | 0.442 | 0.000 |
| pe | Q | genes with a congener's segment | >0.3 | 253 | 0.300 | 0.486 | 0.482 | 0.285 | 0.391 | 0.000 |
| pe | Q | genes with a congener's segment | none | 432 | 0.419 | 0.440 | 0.340 | 0.299 | 0.514 | 0.000 |
| pe | Q | regime | high | 561 | 0.357 | 0.520 | 0.437 | 0.219 | 0.355 | 0.000 |
| pe | Q | regime | low | 561 | 0.289 | 0.433 | 0.394 | 0.246 | 0.460 | 0.000 |
| pe | Q | regime | none | 567 | 0.411 | 0.446 | 0.360 | 0.291 | 0.511 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 929 | 0.356 | 0.456 | 0.702 | 0.431 | 0.620 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 568 | 0.349 | 0.500 | 0.000 | 0.005 | 0.178 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 181 | 0.304 | 0.381 | 0.099 | 0.127 | 0.326 | 0.000 |
| pe | Q | role | Q_deep | 306 | 0.131 | 0.141 | 0.405 | 0.324 | 0.232 | 0.000 |
| pe | Q | role | Q_ils | 191 | 0.084 | 0.047 | 0.361 | 0.293 | 0.079 | 0.000 |
| pe | Q | role | Q_imp | 393 | 0.305 | 0.486 | 0.491 | 0.293 | 0.455 | 0.000 |
| pe | Q | role | Q_lone | 271 | 0.402 | 0.616 | 0.118 | 0.077 | 0.303 | 0.000 |
| pe | Q | role | Q_near | 329 | 0.422 | 0.559 | 0.702 | 0.410 | 0.717 | 0.000 |
| pe | Q | role | Q_rep | 199 | 0.859 | 0.975 | 0.106 | 0.000 | 0.824 | 0.000 |
| pe | Q | twin genus | 0 | 864 | 0.233 | 0.485 | 0.414 | 0.228 | 0.514 | 0.000 |
| pe | Q | twin genus | 1 | 825 | 0.478 | 0.447 | 0.378 | 0.278 | 0.367 | 0.000 |
| pe | N | depth | 10 | 354 | 0.068 | 0.056 | 0.056 | 0.042 | 0.051 | 0.000 |
| pe | N | depth | 100 | 377 | 0.003 | 0.016 | 0.005 | 0.011 | 0.013 | 0.000 |
| pe | N | depth | 3 | 289 | 0.128 | 0.111 | 0.107 | 0.097 | 0.114 | 0.000 |
| pe | N | depth | 30 | 380 | 0.021 | 0.032 | 0.037 | 0.029 | 0.037 | 0.000 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 823 | 0.049 | 0.052 | 0.044 | 0.043 | 0.055 | 0.000 |
| pe | N | genes with a segment of S or a congener | >0.3 | 236 | 0.051 | 0.030 | 0.072 | 0.055 | 0.017 | 0.000 |
| pe | N | genes with a segment of S or a congener | none | 341 | 0.053 | 0.059 | 0.041 | 0.029 | 0.062 | 0.000 |
| pe | N | regime | high | 468 | 0.062 | 0.049 | 0.060 | 0.053 | 0.034 | 0.000 |
| pe | N | regime | low | 467 | 0.039 | 0.054 | 0.047 | 0.034 | 0.064 | 0.000 |
| pe | N | regime | none | 465 | 0.049 | 0.047 | 0.037 | 0.037 | 0.052 | 0.000 |
| pe | N | role | N_0 | 172 | 0.006 | 0.012 | 0.029 | 0.029 | 0.012 | 0.000 |
| pe | N | role | N_0.5 | 192 | 0.073 | 0.026 | 0.010 | 0.021 | 0.016 | 0.000 |
| pe | N | role | N_0.8 | 210 | 0.038 | 0.038 | 0.062 | 0.062 | 0.043 | 0.000 |
| pe | N | role | N_0.95 | 214 | 0.075 | 0.089 | 0.065 | 0.033 | 0.075 | 0.000 |
| pe | N | role | N_1 | 221 | 0.077 | 0.113 | 0.059 | 0.041 | 0.145 | 0.000 |
| pe | N | role | N_imp0.2 | 194 | 0.010 | 0.026 | 0.031 | 0.046 | 0.026 | 0.000 |
| pe | N | role | N_imp0.4 | 197 | 0.061 | 0.030 | 0.071 | 0.056 | 0.015 | 0.000 |
| pe | N | twin genus | 0 | 745 | 0.000 | 0.030 | 0.031 | 0.039 | 0.038 | 0.000 |
| pe | N | twin genus | 1 | 655 | 0.107 | 0.073 | 0.067 | 0.044 | 0.064 | 0.000 |
| se | Q | MRCA with nearest allele genome | 0.002-0.01 | 255 | 0.204 | 0.278 | 0.365 | 0.263 | 0.318 |  |
| se | Q | MRCA with nearest allele genome | long ago (>=0.01) | 60 | 0.067 | 0.083 | 0.100 | 0.050 | 0.133 |  |
| se | Q | MRCA with nearest allele genome | no allele genome | 534 | 0.238 | 0.395 | 0.000 | 0.000 | 0.187 |  |
| se | Q | MRCA with nearest allele genome | recent (<0.002) | 762 | 0.256 | 0.427 | 0.672 | 0.381 | 0.521 |  |
| se | Q | depth | 10 | 415 | 0.243 | 0.366 | 0.359 | 0.224 | 0.371 |  |
| se | Q | depth | 100 | 426 | 0.214 | 0.354 | 0.415 | 0.251 | 0.277 |  |
| se | Q | depth | 3 | 347 | 0.256 | 0.455 | 0.323 | 0.210 | 0.539 |  |
| se | Q | depth | 30 | 423 | 0.229 | 0.357 | 0.409 | 0.206 | 0.300 |  |
| se | Q | distance to rep | 0.005-0.015 | 700 | 0.174 | 0.357 | 0.464 | 0.214 | 0.311 |  |
| se | Q | distance to rep | 0.015-0.03 | 449 | 0.024 | 0.140 | 0.483 | 0.321 | 0.241 |  |
| se | Q | distance to rep | <0.005 | 310 | 0.790 | 0.965 | 0.055 | 0.061 | 0.748 |  |
| se | Q | distance to rep | >=0.03 | 152 | 0.000 | 0.000 | 0.342 | 0.309 | 0.184 |  |
| se | Q | genes with a congener's segment | 0.1-0.3 | 526 | 0.196 | 0.373 | 0.413 | 0.232 | 0.325 |  |
| se | Q | genes with a congener's segment | <=0.1 | 433 | 0.206 | 0.344 | 0.351 | 0.178 | 0.323 |  |
| se | Q | genes with a congener's segment | >0.3 | 244 | 0.209 | 0.406 | 0.439 | 0.221 | 0.352 |  |
| se | Q | genes with a congener's segment | none | 408 | 0.331 | 0.412 | 0.331 | 0.262 | 0.463 |  |
| se | Q | regime | high | 538 | 0.177 | 0.392 | 0.411 | 0.191 | 0.301 |  |
| se | Q | regime | low | 540 | 0.204 | 0.341 | 0.380 | 0.215 | 0.333 |  |
| se | Q | regime | none | 533 | 0.325 | 0.407 | 0.347 | 0.265 | 0.458 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 890 | 0.224 | 0.381 | 0.663 | 0.381 | 0.483 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 540 | 0.230 | 0.378 | 0.000 | 0.004 | 0.172 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | partly | 170 | 0.259 | 0.341 | 0.124 | 0.112 | 0.306 |  |
| se | Q | role | Q_deep | 285 | 0.056 | 0.063 | 0.382 | 0.305 | 0.168 |  |
| se | Q | role | Q_ils | 182 | 0.044 | 0.049 | 0.346 | 0.242 | 0.066 |  |
| se | Q | role | Q_imp | 377 | 0.225 | 0.395 | 0.459 | 0.255 | 0.382 |  |
| se | Q | role | Q_lone | 256 | 0.246 | 0.449 | 0.113 | 0.074 | 0.246 |  |
| se | Q | role | Q_near | 314 | 0.242 | 0.462 | 0.704 | 0.363 | 0.567 |  |
| se | Q | role | Q_rep | 197 | 0.660 | 0.893 | 0.081 | 0.000 | 0.716 |  |
| se | Q | twin genus | 0 | 854 | 0.145 | 0.361 | 0.417 | 0.219 | 0.410 |  |
| se | Q | twin genus | 1 | 757 | 0.336 | 0.402 | 0.337 | 0.229 | 0.312 |  |
| se | N | depth | 10 | 322 | 0.043 | 0.053 | 0.050 | 0.078 | 0.043 |  |
| se | N | depth | 100 | 377 | 0.003 | 0.003 | 0.032 | 0.019 | 0.000 |  |
| se | N | depth | 3 | 228 | 0.048 | 0.162 | 0.088 | 0.083 | 0.193 |  |
| se | N | depth | 30 | 376 | 0.019 | 0.013 | 0.040 | 0.040 | 0.021 |  |
| se | N | genes with a segment of S or a congener | <=0.3 | 770 | 0.029 | 0.055 | 0.044 | 0.048 | 0.062 |  |
| se | N | genes with a segment of S or a congener | >0.3 | 215 | 0.023 | 0.033 | 0.060 | 0.060 | 0.033 |  |
| se | N | genes with a segment of S or a congener | none | 318 | 0.019 | 0.035 | 0.050 | 0.050 | 0.035 |  |
| se | N | regime | high | 432 | 0.023 | 0.053 | 0.035 | 0.049 | 0.053 |  |
| se | N | regime | low | 439 | 0.039 | 0.048 | 0.064 | 0.052 | 0.062 |  |
| se | N | regime | none | 432 | 0.014 | 0.037 | 0.046 | 0.051 | 0.037 |  |
| se | N | role | N_0 | 157 | 0.006 | 0.025 | 0.019 | 0.070 | 0.032 |  |
| se | N | role | N_0.5 | 180 | 0.033 | 0.028 | 0.017 | 0.033 | 0.033 |  |
| se | N | role | N_0.8 | 192 | 0.021 | 0.036 | 0.052 | 0.057 | 0.052 |  |
| se | N | role | N_0.95 | 206 | 0.044 | 0.073 | 0.049 | 0.034 | 0.063 |  |
| se | N | role | N_1 | 207 | 0.034 | 0.092 | 0.092 | 0.058 | 0.097 |  |
| se | N | role | N_imp0.2 | 182 | 0.011 | 0.022 | 0.033 | 0.038 | 0.033 |  |
| se | N | role | N_imp0.4 | 179 | 0.022 | 0.034 | 0.067 | 0.067 | 0.034 |  |
| se | N | twin genus | 0 | 720 | 0.000 | 0.018 | 0.036 | 0.049 | 0.032 |  |
| se | N | twin genus | 1 | 583 | 0.057 | 0.081 | 0.063 | 0.053 | 0.074 |  |
| hifi | Q | MRCA with nearest allele genome | 0.002-0.01 | 97 | 0.340 | 0.412 | 0.546 | 0.423 | 0.505 |  |
| hifi | Q | MRCA with nearest allele genome | long ago (>=0.01) | 20 | 0.200 | 0.200 | 0.000 | 0.200 | 0.200 |  |
| hifi | Q | MRCA with nearest allele genome | no allele genome | 192 | 0.401 | 0.583 | 0.000 | 0.000 | 0.120 |  |
| hifi | Q | MRCA with nearest allele genome | recent (<0.002) | 267 | 0.566 | 0.640 | 0.963 | 0.536 | 0.794 |  |
| hifi | Q | depth | 0.5 | 143 | 0.476 | 0.559 | 0.545 | 0.329 | 0.483 |  |
| hifi | Q | depth | 1 | 159 | 0.453 | 0.553 | 0.541 | 0.340 | 0.491 |  |
| hifi | Q | depth | 2 | 136 | 0.493 | 0.588 | 0.529 | 0.309 | 0.493 |  |
| hifi | Q | depth | 4 | 138 | 0.420 | 0.572 | 0.536 | 0.326 | 0.536 |  |
| hifi | Q | distance to rep | 0.005-0.015 | 270 | 0.607 | 0.574 | 0.607 | 0.348 | 0.507 |  |
| hifi | Q | distance to rep | 0.015-0.03 | 158 | 0.000 | 0.449 | 0.487 | 0.424 | 0.323 |  |
| hifi | Q | distance to rep | <0.005 | 101 | 1.000 | 1.000 | 0.515 | 0.069 | 0.921 |  |
| hifi | Q | distance to rep | >=0.03 | 47 | 0.000 | 0.000 | 0.362 | 0.426 | 0.149 |  |
| hifi | Q | genes with a congener's segment | 0.1-0.3 | 189 | 0.466 | 0.624 | 0.608 | 0.360 | 0.540 |  |
| hifi | Q | genes with a congener's segment | <=0.1 | 156 | 0.442 | 0.474 | 0.506 | 0.295 | 0.494 |  |
| hifi | Q | genes with a congener's segment | >0.3 | 83 | 0.434 | 0.639 | 0.518 | 0.265 | 0.217 |  |
| hifi | Q | genes with a congener's segment | none | 148 | 0.486 | 0.554 | 0.493 | 0.351 | 0.615 |  |
| hifi | Q | regime | high | 192 | 0.505 | 0.677 | 0.583 | 0.297 | 0.443 |  |
| hifi | Q | regime | low | 192 | 0.380 | 0.469 | 0.526 | 0.344 | 0.469 |  |
| hifi | Q | regime | none | 192 | 0.495 | 0.557 | 0.505 | 0.339 | 0.589 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 322 | 0.512 | 0.562 | 0.925 | 0.540 | 0.776 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 197 | 0.416 | 0.594 | 0.000 | 0.000 | 0.112 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | partly | 53 | 0.264 | 0.472 | 0.226 | 0.264 | 0.226 |  |
| hifi | Q | role | Q_deep | 101 | 0.198 | 0.317 | 0.455 | 0.475 | 0.386 |  |
| hifi | Q | role | Q_ils | 74 | 0.122 | 0.176 | 0.486 | 0.459 | 0.000 |  |
| hifi | Q | role | Q_imp | 136 | 0.397 | 0.618 | 0.581 | 0.309 | 0.375 |  |
| hifi | Q | role | Q_lone | 92 | 0.511 | 0.663 | 0.207 | 0.141 | 0.359 |  |
| hifi | Q | role | Q_near | 100 | 0.640 | 0.650 | 1.000 | 0.510 | 0.990 |  |
| hifi | Q | role | Q_rep | 73 | 0.973 | 0.986 | 0.411 | 0.000 | 0.904 |  |
| hifi | Q | twin genus | 0 | 288 | 0.316 | 0.646 | 0.524 | 0.281 | 0.483 |  |
| hifi | Q | twin genus | 1 | 288 | 0.604 | 0.490 | 0.552 | 0.372 | 0.517 |  |
| hifi | N | depth | 0.5 | 129 | 0.047 | 0.062 | 0.047 | 0.054 | 0.039 |  |
| hifi | N | depth | 1 | 135 | 0.037 | 0.037 | 0.022 | 0.067 | 0.052 |  |
| hifi | N | depth | 2 | 122 | 0.066 | 0.066 | 0.090 | 0.033 | 0.066 |  |
| hifi | N | depth | 4 | 127 | 0.055 | 0.039 | 0.047 | 0.047 | 0.047 |  |
| hifi | N | genes with a segment of S or a congener | <=0.3 | 292 | 0.041 | 0.045 | 0.051 | 0.031 | 0.048 |  |
| hifi | N | genes with a segment of S or a congener | >0.3 | 93 | 0.000 | 0.000 | 0.022 | 0.065 | 0.000 |  |
| hifi | N | genes with a segment of S or a congener | none | 128 | 0.109 | 0.102 | 0.070 | 0.086 | 0.094 |  |
| hifi | N | regime | high | 166 | 0.042 | 0.042 | 0.042 | 0.054 | 0.030 |  |
| hifi | N | regime | low | 168 | 0.030 | 0.036 | 0.060 | 0.012 | 0.054 |  |
| hifi | N | regime | none | 179 | 0.078 | 0.073 | 0.050 | 0.084 | 0.067 |  |
| hifi | N | role | N_0 | 75 | 0.000 | 0.000 | 0.013 | 0.040 | 0.000 |  |
| hifi | N | role | N_0.5 | 74 | 0.000 | 0.000 | 0.000 | 0.041 | 0.000 |  |
| hifi | N | role | N_0.8 | 65 | 0.046 | 0.000 | 0.031 | 0.031 | 0.000 |  |
| hifi | N | role | N_0.95 | 73 | 0.096 | 0.123 | 0.082 | 0.041 | 0.041 |  |
| hifi | N | role | N_1 | 72 | 0.222 | 0.236 | 0.208 | 0.069 | 0.319 |  |
| hifi | N | role | N_imp0.2 | 77 | 0.000 | 0.000 | 0.000 | 0.052 | 0.000 |  |
| hifi | N | role | N_imp0.4 | 77 | 0.000 | 0.000 | 0.026 | 0.078 | 0.000 |  |
| hifi | N | twin genus | 0 | 259 | 0.008 | 0.042 | 0.012 | 0.081 | 0.027 |  |
| hifi | N | twin genus | 1 | 254 | 0.094 | 0.059 | 0.091 | 0.020 | 0.075 |  |

## L3: the presence models (trained on the train world, called on the test world)

Models by read type and feature set (base: the default set without ancestry, alleles and polymorphic; _on_noallele / _on_shuffled: the default model on the test world's runs without allele scores and with the shuffled table). Q and N of 'all test rows': present and absent taxa of every row.
| model | stratum | Q | N | auc | tpr_at_5pct_fpr | sensitivity_at_knob | fpr_at_knob | f1_at_knob |
|---|---|---|---|---|---|---|---|---|
| hifi_base | target rows | 576 | 513 | 0.888 | 0.561 | 0.715 | 0.148 |  |
| hifi_base | hard: Q_deep vs N_0.8/0.95 | 101 | 138 | 0.695 | 0.208 | 0.545 | 0.283 |  |
| hifi_base | twin genera | 288 | 254 | 0.860 | 0.444 | 0.701 | 0.185 |  |
| hifi_base | no allele genomes | 192 | 172 | 0.861 | 0.505 | 0.646 | 0.140 |  |
| hifi_base | recombination high | 192 | 166 | 0.945 | 0.682 | 0.875 | 0.157 |  |
| hifi_base | all test rows | 896 | 873 | 0.954 | - | 0.817 | 0.087 | 0.859 |
| hifi_default | target rows | 576 | 513 | 0.948 | 0.797 | 0.837 | 0.088 |  |
| hifi_default | hard: Q_deep vs N_0.8/0.95 | 101 | 138 | 0.856 | 0.653 | 0.693 | 0.174 |  |
| hifi_default | twin genera | 288 | 254 | 0.949 | 0.795 | 0.875 | 0.126 |  |
| hifi_default | no allele genomes | 192 | 172 | 0.858 | 0.479 | 0.682 | 0.134 |  |
| hifi_default | recombination high | 192 | 166 | 0.962 | 0.875 | 0.896 | 0.084 |  |
| hifi_default | all test rows | 896 | 873 | 0.979 | - | 0.895 | 0.052 | 0.920 |
| pe_alleles | target rows | 1689 | 1400 | 0.912 | 0.708 | 0.827 | 0.156 |  |
| pe_alleles | hard: Q_deep vs N_0.8/0.95 | 306 | 424 | 0.773 | 0.454 | 0.663 | 0.229 |  |
| pe_alleles | twin genera | 825 | 655 | 0.904 | 0.695 | 0.842 | 0.205 |  |
| pe_alleles | no allele genomes | 559 | 472 | 0.819 | 0.462 | 0.689 | 0.220 |  |
| pe_alleles | recombination high | 561 | 468 | 0.931 | 0.742 | 0.868 | 0.188 |  |
| pe_alleles | all test rows | 2649 | 3376 | 0.972 | - | 0.889 | 0.069 | 0.900 |
| pe_ancestry+alleles | target rows | 1689 | 1400 | 0.909 | 0.713 | 0.808 | 0.139 |  |
| pe_ancestry+alleles | hard: Q_deep vs N_0.8/0.95 | 306 | 424 | 0.763 | 0.448 | 0.641 | 0.196 |  |
| pe_ancestry+alleles | twin genera | 825 | 655 | 0.898 | 0.703 | 0.822 | 0.182 |  |
| pe_ancestry+alleles | no allele genomes | 559 | 472 | 0.824 | 0.478 | 0.637 | 0.167 |  |
| pe_ancestry+alleles | recombination high | 561 | 468 | 0.930 | 0.745 | 0.850 | 0.188 |  |
| pe_ancestry+alleles | all test rows | 2649 | 3376 | 0.972 | - | 0.877 | 0.059 | 0.898 |
| pe_ancestry | target rows | 1689 | 1400 | 0.849 | 0.499 | 0.717 | 0.196 |  |
| pe_ancestry | hard: Q_deep vs N_0.8/0.95 | 306 | 424 | 0.649 | 0.193 | 0.539 | 0.323 |  |
| pe_ancestry | twin genera | 825 | 655 | 0.793 | 0.436 | 0.659 | 0.249 |  |
| pe_ancestry | no allele genomes | 559 | 472 | 0.824 | 0.474 | 0.671 | 0.216 |  |
| pe_ancestry | recombination high | 561 | 468 | 0.879 | 0.535 | 0.775 | 0.216 |  |
| pe_ancestry | all test rows | 2649 | 3376 | 0.956 | - | 0.819 | 0.084 | 0.850 |
| pe_base | target rows | 1689 | 1400 | 0.849 | 0.504 | 0.727 | 0.196 |  |
| pe_base | hard: Q_deep vs N_0.8/0.95 | 306 | 424 | 0.659 | 0.203 | 0.552 | 0.316 |  |
| pe_base | twin genera | 825 | 655 | 0.794 | 0.440 | 0.665 | 0.241 |  |
| pe_base | no allele genomes | 559 | 472 | 0.827 | 0.481 | 0.682 | 0.193 |  |
| pe_base | recombination high | 561 | 468 | 0.883 | 0.554 | 0.793 | 0.207 |  |
| pe_base | all test rows | 2649 | 3376 | 0.954 | - | 0.826 | 0.087 | 0.853 |
| pe_default | target rows | 1689 | 1400 | 0.914 | 0.718 | 0.831 | 0.156 |  |
| pe_default | hard: Q_deep vs N_0.8/0.95 | 306 | 424 | 0.778 | 0.464 | 0.673 | 0.222 |  |
| pe_default | twin genera | 825 | 655 | 0.908 | 0.710 | 0.859 | 0.218 |  |
| pe_default | no allele genomes | 559 | 472 | 0.825 | 0.476 | 0.696 | 0.235 |  |
| pe_default | recombination high | 561 | 468 | 0.932 | 0.765 | 0.875 | 0.190 |  |
| pe_default | all test rows | 2649 | 3376 | 0.974 | - | 0.892 | 0.065 | 0.903 |
| pe_default_on_noallele | target rows | 1719 | 1476 | 0.918 | 0.717 | 0.823 | 0.138 |  |
| pe_default_on_noallele | hard: Q_deep vs N_0.8/0.95 | 311 | 439 | 0.784 | 0.453 | 0.659 | 0.205 |  |
| pe_default_on_noallele | twin genera | 855 | 724 | 0.921 | 0.720 | 0.857 | 0.174 |  |
| pe_default_on_noallele | no allele genomes | 573 | 498 | 0.826 | 0.482 | 0.682 | 0.211 |  |
| pe_default_on_noallele | recombination high | 574 | 487 | 0.931 | 0.733 | 0.852 | 0.162 |  |
| pe_default_on_noallele | all test rows | 2679 | 3858 | 0.976 | - | 0.885 | 0.054 | 0.902 |
| pe_default_on_shuffled | target rows | 1685 | 1396 | 0.834 | 0.471 | 0.602 | 0.138 |  |
| pe_default_on_shuffled | hard: Q_deep vs N_0.8/0.95 | 306 | 424 | 0.621 | 0.127 | 0.340 | 0.196 |  |
| pe_default_on_shuffled | twin genera | 821 | 650 | 0.791 | 0.445 | 0.616 | 0.211 |  |
| pe_default_on_shuffled | no allele genomes | 559 | 472 | 0.825 | 0.469 | 0.694 | 0.239 |  |
| pe_default_on_shuffled | recombination high | 561 | 464 | 0.862 | 0.538 | 0.668 | 0.151 |  |
| pe_default_on_shuffled | all test rows | 2645 | 3380 | 0.951 | - | 0.746 | 0.057 | 0.820 |
| pe_weights | target rows | 1689 | 1400 | 0.849 | 0.506 | 0.721 | 0.201 |  |
| pe_weights | hard: Q_deep vs N_0.8/0.95 | 306 | 424 | 0.654 | 0.193 | 0.539 | 0.328 |  |
| pe_weights | twin genera | 825 | 655 | 0.792 | 0.435 | 0.657 | 0.247 |  |
| pe_weights | no allele genomes | 559 | 472 | 0.823 | 0.478 | 0.678 | 0.208 |  |
| pe_weights | recombination high | 561 | 468 | 0.883 | 0.537 | 0.799 | 0.231 |  |
| pe_weights | all test rows | 2649 | 3376 | 0.955 | - | 0.822 | 0.087 | 0.851 |
| pe_without weights | target rows | 1689 | 1400 | 0.913 | 0.713 | 0.827 | 0.157 |  |
| pe_without weights | hard: Q_deep vs N_0.8/0.95 | 306 | 424 | 0.771 | 0.454 | 0.660 | 0.226 |  |
| pe_without weights | twin genera | 825 | 655 | 0.906 | 0.696 | 0.852 | 0.215 |  |
| pe_without weights | no allele genomes | 559 | 472 | 0.826 | 0.478 | 0.692 | 0.225 |  |
| pe_without weights | recombination high | 561 | 468 | 0.931 | 0.747 | 0.872 | 0.192 |  |
| pe_without weights | all test rows | 2649 | 3376 | 0.974 | - | 0.889 | 0.068 | 0.900 |

The models' calls on the target rows by bin: 'called' the share called at the model's knob (Q: sensitivity; N: false-positive rate), 'at 5% FPR' the share above the score 5% of the model's N rows exceed.
| set | class | by | level | rows | base called | ancestry called | alleles called | default called | base at 5% FPR | ancestry at 5% FPR | alleles at 5% FPR | default at 5% FPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 271 | 0.661 | 0.661 | 0.768 | 0.786 | 0.446 | 0.454 | 0.605 | 0.631 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.391 | 0.406 | 0.375 | 0.375 | 0.172 | 0.141 | 0.156 | 0.125 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 559 | 0.682 | 0.671 | 0.689 | 0.696 | 0.447 | 0.453 | 0.490 | 0.504 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 795 | 0.809 | 0.794 | 0.980 | 0.979 | 0.591 | 0.576 | 0.940 | 0.945 |
| pe | Q | depth | 10 | 432 | 0.734 | 0.722 | 0.836 | 0.838 | 0.498 | 0.491 | 0.718 | 0.720 |
| pe | Q | depth | 100 | 426 | 0.737 | 0.732 | 0.850 | 0.862 | 0.545 | 0.540 | 0.756 | 0.763 |
| pe | Q | depth | 3 | 407 | 0.698 | 0.690 | 0.779 | 0.786 | 0.450 | 0.459 | 0.629 | 0.654 |
| pe | Q | depth | 30 | 424 | 0.738 | 0.722 | 0.840 | 0.837 | 0.524 | 0.505 | 0.724 | 0.731 |
| pe | Q | distance to rep | 0.005-0.015 | 747 | 0.731 | 0.716 | 0.885 | 0.897 | 0.458 | 0.459 | 0.731 | 0.743 |
| pe | Q | distance to rep | 0.015-0.03 | 465 | 0.671 | 0.658 | 0.783 | 0.772 | 0.368 | 0.355 | 0.624 | 0.634 |
| pe | Q | distance to rep | <0.005 | 319 | 1.000 | 1.000 | 0.994 | 1.000 | 0.975 | 0.975 | 0.978 | 0.987 |
| pe | Q | distance to rep | >=0.03 | 158 | 0.323 | 0.323 | 0.342 | 0.354 | 0.177 | 0.152 | 0.297 | 0.297 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 551 | 0.793 | 0.777 | 0.884 | 0.891 | 0.554 | 0.534 | 0.771 | 0.780 |
| pe | Q | genes with a congener's segment | <=0.1 | 453 | 0.649 | 0.649 | 0.695 | 0.695 | 0.470 | 0.492 | 0.598 | 0.609 |
| pe | Q | genes with a congener's segment | >0.3 | 253 | 0.759 | 0.743 | 0.870 | 0.877 | 0.494 | 0.498 | 0.731 | 0.759 |
| pe | Q | genes with a congener's segment | none | 432 | 0.706 | 0.697 | 0.866 | 0.870 | 0.484 | 0.463 | 0.727 | 0.727 |
| pe | Q | regime | high | 561 | 0.793 | 0.775 | 0.868 | 0.875 | 0.542 | 0.528 | 0.747 | 0.765 |
| pe | Q | regime | low | 561 | 0.668 | 0.672 | 0.734 | 0.733 | 0.496 | 0.506 | 0.636 | 0.649 |
| pe | Q | regime | none | 567 | 0.720 | 0.704 | 0.877 | 0.885 | 0.476 | 0.464 | 0.739 | 0.739 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 929 | 0.767 | 0.758 | 0.962 | 0.963 | 0.550 | 0.540 | 0.909 | 0.919 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 568 | 0.683 | 0.673 | 0.685 | 0.695 | 0.438 | 0.442 | 0.474 | 0.489 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 181 | 0.641 | 0.630 | 0.564 | 0.569 | 0.448 | 0.436 | 0.392 | 0.381 |
| pe | Q | role | Q_deep | 306 | 0.552 | 0.539 | 0.663 | 0.673 | 0.255 | 0.222 | 0.467 | 0.467 |
| pe | Q | role | Q_ils | 191 | 0.450 | 0.440 | 0.592 | 0.613 | 0.136 | 0.152 | 0.361 | 0.372 |
| pe | Q | role | Q_imp | 393 | 0.771 | 0.751 | 0.878 | 0.885 | 0.524 | 0.517 | 0.768 | 0.784 |
| pe | Q | role | Q_lone | 271 | 0.756 | 0.756 | 0.790 | 0.782 | 0.539 | 0.546 | 0.646 | 0.657 |
| pe | Q | role | Q_near | 329 | 0.815 | 0.805 | 0.994 | 0.988 | 0.623 | 0.620 | 0.960 | 0.967 |
| pe | Q | role | Q_rep | 199 | 0.990 | 0.990 | 0.975 | 0.985 | 0.960 | 0.960 | 0.955 | 0.975 |
| pe | Q | twin genus | 0 | 864 | 0.786 | 0.772 | 0.811 | 0.804 | 0.562 | 0.554 | 0.703 | 0.704 |
| pe | Q | twin genus | 1 | 825 | 0.665 | 0.659 | 0.842 | 0.859 | 0.444 | 0.441 | 0.713 | 0.732 |
| pe | N | depth | 10 | 354 | 0.203 | 0.198 | 0.175 | 0.153 | 0.051 | 0.059 | 0.062 | 0.068 |
| pe | N | depth | 100 | 377 | 0.196 | 0.202 | 0.138 | 0.151 | 0.056 | 0.053 | 0.042 | 0.050 |
| pe | N | depth | 3 | 289 | 0.208 | 0.208 | 0.159 | 0.159 | 0.048 | 0.038 | 0.062 | 0.059 |
| pe | N | depth | 30 | 380 | 0.182 | 0.179 | 0.153 | 0.163 | 0.045 | 0.047 | 0.037 | 0.026 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 823 | 0.204 | 0.210 | 0.154 | 0.163 | 0.064 | 0.062 | 0.045 | 0.051 |
| pe | N | genes with a segment of S or a congener | >0.3 | 236 | 0.097 | 0.093 | 0.119 | 0.097 | 0.004 | 0.004 | 0.047 | 0.038 |
| pe | N | genes with a segment of S or a congener | none | 341 | 0.246 | 0.232 | 0.185 | 0.182 | 0.047 | 0.053 | 0.065 | 0.056 |
| pe | N | regime | high | 468 | 0.207 | 0.216 | 0.188 | 0.190 | 0.047 | 0.049 | 0.058 | 0.056 |
| pe | N | regime | low | 467 | 0.169 | 0.171 | 0.122 | 0.128 | 0.066 | 0.060 | 0.039 | 0.049 |
| pe | N | regime | none | 465 | 0.213 | 0.200 | 0.157 | 0.151 | 0.037 | 0.041 | 0.054 | 0.045 |
| pe | N | role | N_0 | 172 | 0.006 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 192 | 0.036 | 0.031 | 0.031 | 0.042 | 0.010 | 0.005 | 0.000 | 0.005 |
| pe | N | role | N_0.8 | 210 | 0.219 | 0.252 | 0.171 | 0.171 | 0.024 | 0.048 | 0.019 | 0.038 |
| pe | N | role | N_0.95 | 214 | 0.411 | 0.393 | 0.285 | 0.271 | 0.121 | 0.107 | 0.089 | 0.084 |
| pe | N | role | N_1 | 221 | 0.443 | 0.443 | 0.371 | 0.385 | 0.163 | 0.158 | 0.154 | 0.149 |
| pe | N | role | N_imp0.2 | 194 | 0.067 | 0.067 | 0.031 | 0.052 | 0.000 | 0.000 | 0.010 | 0.010 |
| pe | N | role | N_imp0.4 | 197 | 0.112 | 0.102 | 0.137 | 0.107 | 0.005 | 0.005 | 0.056 | 0.041 |
| pe | N | twin genus | 0 | 745 | 0.157 | 0.149 | 0.113 | 0.102 | 0.048 | 0.047 | 0.036 | 0.036 |
| pe | N | twin genus | 1 | 655 | 0.241 | 0.249 | 0.205 | 0.218 | 0.052 | 0.053 | 0.066 | 0.066 |
| hifi | Q | MRCA with nearest allele genome | 0.002-0.01 | 97 | 0.629 |  |  | 0.835 | 0.474 |  |  | 0.814 |
| hifi | Q | MRCA with nearest allele genome | long ago (>=0.01) | 20 | 0.550 |  |  | 0.200 | 0.250 |  |  | 0.200 |
| hifi | Q | MRCA with nearest allele genome | no allele genome | 192 | 0.646 |  |  | 0.682 | 0.500 |  |  | 0.583 |
| hifi | Q | MRCA with nearest allele genome | recent (<0.002) | 267 | 0.809 |  |  | 0.996 | 0.659 |  |  | 0.989 |
| hifi | Q | depth | 0.5 | 143 | 0.699 |  |  | 0.825 | 0.524 |  |  | 0.776 |
| hifi | Q | depth | 1 | 159 | 0.686 |  |  | 0.818 | 0.516 |  |  | 0.786 |
| hifi | Q | depth | 2 | 136 | 0.757 |  |  | 0.868 | 0.647 |  |  | 0.831 |
| hifi | Q | depth | 4 | 138 | 0.725 |  |  | 0.841 | 0.565 |  |  | 0.797 |
| hifi | Q | distance to rep | 0.005-0.015 | 270 | 0.785 |  |  | 0.922 | 0.607 |  |  | 0.863 |
| hifi | Q | distance to rep | 0.015-0.03 | 158 | 0.551 |  |  | 0.728 | 0.329 |  |  | 0.684 |
| hifi | Q | distance to rep | <0.005 | 101 | 1.000 |  |  | 1.000 | 1.000 |  |  | 1.000 |
| hifi | Q | distance to rep | >=0.03 | 47 | 0.255 |  |  | 0.362 | 0.128 |  |  | 0.362 |
| hifi | Q | genes with a congener's segment | 0.1-0.3 | 189 | 0.815 |  |  | 0.899 | 0.640 |  |  | 0.852 |
| hifi | Q | genes with a congener's segment | <=0.1 | 156 | 0.558 |  |  | 0.692 | 0.455 |  |  | 0.667 |
| hifi | Q | genes with a congener's segment | >0.3 | 83 | 0.795 |  |  | 0.867 | 0.639 |  |  | 0.807 |
| hifi | Q | genes with a congener's segment | none | 148 | 0.709 |  |  | 0.892 | 0.527 |  |  | 0.858 |
| hifi | Q | regime | high | 192 | 0.875 |  |  | 0.896 | 0.714 |  |  | 0.859 |
| hifi | Q | regime | low | 192 | 0.552 |  |  | 0.708 | 0.448 |  |  | 0.682 |
| hifi | Q | regime | none | 192 | 0.719 |  |  | 0.906 | 0.521 |  |  | 0.849 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 322 | 0.764 |  |  | 0.984 | 0.612 |  |  | 0.978 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 |  |  | 1.000 | 1.000 |  |  | 1.000 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 197 | 0.655 |  |  | 0.690 | 0.513 |  |  | 0.594 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | partly | 53 | 0.623 |  |  | 0.472 | 0.396 |  |  | 0.434 |
| hifi | Q | role | Q_deep | 101 | 0.545 |  |  | 0.693 | 0.347 |  |  | 0.683 |
| hifi | Q | role | Q_ils | 74 | 0.500 |  |  | 0.703 | 0.203 |  |  | 0.595 |
| hifi | Q | role | Q_imp | 136 | 0.743 |  |  | 0.860 | 0.588 |  |  | 0.809 |
| hifi | Q | role | Q_lone | 92 | 0.739 |  |  | 0.793 | 0.587 |  |  | 0.728 |
| hifi | Q | role | Q_near | 100 | 0.780 |  |  | 1.000 | 0.670 |  |  | 0.990 |
| hifi | Q | role | Q_rep | 73 | 1.000 |  |  | 0.959 | 0.986 |  |  | 0.959 |
| hifi | Q | twin genus | 0 | 288 | 0.729 |  |  | 0.799 | 0.569 |  |  | 0.760 |
| hifi | Q | twin genus | 1 | 288 | 0.701 |  |  | 0.875 | 0.552 |  |  | 0.833 |
| hifi | N | depth | 0.5 | 129 | 0.155 |  |  | 0.070 | 0.054 |  |  | 0.039 |
| hifi | N | depth | 1 | 135 | 0.096 |  |  | 0.074 | 0.037 |  |  | 0.037 |
| hifi | N | depth | 2 | 122 | 0.197 |  |  | 0.107 | 0.066 |  |  | 0.066 |
| hifi | N | depth | 4 | 127 | 0.150 |  |  | 0.102 | 0.047 |  |  | 0.063 |
| hifi | N | genes with a segment of S or a congener | <=0.3 | 292 | 0.158 |  |  | 0.079 | 0.051 |  |  | 0.034 |
| hifi | N | genes with a segment of S or a congener | >0.3 | 93 | 0.011 |  |  | 0.011 | 0.000 |  |  | 0.000 |
| hifi | N | genes with a segment of S or a congener | none | 128 | 0.227 |  |  | 0.164 | 0.086 |  |  | 0.125 |
| hifi | N | regime | high | 166 | 0.157 |  |  | 0.084 | 0.054 |  |  | 0.042 |
| hifi | N | regime | low | 168 | 0.113 |  |  | 0.054 | 0.036 |  |  | 0.018 |
| hifi | N | regime | none | 179 | 0.173 |  |  | 0.123 | 0.061 |  |  | 0.089 |
| hifi | N | role | N_0 | 75 | 0.000 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_0.5 | 74 | 0.000 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_0.8 | 65 | 0.169 |  |  | 0.077 | 0.031 |  |  | 0.046 |
| hifi | N | role | N_0.95 | 73 | 0.384 |  |  | 0.260 | 0.137 |  |  | 0.151 |
| hifi | N | role | N_1 | 72 | 0.472 |  |  | 0.278 | 0.194 |  |  | 0.167 |
| hifi | N | role | N_imp0.2 | 77 | 0.026 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_imp0.4 | 77 | 0.013 |  |  | 0.013 | 0.000 |  |  | 0.000 |
| hifi | N | twin genus | 0 | 259 | 0.112 |  |  | 0.050 | 0.012 |  |  | 0.015 |
| hifi | N | twin genus | 1 | 254 | 0.185 |  |  | 0.126 | 0.091 |  |  | 0.087 |
