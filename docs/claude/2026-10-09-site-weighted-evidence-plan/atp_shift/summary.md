# Ancestry true-positive test: /home/falk/atp_shift

Worlds: test (seed 11), train (seed 23); 36 test genera. See docs/claude/2026-10-09-ancestry-true-positive-test/README.md for the design and the predictions.

## Verdicts

| hypothesis | prediction | outcome | verdict |
|---|---|---|---|
| H1a table | stored alleles are the allele genomes' differences (>= 0.95 exact) | min share exact 1.000 | PASS |
| H3 sites | consensus sites = stem + representative's lineage (no recombination, >= 3 congeners): >= 0.98 of them, recall >= 0.95 | share 0.864, recall 0.651 | FAIL |
| H3 fallback | with 2 congeners the nearest congener's own derived states enter (share in history clearly below 1) | mean share 0.838 | FAIL |
| H7 novel reads | the allele scores draw some of a novel species' reads onto the target (it carries the alleles' ancestral base at the representative's own sites) | 5714 of 116100 reads (0.049) | INFO |
| H7 settling | the strains' reads the allele scores move go to their species (>= 0.9) | 0.967 of 3313 moved | PASS |
| H1b/H2b oracle | protal's ancestry, allele and polymorphic features equal the oracle's (error-free reads, every row within 1e-6) | worst share within 1e-6: 1.000 | PASS |
| H1 leak | an allele genome's own reads: allele_explained_share >= 0.95, polymorphic_novel_share <= 0.02 | 0.959, 0.000 | PASS |
| H2 calibration | planted A_a: ancestry_agreement = a (+-0.03) | worst |diff| 0.011 | PASS |
| H1 calibration | planted P_f: allele_explained_share rises with f (a step near 0.5: a read is explained once it shares more than half of the allele's edits in its span) | f 0: 0.000, f 0.25: 0.071, f 0.5: 0.265, f 0.75: 0.534, f 1: 0.707 | PASS |
| H4 ancestry, clean strata | the same with >= 3 congeners and no recombination (the consensus applies) | b 0: 0.976, b 0.5: 0.864 | INFO |
| H5 fixed sites | deep strains with their representative-lineage sites covered (>= 0.8) against late novel species (b 0.8, 0.95), genera with alleles: ancestry_fixed_gain adds to ancestry_agreement (CI above 0) | AUC 0.634 -> 0.860, gain 0.226 [0.148, 0.307]; ancestry_fixed_agreement alone 0.813, the true fixed-site agreement alone 0.813 (160 Q, 291 N rows) | PASS |
| H4 ancestry | ancestry_agreement separates deep strains from novel species that left the stem early (b <= 0.5): AUC >= 0.9 | b 0: 0.965, b 0.5: 0.865 | FAIL |
| H6 confound | allele_explained_share favours strains with a stored relative (Q_near) over late twins, but not lone strains | Q_near 0.886, Q_lone 0.336 | PASS |
| H9 decoy | the decoy (a tip of the deep strains' clade labelled a species) against Q_deep: no feature separates (AUC 0.5 +- 0.1) | largest |AUC - 0.5|: uniqueness 0.543 (307 Q, 196 decoy rows) | PASS |
| H9b crown lineage | N_1, a novel lineage leaving at S's crown, against Q_deep: the ancestry features cannot tell them apart; the allele features can where Q_deep shares its side of the crown with allele genomes | ancestry_agreement 0.522, allele_explained_share 0.632 | INFO |
| H4/H5 gain | the three groups add separation over identity and depth (pe, CI above 0) | 0.129 [0.093, 0.166] | PASS |
| H5 ceiling | the true fixed-site agreement alone (the oracle's) adds separation; compare with '+all three' for the gap | 0.001 [-0.005, 0.008] | INFO |
| H5 ceiling, genera with alleles | the true fixed-site agreement (the oracle's) over identity and depth, where the fixed sites exist; compare with '+all three' | 0.059 [0.044, 0.075] | INFO |
| H8 calls pe target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.063 [0.036, 0.088] | PASS |
| H8 calls pe hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.124 [0.060, 0.177] | PASS |
| H8 calls hifi target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.060 [0.033, 0.086] | PASS |
| H8 calls hifi hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.161 [0.085, 0.224] | PASS |
| H9 shuffled | shuffled alleles cost the default model separation (the alleles carry information) | AUC 0.917 -> 0.835 | PASS |

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
| Q | 0 | 66509 | 0.939 | 0.945 | 487 | 453 | 0.930 |
| Q | 1 | 65698 | 0.889 | 0.929 | 2826 | 2752 | 0.974 |
| N | 0 | 57611 | 0.830 | 0.857 | 1794 | 1642 | 0.915 |
| N | 1 | 58489 | 0.757 | 0.824 | 4270 | 4072 | 0.954 |
| class | twin | from | to | reads |
|---|---|---|---|---|
| D | 0 | congener | congener | 2 |
| D | 0 | congener | target | 15 |
| D | 0 | sister | other | 1 |
| D | 0 | sister | target | 91 |
| D | 0 | target | congener | 2 |
| D | 0 | target | sister | 8 |
| D | 1 | congener | target | 16 |
| D | 1 | sister | target | 371 |
| D | 1 | target | sister | 14 |
| N | 0 | congener | congener | 6 |
| N | 0 | congener | sister | 2 |
| N | 0 | congener | target | 166 |
| N | 0 | other | other | 4 |
| N | 0 | other | target | 9 |
| N | 0 | sister | congener | 10 |
| N | 0 | sister | other | 1 |
| N | 0 | sister | target | 1467 |
| N | 0 | target | congener | 17 |
| N | 0 | target | other | 2 |
| N | 0 | target | sister | 110 |
| N | 1 | congener | sister | 5 |
| N | 1 | congener | target | 77 |
| N | 1 | other | other | 1 |
| N | 1 | other | target | 1 |
| N | 1 | sister | congener | 2 |
| N | 1 | sister | target | 3994 |
| N | 1 | target | congener | 4 |
| N | 1 | target | sister | 186 |
| Q | 0 | congener | target | 70 |
| Q | 0 | other | other | 3 |
| Q | 0 | other | target | 8 |
| Q | 0 | sister | congener | 5 |
| Q | 0 | sister | other | 1 |
| Q | 0 | sister | target | 375 |
| Q | 0 | target | congener | 8 |
| Q | 0 | target | other | 3 |
| Q | 0 | target | sister | 14 |
| Q | 1 | congener | congener | 4 |
| Q | 1 | congener | target | 93 |
| Q | 1 | sister | target | 2659 |
(2 more rows in l1_moves.tsv)
The protal log's 'strain alleles:' lines (every sample of the default run): 343414 unsure reads took the shifts, 24616 of them to another species.

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
| default | pe | ancestry_sites_per_record | 3320 | 0.000 | 1.000 |
| default | pe | ancestry_agreement | 3320 | 0.000 | 1.000 |
| default | pe | ancestry_congener_share | 3320 | 0.000 | 1.000 |
| default | pe | ancestry_indel_sites_per_record | 3320 | 0.000 | 1.000 |
| default | pe | ancestry_indel_congener_share | 3320 | 0.000 | 1.000 |
| default | pe | allele_explained_share | 3320 | 0.000 | 1.000 |
| default | pe | allele_identity_gain | 3320 | 0.000 | 1.000 |
| default | pe | polymorphic_known_share | 3320 | 0.000 | 1.000 |
| default | pe | polymorphic_novel_share | 3320 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_gain | 3320 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_agreement | 3320 | 0.000 | 1.000 |
| default | pe | conserved_mismatch_ratio | 3320 | 0.000 | 1.000 |
| default | pe | conserved_mismatch_rate | 3320 | 0.000 | 1.000 |
| default | pe | ancestry_agreement_weighted | 3320 | 0.000 | 1.000 |
| default | pe | nonsynonymous_share | 3320 | 0.000 | 1.000 |
| default | pe | nonsynonymous_conserved_rate | 3320 | 0.000 | 1.000 |
| default | pe | allele_copy_share | 3320 | 0.000 | 1.000 |
| default | pe | allele_sites_per_kb | 3320 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_share | 3320 | 0.000 | 1.000 |
| default | pe | column_weight_coverage | 3320 | 0.000 | 1.000 |
(40 more rows in l2_oracle.tsv)

## L2b: the planted genomes and the leaky control (error-free reads)

Mean features of the planted genomes (P_f: a share f of the farthest allele genome's differences; A_a: the congeners' base at a share 1 - a of the consensus sites) and of Q_leak (an allele genome itself). truth_agreement: the genome's agreement at the consensus sites by position.
| role | rows | truth_agreement | ancestry_agreement | allele_explained_share | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | identity |
|---|---|---|---|---|---|---|---|---|
| A_a0 | 7 | 0.002 | 0.003 | 0.128 | 0.298 | 0.054 | 0.000 | 0.978 |
| A_a0.25 | 7 | 0.253 | 0.261 | 0.083 | 0.222 | 0.036 | 0.012 | 0.979 |
| A_a0.5 | 6 | 0.502 | 0.498 | 0.045 | 0.148 | 0.023 | -0.006 | 0.980 |
| A_a0.75 | 7 | 0.752 | 0.749 | 0.015 | 0.073 | 0.012 | -0.003 | 0.980 |
| A_a1 | 6 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0 | 7 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0.25 | 6 | 0.873 | 0.920 | 0.071 | 0.206 | 0.000 | 0.076 | 0.980 |
| P_f0.5 | 7 | 0.925 | 0.839 | 0.265 | 0.410 | 0.000 | 0.157 | 0.980 |
| P_f0.75 | 6 | 0.817 | 0.788 | 0.534 | 0.585 | 0.000 | 0.205 | 0.980 |
| P_f1 | 6 | 0.513 | 0.732 | 0.707 | 0.797 | 0.000 | 0.266 | 0.978 |
| Q_leak | 110 | 0.855 | 0.753 | 0.959 | 0.521 | 0.000 | 0.222 | 0.989 |

## L2c: separation, pe (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1707 | 1417 | 0.771 | 0.685 | 0.789 | 0.507 | 0.495 | 0.836 | 0.175 | 0.626 | 0.565 | 0.548 | 0.391 | 0.506 | 0.652 | 0.667 | 0.500 | 0.838 | 0.505 | 0.500 | 0.846 | 0.651 |
| Q_deep vs N_0 | 307 | 175 | 0.817 | 0.811 | 0.876 | 0.409 | 0.562 | 0.965 | 0.046 | 0.693 | 0.649 | 0.658 | 0.422 | 0.674 | 0.758 | 0.588 | 0.500 | 0.970 | 0.521 | 0.500 | 0.959 | 0.755 |
| Q_deep vs N_0.5 | 307 | 194 | 0.689 | 0.684 | 0.766 | 0.527 | 0.531 | 0.865 | 0.156 | 0.692 | 0.658 | 0.659 | 0.451 | 0.638 | 0.717 | 0.575 | 0.500 | 0.874 | 0.489 | 0.500 | 0.886 | 0.715 |
| Q_deep vs N_0.8 | 307 | 212 | 0.602 | 0.612 | 0.665 | 0.516 | 0.511 | 0.681 | 0.335 | 0.643 | 0.618 | 0.614 | 0.441 | 0.563 | 0.653 | 0.527 | 0.500 | 0.689 | 0.524 | 0.500 | 0.720 | 0.652 |
| Q_deep vs N_0.95 | 307 | 217 | 0.545 | 0.555 | 0.551 | 0.522 | 0.498 | 0.583 | 0.428 | 0.644 | 0.629 | 0.640 | 0.459 | 0.561 | 0.610 | 0.558 | 0.500 | 0.587 | 0.498 | 0.500 | 0.642 | 0.610 |
| Q_deep vs N_1 | 307 | 222 | 0.507 | 0.522 | 0.518 | 0.520 | 0.508 | 0.522 | 0.491 | 0.632 | 0.622 | 0.630 | 0.459 | 0.543 | 0.574 | 0.511 | 0.500 | 0.526 | 0.489 | 0.500 | 0.568 | 0.573 |
| Q_near vs N_0.95 | 332 | 217 | 0.759 | 0.669 | 0.724 | 0.544 | 0.478 | 0.782 | 0.223 | 0.886 | 0.772 | 0.736 | 0.411 | 0.611 | 0.891 | 0.682 | 0.500 | 0.780 | 0.499 | 0.500 | 0.937 | 0.891 |
| Q_lone vs N_0.95 | 273 | 217 | 0.768 | 0.679 | 0.781 | 0.488 | 0.458 | 0.798 | 0.206 | 0.336 | 0.314 | 0.307 | 0.332 | 0.303 | 0.353 | 0.705 | 0.500 | 0.796 | 0.497 | 0.500 | 0.660 | 0.353 |
| Q_rep vs N_0 | 200 | 175 | 0.996 | 0.921 | 0.989 | 0.432 | 0.536 | 1.000 | 0.001 | 0.468 | 0.371 | 0.319 | 0.366 | 0.388 | 0.694 | 0.886 | 0.500 | 0.999 | 0.483 | 0.500 | 0.993 | 0.690 |
| Q_imp vs N (b <= 0.95) | 398 | 798 | 0.796 | 0.732 | 0.828 | 0.509 | 0.501 | 0.874 | 0.144 | 0.673 | 0.600 | 0.572 | 0.371 | 0.517 | 0.664 | 0.687 | 0.500 | 0.871 | 0.527 | 0.500 | 0.880 | 0.663 |
| Q (no import) vs N_imp | 1309 | 397 | 0.797 | 0.648 | 0.806 | 0.514 | 0.481 | 0.885 | 0.127 | 0.615 | 0.551 | 0.549 | 0.391 | 0.514 | 0.670 | 0.674 | 0.500 | 0.889 | 0.507 | 0.500 | 0.894 | 0.670 |
| Q_ils vs N_0.95 | 197 | 217 | 0.508 | 0.520 | 0.502 | 0.498 | 0.489 | 0.467 | 0.534 | 0.613 | 0.607 | 0.606 | 0.457 | 0.531 | 0.502 | 0.529 | 0.500 | 0.474 | 0.476 | 0.500 | 0.521 | 0.503 |
| decoy: Q_deep vs decoy | 307 | 196 | 0.497 | 0.478 | 0.501 | 0.543 | 0.511 | 0.491 | 0.512 | 0.522 | 0.521 | 0.518 | 0.500 | 0.512 | 0.529 | 0.485 | 0.500 | 0.495 | 0.502 | 0.500 | 0.515 | 0.529 |
| allele genomes 0 | 561 | 472 | 0.765 | 0.685 | 0.793 | 0.483 | 0.512 | 0.846 | 0.162 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.671 | 0.500 | 0.849 | 0.504 | 0.500 | 0.846 | 0.500 |
| allele genomes 2 | 573 | 471 | 0.774 | 0.685 | 0.785 | 0.515 | 0.489 | 0.831 | 0.183 | 0.786 | 0.669 | 0.639 | 0.315 | 0.515 | 0.827 | 0.683 | 0.500 | 0.832 | 0.527 | 0.500 | 0.888 | 0.827 |
| allele genomes 8 | 573 | 474 | 0.775 | 0.685 | 0.788 | 0.527 | 0.482 | 0.841 | 0.170 | 0.765 | 0.610 | 0.578 | 0.293 | 0.504 | 0.834 | 0.651 | 0.500 | 0.840 | 0.487 | 0.500 | 0.892 | 0.835 |
| congeners 2 | 854 | 694 | 0.780 | 0.675 | 0.791 | 0.507 | 0.494 | 0.838 | 0.177 | 0.636 | 0.580 | 0.543 | 0.382 | 0.506 | 0.643 | 0.656 | 0.500 | 0.841 | 0.503 | 0.500 | 0.848 | 0.642 |
| congeners 6 | 853 | 723 | 0.768 | 0.694 | 0.787 | 0.508 | 0.496 | 0.836 | 0.170 | 0.617 | 0.551 | 0.553 | 0.400 | 0.504 | 0.662 | 0.690 | 0.500 | 0.838 | 0.508 | 0.500 | 0.845 | 0.660 |
| twin 0 | 864 | 744 | 0.876 | 0.766 | 0.866 | 0.602 | 0.486 | 0.873 | 0.133 | 0.636 | 0.551 | 0.513 | 0.327 | 0.492 | 0.688 | 0.662 | 0.500 | 0.877 | 0.503 | 0.500 | 0.865 | 0.687 |
| twin 1 | 843 | 673 | 0.745 | 0.602 | 0.698 | 0.460 | 0.530 | 0.794 | 0.222 | 0.621 | 0.579 | 0.579 | 0.446 | 0.519 | 0.625 | 0.684 | 0.500 | 0.794 | 0.508 | 0.500 | 0.825 | 0.624 |
| regime low | 566 | 472 | 0.754 | 0.683 | 0.787 | 0.500 | 0.498 | 0.829 | 0.183 | 0.644 | 0.571 | 0.567 | 0.392 | 0.519 | 0.656 | 0.672 | 0.500 | 0.835 | 0.500 | 0.500 | 0.811 | 0.656 |
| regime high | 571 | 474 | 0.791 | 0.712 | 0.820 | 0.525 | 0.493 | 0.857 | 0.163 | 0.649 | 0.567 | 0.511 | 0.378 | 0.479 | 0.650 | 0.655 | 0.500 | 0.860 | 0.502 | 0.500 | 0.859 | 0.648 |
| regime none | 570 | 471 | 0.771 | 0.662 | 0.763 | 0.496 | 0.498 | 0.836 | 0.169 | 0.589 | 0.554 | 0.562 | 0.403 | 0.521 | 0.652 | 0.677 | 0.500 | 0.833 | 0.513 | 0.500 | 0.876 | 0.650 |
| depth 3 | 424 | 300 | 0.742 | 0.750 | 0.748 | 0.527 | 0.454 | 0.782 | 0.235 | 0.611 | 0.569 | 0.547 | 0.415 | 0.502 | 0.637 | 0.618 | 0.500 | 0.785 | 0.515 | 0.500 | 0.805 | 0.634 |
| depth 10 | 433 | 359 | 0.764 | 0.745 | 0.784 | 0.484 | 0.496 | 0.844 | 0.171 | 0.624 | 0.571 | 0.549 | 0.409 | 0.527 | 0.650 | 0.675 | 0.500 | 0.844 | 0.482 | 0.500 | 0.853 | 0.649 |
| depth 30 | 424 | 381 | 0.794 | 0.694 | 0.814 | 0.509 | 0.506 | 0.863 | 0.148 | 0.633 | 0.561 | 0.544 | 0.369 | 0.496 | 0.654 | 0.682 | 0.500 | 0.865 | 0.517 | 0.500 | 0.858 | 0.654 |
| depth 100 | 426 | 377 | 0.786 | 0.635 | 0.817 | 0.500 | 0.522 | 0.857 | 0.147 | 0.646 | 0.564 | 0.558 | 0.355 | 0.504 | 0.679 | 0.711 | 0.500 | 0.860 | 0.515 | 0.500 | 0.867 | 0.678 |
| all, within identity bands of 0.005 | 1707 | 1417 | 0.558 | 0.557 | 0.647 | 0.628 | 0.678 | 0.739 | 0.277 | 0.683 | 0.680 | 0.632 | 0.478 | 0.591 | 0.616 | 0.555 | 0.500 | 0.740 | 0.506 | 0.500 | 0.779 | 0.615 |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 157 | 84 | 0.783 | 0.852 | 0.977 | 0.724 | 0.669 | 0.781 | 0.958 | 0.781 |
| Q_deep vs N_0 | congeners 6 | 150 | 91 | 0.845 | 0.897 | 0.956 | 0.662 | 0.646 | 0.737 | 0.958 | 0.731 |
| Q_deep vs N_0 | no twin | 162 | 92 | 0.972 | 0.959 | 0.985 | 0.759 | 0.684 | 0.799 | 0.988 | 0.796 |
| Q_deep vs N_0 | twin | 145 | 83 | 0.916 | 0.840 | 0.948 | 0.635 | 0.627 | 0.723 | 0.929 | 0.720 |
| Q_deep vs N_0 | no allele genomes | 79 | 59 | 0.828 | 0.877 | 0.973 | 0.500 | 0.500 | 0.500 | 0.973 | 0.500 |
| Q_deep vs N_0 | allele genomes | 228 | 116 | 0.817 | 0.879 | 0.965 | 0.807 | 0.740 | 0.930 | 0.965 | 0.930 |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 20 | 0.815 | 0.823 | 0.976 | 0.768 | 0.798 | 0.945 | 0.995 | 0.946 |
| Q_deep vs N_0 | recombination none | 105 | 58 | 0.822 | 0.856 | 0.988 | 0.689 | 0.700 | 0.766 | 0.974 | 0.763 |
| Q_deep vs N_0 | recombination high | 102 | 57 | 0.808 | 0.854 | 0.940 | 0.718 | 0.643 | 0.773 | 0.963 | 0.766 |
| Q_deep vs N_0.5 | congeners 2 | 157 | 95 | 0.677 | 0.746 | 0.847 | 0.748 | 0.688 | 0.740 | 0.903 | 0.739 |
| Q_deep vs N_0.5 | congeners 6 | 150 | 99 | 0.704 | 0.785 | 0.878 | 0.637 | 0.631 | 0.694 | 0.877 | 0.692 |
| Q_deep vs N_0.5 | no twin | 162 | 100 | 0.937 | 0.926 | 0.911 | 0.729 | 0.671 | 0.802 | 0.917 | 0.801 |
| Q_deep vs N_0.5 | twin | 145 | 94 | 0.694 | 0.658 | 0.821 | 0.667 | 0.640 | 0.649 | 0.843 | 0.648 |
| Q_deep vs N_0.5 | no allele genomes | 79 | 69 | 0.707 | 0.756 | 0.855 | 0.500 | 0.500 | 0.500 | 0.855 | 0.500 |
| Q_deep vs N_0.5 | allele genomes | 228 | 125 | 0.687 | 0.774 | 0.873 | 0.804 | 0.728 | 0.845 | 0.922 | 0.844 |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 21 | 0.641 | 0.694 | 0.864 | 0.689 | 0.731 | 0.896 | 0.954 | 0.893 |
| Q_deep vs N_0.5 | recombination none | 105 | 63 | 0.719 | 0.781 | 0.886 | 0.675 | 0.684 | 0.747 | 0.920 | 0.742 |
| Q_deep vs N_0.5 | recombination high | 102 | 67 | 0.705 | 0.789 | 0.891 | 0.724 | 0.649 | 0.716 | 0.937 | 0.716 |
| Q_deep vs N_0.8 | congeners 2 | 157 | 108 | 0.609 | 0.680 | 0.679 | 0.685 | 0.622 | 0.689 | 0.755 | 0.688 |
| Q_deep vs N_0.8 | congeners 6 | 150 | 104 | 0.599 | 0.650 | 0.686 | 0.598 | 0.606 | 0.615 | 0.688 | 0.614 |
| Q_deep vs N_0.8 | no twin | 162 | 111 | 0.714 | 0.740 | 0.686 | 0.674 | 0.638 | 0.719 | 0.753 | 0.719 |
| Q_deep vs N_0.8 | twin | 145 | 101 | 0.632 | 0.605 | 0.678 | 0.619 | 0.587 | 0.606 | 0.684 | 0.606 |
| Q_deep vs N_0.8 | no allele genomes | 79 | 66 | 0.616 | 0.663 | 0.690 | 0.500 | 0.500 | 0.500 | 0.690 | 0.500 |
| Q_deep vs N_0.8 | allele genomes | 228 | 146 | 0.601 | 0.669 | 0.687 | 0.735 | 0.678 | 0.742 | 0.776 | 0.743 |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 24 | 0.557 | 0.606 | 0.584 | 0.666 | 0.761 | 0.653 | 0.709 | 0.660 |
| Q_deep vs N_0.8 | recombination none | 105 | 72 | 0.576 | 0.593 | 0.676 | 0.647 | 0.663 | 0.647 | 0.742 | 0.648 |
| Q_deep vs N_0.8 | recombination high | 102 | 74 | 0.615 | 0.736 | 0.713 | 0.651 | 0.597 | 0.671 | 0.766 | 0.672 |
| Q_deep vs N_0.95 | congeners 2 | 157 | 105 | 0.540 | 0.528 | 0.591 | 0.689 | 0.657 | 0.618 | 0.665 | 0.618 |
| Q_deep vs N_0.95 | congeners 6 | 150 | 112 | 0.552 | 0.575 | 0.581 | 0.598 | 0.625 | 0.597 | 0.612 | 0.597 |
| Q_deep vs N_0.95 | no twin | 162 | 114 | 0.565 | 0.588 | 0.585 | 0.654 | 0.653 | 0.659 | 0.663 | 0.659 |
| Q_deep vs N_0.95 | twin | 145 | 103 | 0.566 | 0.503 | 0.585 | 0.630 | 0.620 | 0.581 | 0.626 | 0.581 |
| Q_deep vs N_0.95 | no allele genomes | 79 | 72 | 0.567 | 0.550 | 0.605 | 0.500 | 0.500 | 0.500 | 0.605 | 0.500 |
| Q_deep vs N_0.95 | allele genomes | 228 | 145 | 0.551 | 0.559 | 0.586 | 0.723 | 0.710 | 0.644 | 0.666 | 0.644 |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.559 | 0.652 | 0.579 | 0.672 | 0.788 | 0.591 | 0.583 | 0.582 |
| Q_deep vs N_0.95 | recombination none | 105 | 70 | 0.560 | 0.559 | 0.609 | 0.657 | 0.691 | 0.614 | 0.650 | 0.613 |
| Q_deep vs N_0.95 | recombination high | 102 | 72 | 0.558 | 0.599 | 0.625 | 0.657 | 0.633 | 0.644 | 0.713 | 0.645 |

## L2c: separation, se (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1651 | 1334 | 0.750 | 0.685 | 0.742 | 0.491 | 0.507 | 0.807 | 0.208 | 0.625 | 0.569 | 0.546 | 0.409 | 0.498 | 0.634 | 0.629 | 0.500 | 0.808 | 0.506 | 0.500 | 0.820 | - |
| Q_deep vs N_0 | 297 | 161 | 0.771 | 0.781 | 0.799 | 0.428 | 0.592 | 0.926 | 0.084 | 0.704 | 0.666 | 0.651 | 0.457 | 0.640 | 0.703 | 0.576 | 0.500 | 0.935 | 0.519 | 0.500 | 0.919 | - |
| Q_deep vs N_0.5 | 297 | 182 | 0.661 | 0.669 | 0.707 | 0.508 | 0.539 | 0.817 | 0.210 | 0.692 | 0.660 | 0.652 | 0.453 | 0.610 | 0.702 | 0.563 | 0.500 | 0.826 | 0.471 | 0.500 | 0.847 | - |
| Q_deep vs N_0.8 | 297 | 196 | 0.578 | 0.581 | 0.613 | 0.521 | 0.496 | 0.639 | 0.368 | 0.648 | 0.633 | 0.619 | 0.445 | 0.576 | 0.644 | 0.559 | 0.500 | 0.647 | 0.505 | 0.500 | 0.703 | - |
| Q_deep vs N_0.95 | 297 | 212 | 0.514 | 0.527 | 0.518 | 0.522 | 0.505 | 0.537 | 0.463 | 0.662 | 0.647 | 0.640 | 0.482 | 0.566 | 0.615 | 0.526 | 0.500 | 0.540 | 0.502 | 0.500 | 0.619 | - |
| Q_deep vs N_1 | 297 | 212 | 0.513 | 0.506 | 0.493 | 0.533 | 0.508 | 0.482 | 0.522 | 0.640 | 0.625 | 0.617 | 0.473 | 0.541 | 0.559 | 0.520 | 0.500 | 0.481 | 0.492 | 0.500 | 0.542 | - |
| Q_near vs N_0.95 | 320 | 212 | 0.703 | 0.652 | 0.672 | 0.525 | 0.494 | 0.743 | 0.269 | 0.877 | 0.775 | 0.730 | 0.435 | 0.603 | 0.863 | 0.621 | 0.500 | 0.743 | 0.511 | 0.500 | 0.904 | - |
| Q_lone vs N_0.95 | 260 | 212 | 0.740 | 0.667 | 0.752 | 0.486 | 0.470 | 0.779 | 0.228 | 0.371 | 0.346 | 0.328 | 0.366 | 0.308 | 0.348 | 0.630 | 0.500 | 0.778 | 0.515 | 0.500 | 0.640 | - |
| Q_rep vs N_0 | 197 | 161 | 0.979 | 0.931 | 0.969 | 0.429 | 0.569 | 0.984 | 0.016 | 0.480 | 0.404 | 0.344 | 0.389 | 0.387 | 0.642 | 0.818 | 0.500 | 0.984 | 0.486 | 0.500 | 0.973 | - |
| Q_imp vs N (b <= 0.95) | 388 | 751 | 0.773 | 0.714 | 0.778 | 0.468 | 0.512 | 0.846 | 0.175 | 0.667 | 0.605 | 0.568 | 0.388 | 0.504 | 0.649 | 0.631 | 0.500 | 0.847 | 0.524 | 0.500 | 0.853 | - |
| Q (no import) vs N_imp | 1263 | 371 | 0.772 | 0.672 | 0.748 | 0.494 | 0.504 | 0.851 | 0.164 | 0.612 | 0.554 | 0.546 | 0.416 | 0.516 | 0.654 | 0.640 | 0.500 | 0.854 | 0.509 | 0.500 | 0.876 | - |
| Q_ils vs N_0.95 | 189 | 212 | 0.521 | 0.511 | 0.489 | 0.498 | 0.478 | 0.454 | 0.548 | 0.610 | 0.596 | 0.595 | 0.451 | 0.531 | 0.525 | 0.486 | 0.500 | 0.456 | 0.491 | 0.500 | 0.500 | - |
| decoy: Q_deep vs decoy | 297 | 189 | 0.484 | 0.478 | 0.510 | 0.546 | 0.503 | 0.488 | 0.511 | 0.516 | 0.522 | 0.505 | 0.529 | 0.517 | 0.528 | 0.530 | 0.500 | 0.487 | 0.508 | 0.500 | 0.521 | - |
| allele genomes 0 | 535 | 433 | 0.747 | 0.685 | 0.758 | 0.468 | 0.517 | 0.828 | 0.181 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.619 | 0.500 | 0.831 | 0.513 | 0.500 | 0.828 | - |
| allele genomes 2 | 560 | 452 | 0.755 | 0.679 | 0.731 | 0.484 | 0.511 | 0.795 | 0.226 | 0.760 | 0.650 | 0.618 | 0.350 | 0.494 | 0.787 | 0.650 | 0.500 | 0.797 | 0.525 | 0.500 | 0.845 | - |
| allele genomes 8 | 556 | 449 | 0.750 | 0.690 | 0.739 | 0.522 | 0.494 | 0.808 | 0.206 | 0.765 | 0.630 | 0.590 | 0.347 | 0.506 | 0.792 | 0.619 | 0.500 | 0.805 | 0.481 | 0.500 | 0.861 | - |
| congeners 2 | 828 | 654 | 0.760 | 0.682 | 0.741 | 0.490 | 0.512 | 0.806 | 0.211 | 0.639 | 0.583 | 0.544 | 0.404 | 0.495 | 0.627 | 0.624 | 0.500 | 0.808 | 0.506 | 0.500 | 0.815 | - |
| congeners 6 | 823 | 680 | 0.746 | 0.687 | 0.744 | 0.491 | 0.504 | 0.809 | 0.202 | 0.611 | 0.554 | 0.546 | 0.413 | 0.500 | 0.641 | 0.639 | 0.500 | 0.811 | 0.506 | 0.500 | 0.826 | - |
| twin 0 | 857 | 723 | 0.854 | 0.750 | 0.827 | 0.533 | 0.510 | 0.848 | 0.160 | 0.628 | 0.554 | 0.512 | 0.347 | 0.485 | 0.672 | 0.620 | 0.500 | 0.852 | 0.514 | 0.500 | 0.847 | - |
| twin 1 | 794 | 611 | 0.688 | 0.608 | 0.645 | 0.488 | 0.549 | 0.758 | 0.261 | 0.624 | 0.587 | 0.576 | 0.468 | 0.514 | 0.604 | 0.645 | 0.500 | 0.758 | 0.498 | 0.500 | 0.788 | - |
| regime low | 552 | 454 | 0.734 | 0.691 | 0.743 | 0.473 | 0.517 | 0.802 | 0.212 | 0.639 | 0.576 | 0.558 | 0.419 | 0.504 | 0.637 | 0.633 | 0.500 | 0.804 | 0.498 | 0.500 | 0.790 | - |
| regime high | 551 | 439 | 0.765 | 0.691 | 0.762 | 0.512 | 0.489 | 0.818 | 0.203 | 0.654 | 0.575 | 0.524 | 0.388 | 0.481 | 0.633 | 0.617 | 0.500 | 0.819 | 0.503 | 0.500 | 0.823 | - |
| regime none | 548 | 441 | 0.751 | 0.673 | 0.724 | 0.491 | 0.518 | 0.809 | 0.202 | 0.585 | 0.554 | 0.552 | 0.418 | 0.511 | 0.632 | 0.639 | 0.500 | 0.810 | 0.517 | 0.500 | 0.853 | - |
| depth 3 | 375 | 250 | 0.736 | 0.747 | 0.710 | 0.500 | 0.484 | 0.724 | 0.299 | 0.612 | 0.581 | 0.539 | 0.425 | 0.502 | 0.604 | 0.547 | 0.500 | 0.726 | 0.483 | 0.500 | 0.746 | - |
| depth 10 | 426 | 328 | 0.734 | 0.756 | 0.717 | 0.470 | 0.484 | 0.809 | 0.211 | 0.620 | 0.575 | 0.549 | 0.426 | 0.519 | 0.635 | 0.635 | 0.500 | 0.808 | 0.490 | 0.500 | 0.824 | - |
| depth 30 | 424 | 379 | 0.768 | 0.724 | 0.776 | 0.498 | 0.519 | 0.846 | 0.165 | 0.641 | 0.570 | 0.547 | 0.389 | 0.493 | 0.641 | 0.654 | 0.500 | 0.845 | 0.544 | 0.500 | 0.845 | - |
| depth 100 | 426 | 377 | 0.762 | 0.634 | 0.780 | 0.490 | 0.531 | 0.844 | 0.162 | 0.645 | 0.562 | 0.558 | 0.395 | 0.491 | 0.659 | 0.695 | 0.500 | 0.848 | 0.509 | 0.500 | 0.855 | - |
| all, within identity bands of 0.005 | 1651 | 1334 | 0.540 | 0.571 | 0.579 | 0.586 | 0.672 | 0.713 | 0.301 | 0.693 | 0.690 | 0.640 | 0.493 | 0.597 | 0.610 | 0.551 | 0.500 | 0.715 | 0.507 | 0.500 | 0.759 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 154 | 78 | 0.717 | 0.739 | 0.907 | 0.726 | 0.692 | 0.712 | 0.899 | - |
| Q_deep vs N_0 | congeners 6 | 143 | 83 | 0.818 | 0.855 | 0.945 | 0.673 | 0.611 | 0.693 | 0.942 | - |
| Q_deep vs N_0 | no twin | 159 | 89 | 0.940 | 0.890 | 0.952 | 0.757 | 0.681 | 0.780 | 0.951 | - |
| Q_deep vs N_0 | twin | 138 | 72 | 0.779 | 0.750 | 0.897 | 0.656 | 0.619 | 0.642 | 0.883 | - |
| Q_deep vs N_0 | no allele genomes | 77 | 51 | 0.795 | 0.815 | 0.930 | 0.500 | 0.500 | 0.500 | 0.930 | - |
| Q_deep vs N_0 | allele genomes | 220 | 110 | 0.767 | 0.796 | 0.926 | 0.815 | 0.733 | 0.831 | 0.929 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 33 | 18 | 0.670 | 0.753 | 0.956 | 0.770 | 0.705 | 0.818 | 0.985 | - |
| Q_deep vs N_0 | recombination none | 101 | 52 | 0.732 | 0.782 | 0.978 | 0.683 | 0.680 | 0.696 | 0.962 | - |
| Q_deep vs N_0 | recombination high | 99 | 52 | 0.789 | 0.783 | 0.879 | 0.747 | 0.658 | 0.747 | 0.906 | - |
| Q_deep vs N_0.5 | congeners 2 | 154 | 88 | 0.658 | 0.688 | 0.780 | 0.727 | 0.675 | 0.724 | 0.835 | - |
| Q_deep vs N_0.5 | congeners 6 | 143 | 94 | 0.675 | 0.730 | 0.857 | 0.655 | 0.629 | 0.678 | 0.864 | - |
| Q_deep vs N_0.5 | no twin | 159 | 98 | 0.840 | 0.858 | 0.865 | 0.729 | 0.663 | 0.767 | 0.884 | - |
| Q_deep vs N_0.5 | twin | 138 | 84 | 0.617 | 0.607 | 0.767 | 0.662 | 0.637 | 0.667 | 0.794 | - |
| Q_deep vs N_0.5 | no allele genomes | 77 | 62 | 0.677 | 0.751 | 0.816 | 0.500 | 0.500 | 0.500 | 0.816 | - |
| Q_deep vs N_0.5 | allele genomes | 220 | 120 | 0.658 | 0.692 | 0.827 | 0.800 | 0.722 | 0.826 | 0.879 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 33 | 20 | 0.586 | 0.564 | 0.861 | 0.753 | 0.809 | 0.820 | 0.933 | - |
| Q_deep vs N_0.5 | recombination none | 101 | 62 | 0.678 | 0.709 | 0.864 | 0.688 | 0.691 | 0.722 | 0.893 | - |
| Q_deep vs N_0.5 | recombination high | 99 | 59 | 0.675 | 0.737 | 0.833 | 0.734 | 0.651 | 0.692 | 0.907 | - |
| Q_deep vs N_0.8 | congeners 2 | 154 | 98 | 0.598 | 0.617 | 0.633 | 0.688 | 0.645 | 0.682 | 0.732 | - |
| Q_deep vs N_0.8 | congeners 6 | 143 | 98 | 0.567 | 0.615 | 0.647 | 0.603 | 0.594 | 0.605 | 0.679 | - |
| Q_deep vs N_0.8 | no twin | 159 | 106 | 0.660 | 0.673 | 0.658 | 0.684 | 0.645 | 0.690 | 0.728 | - |
| Q_deep vs N_0.8 | twin | 138 | 90 | 0.541 | 0.573 | 0.619 | 0.616 | 0.592 | 0.621 | 0.676 | - |
| Q_deep vs N_0.8 | no allele genomes | 77 | 61 | 0.627 | 0.666 | 0.681 | 0.500 | 0.500 | 0.500 | 0.681 | - |
| Q_deep vs N_0.8 | allele genomes | 220 | 135 | 0.564 | 0.595 | 0.629 | 0.739 | 0.691 | 0.729 | 0.734 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 33 | 21 | 0.495 | 0.579 | 0.602 | 0.666 | 0.768 | 0.651 | 0.801 | - |
| Q_deep vs N_0.8 | recombination none | 101 | 65 | 0.557 | 0.576 | 0.672 | 0.630 | 0.661 | 0.632 | 0.767 | - |
| Q_deep vs N_0.8 | recombination high | 99 | 68 | 0.567 | 0.619 | 0.623 | 0.683 | 0.618 | 0.639 | 0.699 | - |
| Q_deep vs N_0.95 | congeners 2 | 154 | 106 | 0.502 | 0.498 | 0.542 | 0.702 | 0.664 | 0.638 | 0.661 | - |
| Q_deep vs N_0.95 | congeners 6 | 143 | 106 | 0.528 | 0.538 | 0.534 | 0.620 | 0.613 | 0.590 | 0.574 | - |
| Q_deep vs N_0.95 | no twin | 159 | 115 | 0.536 | 0.542 | 0.556 | 0.672 | 0.659 | 0.664 | 0.658 | - |
| Q_deep vs N_0.95 | twin | 138 | 97 | 0.478 | 0.481 | 0.512 | 0.653 | 0.621 | 0.574 | 0.578 | - |
| Q_deep vs N_0.95 | no allele genomes | 77 | 72 | 0.527 | 0.564 | 0.556 | 0.500 | 0.500 | 0.500 | 0.556 | - |
| Q_deep vs N_0.95 | allele genomes | 220 | 140 | 0.530 | 0.515 | 0.542 | 0.748 | 0.708 | 0.653 | 0.643 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 33 | 22 | 0.574 | 0.562 | 0.523 | 0.679 | 0.738 | 0.635 | 0.544 | - |
| Q_deep vs N_0.95 | recombination none | 101 | 67 | 0.543 | 0.502 | 0.554 | 0.655 | 0.657 | 0.633 | 0.633 | - |
| Q_deep vs N_0.95 | recombination high | 99 | 70 | 0.514 | 0.546 | 0.551 | 0.693 | 0.648 | 0.645 | 0.684 | - |

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
| pe | identity, depth +ancestry | 0.782 | 0.832 | 0.050 | 0.032 | 0.063 |
| pe | identity, depth +alleles | 0.782 | 0.855 | 0.073 | 0.038 | 0.117 |
| pe | identity, depth +polymorphic | 0.782 | 0.815 | 0.033 | 0.006 | 0.063 |
| pe | identity, depth +weights | 0.782 | 0.837 | 0.055 | 0.040 | 0.066 |
| pe | identity, depth +all three | 0.782 | 0.911 | 0.129 | 0.093 | 0.166 |
| pe | identity, depth +ancestry, polymorphic | 0.782 | 0.889 | 0.107 | 0.077 | 0.137 |
| pe | identity, depth +all four | 0.782 | 0.913 | 0.131 | 0.093 | 0.168 |
| pe | identity, depth +oracle fixed agreement | 0.782 | 0.783 | 0.001 | -0.005 | 0.008 |
| pe, genera with alleles | identity, depth +ancestry | 0.786 | 0.833 | 0.047 | 0.026 | 0.066 |
| pe, genera with alleles | identity, depth +alleles | 0.786 | 0.917 | 0.131 | 0.076 | 0.191 |
| pe, genera with alleles | identity, depth +polymorphic | 0.786 | 0.881 | 0.095 | 0.055 | 0.146 |
| pe, genera with alleles | identity, depth +weights | 0.786 | 0.838 | 0.052 | 0.034 | 0.070 |
| pe, genera with alleles | identity, depth +all three | 0.786 | 0.944 | 0.158 | 0.100 | 0.214 |
| pe, genera with alleles | identity, depth +ancestry, polymorphic | 0.786 | 0.917 | 0.131 | 0.089 | 0.183 |
| pe, genera with alleles | identity, depth +all four | 0.786 | 0.945 | 0.159 | 0.101 | 0.217 |
| pe, genera with alleles | identity, depth +oracle fixed agreement | 0.786 | 0.845 | 0.059 | 0.044 | 0.075 |
| se | identity, depth +ancestry | 0.765 | 0.810 | 0.045 | 0.028 | 0.059 |
| se | identity, depth +alleles | 0.765 | 0.848 | 0.082 | 0.047 | 0.125 |
| se | identity, depth +polymorphic | 0.765 | 0.803 | 0.038 | 0.013 | 0.062 |
| se | identity, depth +weights | 0.765 | 0.812 | 0.046 | 0.034 | 0.060 |
| se | identity, depth +all three | 0.765 | 0.900 | 0.135 | 0.102 | 0.169 |
| se | identity, depth +ancestry, polymorphic | 0.765 | 0.866 | 0.101 | 0.072 | 0.130 |
| se | identity, depth +all four | 0.765 | 0.901 | 0.136 | 0.103 | 0.167 |
| se | identity, depth +oracle fixed agreement | 0.765 | 0.765 | 0.000 | 0.000 | 0.000 |
| se, genera with alleles | identity, depth +ancestry | 0.768 | 0.804 | 0.036 | 0.019 | 0.053 |
| se, genera with alleles | identity, depth +alleles | 0.768 | 0.908 | 0.141 | 0.091 | 0.196 |
| se, genera with alleles | identity, depth +polymorphic | 0.768 | 0.844 | 0.077 | 0.039 | 0.113 |
| se, genera with alleles | identity, depth +weights | 0.768 | 0.806 | 0.039 | 0.023 | 0.053 |
| se, genera with alleles | identity, depth +all three | 0.768 | 0.932 | 0.164 | 0.114 | 0.218 |
| se, genera with alleles | identity, depth +ancestry, polymorphic | 0.768 | 0.889 | 0.121 | 0.081 | 0.160 |
| se, genera with alleles | identity, depth +all four | 0.768 | 0.932 | 0.164 | 0.115 | 0.219 |
| se, genera with alleles | identity, depth +oracle fixed agreement | 0.768 | 0.768 | 0.000 | 0.000 | 0.000 |
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
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 276 | 0.272 | 0.351 | 0.435 | 0.333 | 0.496 | 0.000 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.078 | 0.109 | 0.016 | 0.062 | 0.094 | 0.000 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 561 | 0.392 | 0.517 | 0.000 | 0.000 | 0.189 | 0.000 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 806 | 0.402 | 0.527 | 0.730 | 0.460 | 0.769 | 0.000 |
| pe | Q | depth | 10 | 433 | 0.376 | 0.469 | 0.418 | 0.268 | 0.480 | 0.000 |
| pe | Q | depth | 100 | 426 | 0.345 | 0.455 | 0.415 | 0.289 | 0.507 | 0.000 |
| pe | Q | depth | 3 | 424 | 0.392 | 0.528 | 0.410 | 0.281 | 0.547 | 0.000 |
| pe | Q | depth | 30 | 424 | 0.349 | 0.467 | 0.417 | 0.257 | 0.502 | 0.000 |
| pe | Q | distance to rep | 0.005-0.015 | 763 | 0.383 | 0.491 | 0.533 | 0.301 | 0.467 | 0.000 |
| pe | Q | distance to rep | 0.015-0.03 | 467 | 0.051 | 0.270 | 0.493 | 0.358 | 0.409 | 0.000 |
| pe | Q | distance to rep | <0.005 | 319 | 0.966 | 0.994 | 0.078 | 0.075 | 0.906 | 0.000 |
| pe | Q | distance to rep | >=0.03 | 158 | 0.000 | 0.006 | 0.297 | 0.291 | 0.209 | 0.000 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 562 | 0.368 | 0.511 | 0.457 | 0.278 | 0.496 | 0.000 |
| pe | Q | genes with a congener's segment | <=0.1 | 457 | 0.330 | 0.440 | 0.383 | 0.208 | 0.488 | 0.000 |
| pe | Q | genes with a congener's segment | >0.3 | 255 | 0.322 | 0.510 | 0.494 | 0.302 | 0.482 | 0.000 |
| pe | Q | genes with a congener's segment | none | 433 | 0.425 | 0.464 | 0.349 | 0.321 | 0.564 | 0.000 |
| pe | Q | regime | high | 571 | 0.366 | 0.529 | 0.462 | 0.240 | 0.443 | 0.000 |
| pe | Q | regime | low | 566 | 0.313 | 0.443 | 0.417 | 0.258 | 0.514 | 0.000 |
| pe | Q | regime | none | 570 | 0.418 | 0.467 | 0.367 | 0.323 | 0.570 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 945 | 0.354 | 0.469 | 0.731 | 0.466 | 0.740 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 570 | 0.389 | 0.507 | 0.000 | 0.005 | 0.174 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 181 | 0.309 | 0.420 | 0.099 | 0.133 | 0.331 | 0.000 |
| pe | Q | role | Q_deep | 307 | 0.143 | 0.173 | 0.414 | 0.352 | 0.296 | 0.000 |
| pe | Q | role | Q_ils | 197 | 0.091 | 0.066 | 0.411 | 0.320 | 0.122 | 0.000 |
| pe | Q | role | Q_imp | 398 | 0.327 | 0.503 | 0.500 | 0.324 | 0.538 | 0.000 |
| pe | Q | role | Q_lone | 273 | 0.436 | 0.612 | 0.121 | 0.092 | 0.322 | 0.000 |
| pe | Q | role | Q_near | 332 | 0.416 | 0.572 | 0.738 | 0.428 | 0.852 | 0.000 |
| pe | Q | role | Q_rep | 200 | 0.875 | 0.980 | 0.120 | 0.000 | 0.845 | 0.000 |
| pe | Q | twin genus | 0 | 864 | 0.257 | 0.507 | 0.407 | 0.222 | 0.541 | 0.000 |
| pe | Q | twin genus | 1 | 843 | 0.477 | 0.452 | 0.423 | 0.326 | 0.477 | 0.000 |
| pe | N | depth | 10 | 359 | 0.064 | 0.061 | 0.045 | 0.047 | 0.050 | 0.000 |
| pe | N | depth | 100 | 377 | 0.008 | 0.027 | 0.005 | 0.016 | 0.024 | 0.000 |
| pe | N | depth | 3 | 300 | 0.123 | 0.090 | 0.123 | 0.107 | 0.103 | 0.000 |
| pe | N | depth | 30 | 381 | 0.021 | 0.031 | 0.039 | 0.042 | 0.034 | 0.000 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 830 | 0.052 | 0.053 | 0.049 | 0.051 | 0.057 | 0.000 |
| pe | N | genes with a segment of S or a congener | >0.3 | 241 | 0.037 | 0.017 | 0.071 | 0.062 | 0.017 | 0.000 |
| pe | N | genes with a segment of S or a congener | none | 346 | 0.055 | 0.066 | 0.035 | 0.040 | 0.058 | 0.000 |
| pe | N | regime | high | 474 | 0.059 | 0.053 | 0.057 | 0.059 | 0.036 | 0.000 |
| pe | N | regime | low | 472 | 0.040 | 0.044 | 0.059 | 0.049 | 0.066 | 0.000 |
| pe | N | regime | none | 471 | 0.051 | 0.053 | 0.032 | 0.042 | 0.049 | 0.000 |
| pe | N | role | N_0 | 175 | 0.006 | 0.006 | 0.029 | 0.040 | 0.006 | 0.000 |
| pe | N | role | N_0.5 | 194 | 0.077 | 0.015 | 0.021 | 0.052 | 0.010 | 0.000 |
| pe | N | role | N_0.8 | 212 | 0.047 | 0.033 | 0.061 | 0.075 | 0.028 | 0.000 |
| pe | N | role | N_0.95 | 217 | 0.074 | 0.088 | 0.065 | 0.051 | 0.078 | 0.000 |
| pe | N | role | N_1 | 222 | 0.077 | 0.153 | 0.054 | 0.032 | 0.171 | 0.000 |
| pe | N | role | N_imp0.2 | 198 | 0.015 | 0.020 | 0.045 | 0.040 | 0.020 | 0.000 |
| pe | N | role | N_imp0.4 | 199 | 0.045 | 0.015 | 0.065 | 0.060 | 0.015 | 0.000 |
| pe | N | twin genus | 0 | 744 | 0.000 | 0.039 | 0.022 | 0.035 | 0.038 | 0.000 |
| pe | N | twin genus | 1 | 673 | 0.105 | 0.062 | 0.080 | 0.067 | 0.064 | 0.000 |
| se | Q | MRCA with nearest allele genome | 0.002-0.01 | 266 | 0.207 | 0.305 | 0.414 | 0.192 | 0.342 |  |
| se | Q | MRCA with nearest allele genome | long ago (>=0.01) | 62 | 0.081 | 0.097 | 0.032 | 0.032 | 0.129 |  |
| se | Q | MRCA with nearest allele genome | no allele genome | 535 | 0.277 | 0.477 | 0.000 | 0.000 | 0.140 |  |
| se | Q | MRCA with nearest allele genome | recent (<0.002) | 788 | 0.246 | 0.463 | 0.712 | 0.364 | 0.539 |  |
| se | Q | depth | 10 | 426 | 0.246 | 0.418 | 0.406 | 0.209 | 0.369 |  |
| se | Q | depth | 100 | 426 | 0.202 | 0.392 | 0.427 | 0.204 | 0.272 |  |
| se | Q | depth | 3 | 375 | 0.315 | 0.499 | 0.379 | 0.229 | 0.528 |  |
| se | Q | depth | 30 | 424 | 0.219 | 0.413 | 0.415 | 0.184 | 0.302 |  |
| se | Q | distance to rep | 0.005-0.015 | 728 | 0.179 | 0.423 | 0.508 | 0.240 | 0.315 |  |
| se | Q | distance to rep | 0.015-0.03 | 457 | 0.035 | 0.201 | 0.484 | 0.265 | 0.245 |  |
| se | Q | distance to rep | <0.005 | 312 | 0.821 | 0.981 | 0.109 | 0.074 | 0.728 |  |
| se | Q | distance to rep | >=0.03 | 154 | 0.000 | 0.006 | 0.312 | 0.136 | 0.201 |  |
| se | Q | genes with a congener's segment | 0.1-0.3 | 540 | 0.198 | 0.435 | 0.439 | 0.233 | 0.319 |  |
| se | Q | genes with a congener's segment | <=0.1 | 444 | 0.216 | 0.390 | 0.392 | 0.146 | 0.345 |  |
| se | Q | genes with a congener's segment | >0.3 | 250 | 0.204 | 0.464 | 0.468 | 0.232 | 0.332 |  |
| se | Q | genes with a congener's segment | none | 417 | 0.355 | 0.439 | 0.348 | 0.218 | 0.458 |  |
| se | Q | regime | high | 551 | 0.176 | 0.441 | 0.430 | 0.200 | 0.283 |  |
| se | Q | regime | low | 552 | 0.210 | 0.399 | 0.428 | 0.183 | 0.346 |  |
| se | Q | regime | none | 548 | 0.345 | 0.445 | 0.365 | 0.235 | 0.460 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 924 | 0.212 | 0.412 | 0.703 | 0.342 | 0.506 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 0.000 | 0.000 | 0.909 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 541 | 0.266 | 0.464 | 0.000 | 0.004 | 0.128 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | partly | 175 | 0.291 | 0.366 | 0.131 | 0.126 | 0.297 |  |
| se | Q | role | Q_deep | 297 | 0.077 | 0.098 | 0.411 | 0.259 | 0.192 |  |
| se | Q | role | Q_ils | 189 | 0.042 | 0.069 | 0.376 | 0.196 | 0.085 |  |
| se | Q | role | Q_imp | 388 | 0.224 | 0.466 | 0.482 | 0.260 | 0.366 |  |
| se | Q | role | Q_lone | 260 | 0.281 | 0.542 | 0.123 | 0.085 | 0.204 |  |
| se | Q | role | Q_near | 320 | 0.241 | 0.497 | 0.747 | 0.316 | 0.600 |  |
| se | Q | role | Q_rep | 197 | 0.680 | 0.934 | 0.112 | 0.010 | 0.706 |  |
| se | Q | twin genus | 0 | 857 | 0.156 | 0.435 | 0.394 | 0.152 | 0.406 |  |
| se | Q | twin genus | 1 | 794 | 0.338 | 0.421 | 0.422 | 0.264 | 0.316 |  |
| se | N | depth | 10 | 328 | 0.088 | 0.061 | 0.064 | 0.073 | 0.040 |  |
| se | N | depth | 100 | 377 | 0.000 | 0.005 | 0.024 | 0.019 | 0.003 |  |
| se | N | depth | 3 | 250 | 0.108 | 0.156 | 0.096 | 0.092 | 0.184 |  |
| se | N | depth | 30 | 379 | 0.029 | 0.016 | 0.034 | 0.034 | 0.018 |  |
| se | N | genes with a segment of S or a congener | <=0.3 | 786 | 0.053 | 0.056 | 0.042 | 0.055 | 0.061 |  |
| se | N | genes with a segment of S or a congener | >0.3 | 223 | 0.049 | 0.022 | 0.067 | 0.049 | 0.022 |  |
| se | N | genes with a segment of S or a congener | none | 325 | 0.043 | 0.055 | 0.058 | 0.040 | 0.043 |  |
| se | N | regime | high | 439 | 0.048 | 0.059 | 0.041 | 0.055 | 0.059 |  |
| se | N | regime | low | 454 | 0.057 | 0.040 | 0.062 | 0.059 | 0.051 |  |
| se | N | regime | none | 441 | 0.045 | 0.052 | 0.048 | 0.036 | 0.041 |  |
| se | N | role | N_0 | 161 | 0.025 | 0.025 | 0.037 | 0.056 | 0.031 |  |
| se | N | role | N_0.5 | 182 | 0.060 | 0.033 | 0.033 | 0.044 | 0.022 |  |
| se | N | role | N_0.8 | 196 | 0.056 | 0.046 | 0.061 | 0.046 | 0.077 |  |
| se | N | role | N_0.95 | 212 | 0.066 | 0.080 | 0.042 | 0.047 | 0.052 |  |
| se | N | role | N_1 | 212 | 0.047 | 0.113 | 0.071 | 0.066 | 0.108 |  |
| se | N | role | N_imp0.2 | 186 | 0.048 | 0.016 | 0.027 | 0.038 | 0.027 |  |
| se | N | role | N_imp0.4 | 185 | 0.043 | 0.022 | 0.076 | 0.054 | 0.022 |  |
| se | N | twin genus | 0 | 723 | 0.000 | 0.030 | 0.024 | 0.029 | 0.028 |  |
| se | N | twin genus | 1 | 611 | 0.110 | 0.074 | 0.082 | 0.075 | 0.077 |  |
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
| pe_alleles | target rows | 1707 | 1417 | 0.918 | 0.710 | 0.825 | 0.140 |  |
| pe_alleles | hard: Q_deep vs N_0.8/0.95 | 307 | 429 | 0.783 | 0.433 | 0.661 | 0.219 |  |
| pe_alleles | twin genera | 843 | 673 | 0.916 | 0.690 | 0.845 | 0.174 |  |
| pe_alleles | no allele genomes | 561 | 472 | 0.819 | 0.456 | 0.695 | 0.216 |  |
| pe_alleles | recombination high | 571 | 474 | 0.936 | 0.727 | 0.853 | 0.150 |  |
| pe_alleles | all test rows | 2667 | 3342 | 0.973 | - | 0.888 | 0.065 | 0.902 |
| pe_ancestry+alleles | target rows | 1707 | 1417 | 0.916 | 0.721 | 0.794 | 0.121 |  |
| pe_ancestry+alleles | hard: Q_deep vs N_0.8/0.95 | 307 | 429 | 0.780 | 0.459 | 0.612 | 0.182 |  |
| pe_ancestry+alleles | twin genera | 843 | 673 | 0.909 | 0.700 | 0.811 | 0.163 |  |
| pe_ancestry+alleles | no allele genomes | 561 | 472 | 0.824 | 0.472 | 0.610 | 0.167 |  |
| pe_ancestry+alleles | recombination high | 571 | 474 | 0.935 | 0.750 | 0.841 | 0.131 |  |
| pe_ancestry+alleles | all test rows | 2667 | 3342 | 0.974 | - | 0.868 | 0.053 | 0.897 |
| pe_ancestry | target rows | 1707 | 1417 | 0.854 | 0.503 | 0.711 | 0.179 |  |
| pe_ancestry | hard: Q_deep vs N_0.8/0.95 | 307 | 429 | 0.655 | 0.199 | 0.521 | 0.312 |  |
| pe_ancestry | twin genera | 843 | 673 | 0.815 | 0.472 | 0.648 | 0.210 |  |
| pe_ancestry | no allele genomes | 561 | 472 | 0.820 | 0.469 | 0.654 | 0.182 |  |
| pe_ancestry | recombination high | 571 | 474 | 0.888 | 0.559 | 0.785 | 0.188 |  |
| pe_ancestry | all test rows | 2667 | 3342 | 0.956 | - | 0.815 | 0.078 | 0.852 |
| pe_base | target rows | 1707 | 1417 | 0.854 | 0.493 | 0.707 | 0.185 |  |
| pe_base | hard: Q_deep vs N_0.8/0.95 | 307 | 429 | 0.655 | 0.186 | 0.528 | 0.324 |  |
| pe_base | twin genera | 843 | 673 | 0.817 | 0.438 | 0.651 | 0.207 |  |
| pe_base | no allele genomes | 561 | 472 | 0.819 | 0.456 | 0.636 | 0.184 |  |
| pe_base | recombination high | 571 | 474 | 0.889 | 0.557 | 0.764 | 0.184 |  |
| pe_base | all test rows | 2667 | 3342 | 0.954 | - | 0.813 | 0.084 | 0.847 |
| pe_default | target rows | 1707 | 1417 | 0.917 | 0.709 | 0.815 | 0.143 |  |
| pe_default | hard: Q_deep vs N_0.8/0.95 | 307 | 429 | 0.779 | 0.436 | 0.642 | 0.224 |  |
| pe_default | twin genera | 843 | 673 | 0.916 | 0.669 | 0.849 | 0.198 |  |
| pe_default | no allele genomes | 561 | 472 | 0.817 | 0.415 | 0.674 | 0.229 |  |
| pe_default | recombination high | 571 | 474 | 0.934 | 0.716 | 0.837 | 0.137 |  |
| pe_default | all test rows | 2667 | 3342 | 0.974 | - | 0.882 | 0.063 | 0.899 |
| pe_default_on_noallele | target rows | 1719 | 1476 | 0.916 | 0.716 | 0.808 | 0.138 |  |
| pe_default_on_noallele | hard: Q_deep vs N_0.8/0.95 | 311 | 439 | 0.776 | 0.424 | 0.617 | 0.210 |  |
| pe_default_on_noallele | twin genera | 855 | 724 | 0.915 | 0.710 | 0.840 | 0.186 |  |
| pe_default_on_noallele | no allele genomes | 573 | 498 | 0.825 | 0.452 | 0.682 | 0.223 |  |
| pe_default_on_noallele | recombination high | 574 | 487 | 0.929 | 0.709 | 0.824 | 0.144 |  |
| pe_default_on_noallele | all test rows | 2679 | 3858 | 0.977 | - | 0.877 | 0.054 | 0.897 |
| pe_default_on_shuffled | target rows | 1691 | 1399 | 0.835 | 0.458 | 0.588 | 0.127 |  |
| pe_default_on_shuffled | hard: Q_deep vs N_0.8/0.95 | 306 | 423 | 0.609 | 0.114 | 0.324 | 0.189 |  |
| pe_default_on_shuffled | twin genera | 827 | 656 | 0.800 | 0.414 | 0.609 | 0.191 |  |
| pe_default_on_shuffled | no allele genomes | 561 | 472 | 0.815 | 0.408 | 0.679 | 0.225 |  |
| pe_default_on_shuffled | recombination high | 565 | 466 | 0.867 | 0.504 | 0.655 | 0.139 |  |
| pe_default_on_shuffled | all test rows | 2651 | 3364 | 0.950 | - | 0.737 | 0.056 | 0.816 |
| pe_weights | target rows | 1707 | 1417 | 0.855 | 0.514 | 0.703 | 0.190 |  |
| pe_weights | hard: Q_deep vs N_0.8/0.95 | 307 | 429 | 0.656 | 0.186 | 0.511 | 0.343 |  |
| pe_weights | twin genera | 843 | 673 | 0.817 | 0.474 | 0.643 | 0.217 |  |
| pe_weights | no allele genomes | 561 | 472 | 0.822 | 0.476 | 0.647 | 0.184 |  |
| pe_weights | recombination high | 571 | 474 | 0.887 | 0.538 | 0.771 | 0.194 |  |
| pe_weights | all test rows | 2667 | 3342 | 0.955 | - | 0.810 | 0.086 | 0.845 |
| pe_without weights | target rows | 1707 | 1417 | 0.918 | 0.714 | 0.822 | 0.140 |  |
| pe_without weights | hard: Q_deep vs N_0.8/0.95 | 307 | 429 | 0.782 | 0.430 | 0.655 | 0.217 |  |
| pe_without weights | twin genera | 843 | 673 | 0.915 | 0.668 | 0.854 | 0.198 |  |
| pe_without weights | no allele genomes | 561 | 472 | 0.820 | 0.415 | 0.688 | 0.222 |  |
| pe_without weights | recombination high | 571 | 474 | 0.935 | 0.713 | 0.844 | 0.143 |  |
| pe_without weights | all test rows | 2667 | 3342 | 0.975 | - | 0.886 | 0.063 | 0.902 |

The models' calls on the target rows by bin: 'called' the share called at the model's knob (Q: sensitivity; N: false-positive rate), 'at 5% FPR' the share above the score 5% of the model's N rows exceed.
| set | class | by | level | rows | base called | ancestry called | alleles called | default called | base at 5% FPR | ancestry at 5% FPR | alleles at 5% FPR | default at 5% FPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 276 | 0.634 | 0.641 | 0.768 | 0.761 | 0.442 | 0.457 | 0.598 | 0.612 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.391 | 0.406 | 0.297 | 0.281 | 0.125 | 0.125 | 0.141 | 0.156 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 561 | 0.636 | 0.654 | 0.695 | 0.674 | 0.428 | 0.449 | 0.503 | 0.487 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 806 | 0.806 | 0.799 | 0.976 | 0.974 | 0.584 | 0.587 | 0.938 | 0.942 |
| pe | Q | depth | 10 | 433 | 0.704 | 0.732 | 0.829 | 0.831 | 0.497 | 0.497 | 0.727 | 0.718 |
| pe | Q | depth | 100 | 426 | 0.735 | 0.739 | 0.852 | 0.838 | 0.552 | 0.559 | 0.756 | 0.758 |
| pe | Q | depth | 3 | 424 | 0.667 | 0.665 | 0.769 | 0.752 | 0.401 | 0.425 | 0.618 | 0.627 |
| pe | Q | depth | 30 | 424 | 0.722 | 0.708 | 0.849 | 0.837 | 0.521 | 0.533 | 0.738 | 0.733 |
| pe | Q | distance to rep | 0.005-0.015 | 763 | 0.712 | 0.706 | 0.887 | 0.882 | 0.459 | 0.478 | 0.743 | 0.746 |
| pe | Q | distance to rep | 0.015-0.03 | 467 | 0.640 | 0.655 | 0.779 | 0.752 | 0.345 | 0.338 | 0.621 | 0.606 |
| pe | Q | distance to rep | <0.005 | 319 | 0.997 | 0.997 | 0.997 | 0.997 | 0.966 | 0.981 | 0.969 | 0.975 |
| pe | Q | distance to rep | >=0.03 | 158 | 0.297 | 0.323 | 0.310 | 0.310 | 0.139 | 0.146 | 0.291 | 0.304 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 562 | 0.756 | 0.772 | 0.870 | 0.861 | 0.528 | 0.541 | 0.744 | 0.751 |
| pe | Q | genes with a congener's segment | <=0.1 | 457 | 0.626 | 0.628 | 0.687 | 0.696 | 0.468 | 0.481 | 0.597 | 0.606 |
| pe | Q | genes with a congener's segment | >0.3 | 255 | 0.749 | 0.765 | 0.871 | 0.855 | 0.494 | 0.498 | 0.749 | 0.749 |
| pe | Q | genes with a congener's segment | none | 433 | 0.704 | 0.688 | 0.885 | 0.857 | 0.471 | 0.480 | 0.762 | 0.741 |
| pe | Q | regime | high | 571 | 0.764 | 0.785 | 0.853 | 0.837 | 0.525 | 0.539 | 0.732 | 0.736 |
| pe | Q | regime | low | 566 | 0.647 | 0.650 | 0.726 | 0.737 | 0.482 | 0.491 | 0.629 | 0.645 |
| pe | Q | regime | none | 570 | 0.711 | 0.698 | 0.895 | 0.870 | 0.470 | 0.479 | 0.768 | 0.747 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 945 | 0.760 | 0.757 | 0.954 | 0.953 | 0.537 | 0.545 | 0.907 | 0.913 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 570 | 0.639 | 0.656 | 0.691 | 0.670 | 0.425 | 0.446 | 0.486 | 0.468 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 181 | 0.630 | 0.630 | 0.558 | 0.536 | 0.448 | 0.436 | 0.370 | 0.387 |
| pe | Q | role | Q_deep | 307 | 0.528 | 0.521 | 0.661 | 0.642 | 0.231 | 0.241 | 0.466 | 0.479 |
| pe | Q | role | Q_ils | 197 | 0.401 | 0.411 | 0.584 | 0.589 | 0.132 | 0.137 | 0.371 | 0.355 |
| pe | Q | role | Q_imp | 398 | 0.749 | 0.759 | 0.882 | 0.872 | 0.508 | 0.513 | 0.766 | 0.769 |
| pe | Q | role | Q_lone | 273 | 0.714 | 0.736 | 0.799 | 0.766 | 0.538 | 0.542 | 0.674 | 0.645 |
| pe | Q | role | Q_near | 332 | 0.834 | 0.822 | 0.985 | 0.985 | 0.620 | 0.639 | 0.958 | 0.964 |
| pe | Q | role | Q_rep | 200 | 0.980 | 0.985 | 0.970 | 0.975 | 0.945 | 0.970 | 0.945 | 0.960 |
| pe | Q | twin genus | 0 | 864 | 0.762 | 0.773 | 0.806 | 0.781 | 0.542 | 0.545 | 0.703 | 0.691 |
| pe | Q | twin genus | 1 | 843 | 0.651 | 0.648 | 0.845 | 0.849 | 0.442 | 0.460 | 0.718 | 0.728 |
| pe | N | depth | 10 | 359 | 0.192 | 0.189 | 0.156 | 0.153 | 0.058 | 0.064 | 0.061 | 0.058 |
| pe | N | depth | 100 | 377 | 0.180 | 0.164 | 0.109 | 0.122 | 0.050 | 0.045 | 0.045 | 0.048 |
| pe | N | depth | 3 | 300 | 0.200 | 0.193 | 0.167 | 0.180 | 0.040 | 0.037 | 0.057 | 0.063 |
| pe | N | depth | 30 | 381 | 0.171 | 0.171 | 0.136 | 0.126 | 0.050 | 0.052 | 0.039 | 0.034 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 830 | 0.195 | 0.193 | 0.136 | 0.136 | 0.063 | 0.064 | 0.043 | 0.054 |
| pe | N | genes with a segment of S or a congener | >0.3 | 241 | 0.050 | 0.058 | 0.071 | 0.079 | 0.012 | 0.004 | 0.037 | 0.029 |
| pe | N | genes with a segment of S or a congener | none | 346 | 0.254 | 0.228 | 0.199 | 0.205 | 0.046 | 0.049 | 0.075 | 0.055 |
| pe | N | regime | high | 474 | 0.184 | 0.188 | 0.150 | 0.137 | 0.034 | 0.040 | 0.051 | 0.055 |
| pe | N | regime | low | 472 | 0.157 | 0.153 | 0.102 | 0.121 | 0.074 | 0.070 | 0.032 | 0.047 |
| pe | N | regime | none | 471 | 0.214 | 0.195 | 0.170 | 0.172 | 0.042 | 0.040 | 0.068 | 0.049 |
| pe | N | role | N_0 | 175 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 194 | 0.026 | 0.021 | 0.031 | 0.031 | 0.010 | 0.000 | 0.005 | 0.005 |
| pe | N | role | N_0.8 | 212 | 0.245 | 0.255 | 0.165 | 0.156 | 0.033 | 0.038 | 0.042 | 0.033 |
| pe | N | role | N_0.95 | 217 | 0.401 | 0.369 | 0.272 | 0.290 | 0.120 | 0.129 | 0.083 | 0.092 |
| pe | N | role | N_1 | 222 | 0.441 | 0.423 | 0.342 | 0.347 | 0.144 | 0.149 | 0.140 | 0.153 |
| pe | N | role | N_imp0.2 | 198 | 0.045 | 0.040 | 0.035 | 0.035 | 0.005 | 0.005 | 0.020 | 0.010 |
| pe | N | role | N_imp0.4 | 199 | 0.055 | 0.065 | 0.080 | 0.085 | 0.015 | 0.005 | 0.040 | 0.035 |
| pe | N | twin genus | 0 | 744 | 0.165 | 0.151 | 0.110 | 0.094 | 0.048 | 0.051 | 0.032 | 0.028 |
| pe | N | twin genus | 1 | 673 | 0.207 | 0.210 | 0.174 | 0.198 | 0.052 | 0.049 | 0.070 | 0.074 |
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
