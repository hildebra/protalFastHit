# Ancestry true-positive test: /home/falk/atp_rates

Worlds: test (seed 11), train (seed 23); 36 test genera. See docs/claude/2026-10-09-ancestry-true-positive-test/README.md for the design and the predictions.

## Verdicts

| hypothesis | prediction | outcome | verdict |
|---|---|---|---|
| H1a table | stored alleles are the allele genomes' differences (>= 0.95 exact) | min share exact 1.000 | PASS |
| H3 sites | consensus sites = stem + representative's lineage (no recombination, >= 3 congeners): >= 0.98 of them, recall >= 0.95 | share 0.864, recall 0.651 | FAIL |
| H3 fallback | with 2 congeners the nearest congener's own derived states enter (share in history clearly below 1) | mean share 0.838 | FAIL |
| H7 novel reads | the allele scores draw some of a novel species' reads onto the target (it carries the alleles' ancestral base at the representative's own sites) | 6388 of 116100 reads (0.055) | INFO |
| H7 settling | the strains' reads the allele scores move go to their species (>= 0.9) | 0.999 of 3264 moved | PASS |
| H1b/H2b oracle | protal's ancestry, allele and polymorphic features equal the oracle's (error-free reads, every row within 1e-6) | worst share within 1e-6: 1.000 | PASS |
| H1 leak | an allele genome's own reads: allele_explained_share >= 0.95, polymorphic_novel_share <= 0.02 | 0.956, 0.000 | PASS |
| H2 calibration | planted A_a: ancestry_agreement = a (+-0.03) | worst |diff| 0.010 | PASS |
| H1 calibration | planted P_f: allele_explained_share rises with f (a step near 0.5: a read is explained once it shares more than half of the allele's edits in its span) | f 0: 0.000, f 0.25: 0.072, f 0.5: 0.266, f 0.75: 0.533, f 1: 0.705 | PASS |
| H4 ancestry, clean strata | the same with >= 3 congeners and no recombination (the consensus applies) | b 0: 0.978, b 0.5: 0.863 | INFO |
| H5 fixed sites | deep strains with their representative-lineage sites covered (>= 0.8) against late novel species (b 0.8, 0.95), genera with alleles: ancestry_fixed_gain adds to ancestry_agreement (CI above 0) | AUC 0.636 -> 0.879, gain 0.243 [0.165, 0.328]; ancestry_fixed_agreement alone 0.837, the true fixed-site agreement alone 0.837 (161 Q, 300 N rows) | PASS |
| H4 ancestry | ancestry_agreement separates deep strains from novel species that left the stem early (b <= 0.5): AUC >= 0.9 | b 0: 0.967, b 0.5: 0.867 | FAIL |
| H6 confound | allele_explained_share favours strains with a stored relative (Q_near) over late twins, but not lone strains | Q_near 0.893, Q_lone 0.328 | PASS |
| H9 decoy | the decoy (a tip of the deep strains' clade labelled a species) against Q_deep: no feature separates (AUC 0.5 +- 0.1) | largest |AUC - 0.5|: uniqueness 0.537 (312 Q, 199 decoy rows) | PASS |
| H9b crown lineage | N_1, a novel lineage leaving at S's crown, against Q_deep: the ancestry features cannot tell them apart; the allele features can where Q_deep shares its side of the crown with allele genomes | ancestry_agreement 0.519, allele_explained_share 0.643 | INFO |
| H4/H5 gain | the three groups add separation over identity and depth (pe, CI above 0) | 0.134 [0.099, 0.171] | PASS |
| H5 ceiling | the true fixed-site agreement alone (the oracle's) adds separation; compare with '+all three' for the gap | 0.001 [-0.006, 0.009] | INFO |
| H5 ceiling, genera with alleles | the true fixed-site agreement (the oracle's) over identity and depth, where the fixed sites exist; compare with '+all three' | 0.079 [0.054, 0.104] | INFO |
| H8 calls pe target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.068 [0.042, 0.098] | PASS |
| H8 calls pe hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.130 [0.067, 0.190] | PASS |
| H8 calls hifi target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.060 [0.033, 0.086] | PASS |
| H8 calls hifi hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.161 [0.085, 0.224] | PASS |
| H9 shuffled | shuffled alleles cost the default model separation (the alleles carry information) | AUC 0.925 -> 0.837 | PASS |

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

The build's line: Strain alleles: 7538 alleles (102032 edits) of 2762 gene copies of 24 species, up to 4 each (of 38880 full-reference copies of the database's species and genes: 0 of genomes outside --allele_genome_share 1, 26074 identical to the representative's, 2353 repeated, 19 covering less than 0.5 of it or more than 0.1 apart, 386 as far as the nearest congener's copy or farther; 10048 alleles offered): /home/falk/atp_rates/test/db/strain_alleles.tsv

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
| Q | 0 | 66510 | 0.939 | 0.946 | 511 | 507 | 0.992 |
| Q | 1 | 65698 | 0.889 | 0.930 | 2753 | 2753 | 1.000 |
| N | 0 | 57611 | 0.830 | 0.866 | 2035 | 2029 | 0.997 |
| N | 1 | 58489 | 0.757 | 0.832 | 4360 | 4359 | 1.000 |
| class | twin | from | to | reads |
|---|---|---|---|---|
| D | 0 | congener | target | 20 |
| D | 0 | sister | other | 1 |
| D | 0 | sister | target | 98 |
| D | 1 | congener | target | 16 |
| D | 1 | sister | target | 385 |
| N | 0 | congener | target | 214 |
| N | 0 | other | other | 5 |
| N | 0 | other | target | 10 |
| N | 0 | sister | other | 1 |
| N | 0 | sister | target | 1805 |
| N | 1 | congener | target | 80 |
| N | 1 | other | other | 1 |
| N | 1 | other | target | 1 |
| N | 1 | sister | target | 4278 |
| Q | 0 | congener | target | 66 |
| Q | 0 | none | target | 1 |
| Q | 0 | other | other | 3 |
| Q | 0 | other | target | 8 |
| Q | 0 | sister | target | 432 |
| Q | 0 | target | other | 1 |
| Q | 1 | congener | target | 90 |
| Q | 1 | sister | target | 2663 |
The protal log's 'strain alleles:' lines (every sample of the default run): 343414 unsure reads took the shifts, 27332 of them to another species.

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
| default | pe | ancestry_sites_per_record | 3421 | 0.000 | 1.000 |
| default | pe | ancestry_agreement | 3421 | 0.000 | 1.000 |
| default | pe | ancestry_congener_share | 3421 | 0.000 | 1.000 |
| default | pe | ancestry_indel_sites_per_record | 3421 | 0.000 | 1.000 |
| default | pe | ancestry_indel_congener_share | 3421 | 0.000 | 1.000 |
| default | pe | allele_explained_share | 3421 | 0.000 | 1.000 |
| default | pe | allele_identity_gain | 3421 | 0.000 | 1.000 |
| default | pe | polymorphic_known_share | 3421 | 0.000 | 1.000 |
| default | pe | polymorphic_novel_share | 3421 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_gain | 3421 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_agreement | 3421 | 0.000 | 1.000 |
| default | pe | conserved_mismatch_ratio | 3421 | 0.000 | 1.000 |
| default | pe | conserved_mismatch_rate | 3421 | 0.000 | 1.000 |
| default | pe | ancestry_agreement_weighted | 3421 | 0.000 | 1.000 |
| default | pe | nonsynonymous_share | 3421 | 0.000 | 1.000 |
| default | pe | nonsynonymous_conserved_rate | 3421 | 0.000 | 1.000 |
| default | pe | allele_copy_share | 3421 | 0.000 | 1.000 |
| default | pe | allele_sites_per_kb | 3421 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_share | 3421 | 0.000 | 1.000 |
| default | pe | column_weight_coverage | 3421 | 0.000 | 1.000 |
(40 more rows in l2_oracle.tsv)

## L2b: the planted genomes and the leaky control (error-free reads)

Mean features of the planted genomes (P_f: a share f of the farthest allele genome's differences; A_a: the congeners' base at a share 1 - a of the consensus sites) and of Q_leak (an allele genome itself). truth_agreement: the genome's agreement at the consensus sites by position.
| role | rows | truth_agreement | ancestry_agreement | allele_explained_share | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | identity |
|---|---|---|---|---|---|---|---|---|
| A_a0 | 7 | 0.002 | 0.003 | 0.131 | 0.298 | 0.056 | 0.000 | 0.977 |
| A_a0.25 | 7 | 0.253 | 0.260 | 0.088 | 0.228 | 0.039 | 0.011 | 0.979 |
| A_a0.5 | 6 | 0.502 | 0.498 | 0.044 | 0.150 | 0.024 | -0.006 | 0.980 |
| A_a0.75 | 7 | 0.752 | 0.749 | 0.015 | 0.072 | 0.012 | -0.004 | 0.980 |
| A_a1 | 6 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0 | 7 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0.25 | 6 | 0.873 | 0.920 | 0.072 | 0.207 | 0.000 | 0.078 | 0.980 |
| P_f0.5 | 7 | 0.925 | 0.839 | 0.266 | 0.411 | 0.000 | 0.160 | 0.980 |
| P_f0.75 | 6 | 0.817 | 0.788 | 0.533 | 0.582 | 0.001 | 0.208 | 0.980 |
| P_f1 | 6 | 0.513 | 0.731 | 0.705 | 0.798 | 0.000 | 0.268 | 0.978 |
| Q_leak | 110 | 0.855 | 0.752 | 0.956 | 0.519 | 0.000 | 0.241 | 0.989 |

## L2c: separation, pe (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1725 | 1497 | 0.777 | 0.684 | 0.804 | 0.566 | 0.505 | 0.838 | 0.173 | 0.626 | 0.559 | 0.547 | 0.371 | 0.504 | 0.661 | 0.662 | 0.500 | 0.840 | 0.504 | 0.500 | 0.843 | 0.661 |
| Q_deep vs N_0 | 312 | 195 | 0.801 | 0.806 | 0.900 | 0.549 | 0.576 | 0.967 | 0.047 | 0.706 | 0.657 | 0.673 | 0.401 | 0.684 | 0.760 | 0.574 | 0.500 | 0.972 | 0.511 | 0.500 | 0.968 | 0.757 |
| Q_deep vs N_0.5 | 312 | 206 | 0.694 | 0.674 | 0.812 | 0.584 | 0.540 | 0.867 | 0.154 | 0.685 | 0.648 | 0.662 | 0.440 | 0.638 | 0.724 | 0.584 | 0.500 | 0.875 | 0.484 | 0.500 | 0.882 | 0.724 |
| Q_deep vs N_0.8 | 312 | 219 | 0.613 | 0.636 | 0.684 | 0.541 | 0.519 | 0.679 | 0.337 | 0.652 | 0.622 | 0.633 | 0.436 | 0.567 | 0.654 | 0.525 | 0.500 | 0.684 | 0.521 | 0.500 | 0.716 | 0.655 |
| Q_deep vs N_0.95 | 312 | 226 | 0.556 | 0.573 | 0.552 | 0.543 | 0.503 | 0.589 | 0.423 | 0.651 | 0.630 | 0.651 | 0.443 | 0.566 | 0.614 | 0.556 | 0.500 | 0.593 | 0.497 | 0.500 | 0.653 | 0.613 |
| Q_deep vs N_1 | 312 | 227 | 0.503 | 0.527 | 0.529 | 0.531 | 0.511 | 0.519 | 0.491 | 0.643 | 0.631 | 0.649 | 0.464 | 0.553 | 0.570 | 0.508 | 0.500 | 0.521 | 0.485 | 0.500 | 0.580 | 0.570 |
| Q_near vs N_0.95 | 334 | 226 | 0.773 | 0.669 | 0.723 | 0.585 | 0.483 | 0.786 | 0.221 | 0.893 | 0.769 | 0.732 | 0.379 | 0.610 | 0.920 | 0.676 | 0.500 | 0.784 | 0.500 | 0.500 | 0.949 | 0.920 |
| Q_lone vs N_0.95 | 277 | 226 | 0.780 | 0.672 | 0.777 | 0.524 | 0.463 | 0.806 | 0.202 | 0.328 | 0.304 | 0.302 | 0.296 | 0.303 | 0.353 | 0.702 | 0.500 | 0.804 | 0.499 | 0.500 | 0.645 | 0.351 |
| Q_rep vs N_0 | 201 | 195 | 0.995 | 0.907 | 0.988 | 0.585 | 0.552 | 1.000 | 0.001 | 0.474 | 0.370 | 0.328 | 0.326 | 0.392 | 0.696 | 0.879 | 0.500 | 0.999 | 0.477 | 0.500 | 0.999 | 0.693 |
| Q_imp vs N (b <= 0.95) | 401 | 846 | 0.809 | 0.728 | 0.849 | 0.579 | 0.511 | 0.879 | 0.139 | 0.666 | 0.592 | 0.570 | 0.346 | 0.515 | 0.674 | 0.686 | 0.500 | 0.876 | 0.526 | 0.500 | 0.878 | 0.674 |
| Q (no import) vs N_imp | 1324 | 424 | 0.799 | 0.657 | 0.820 | 0.573 | 0.497 | 0.883 | 0.127 | 0.616 | 0.543 | 0.540 | 0.373 | 0.501 | 0.681 | 0.659 | 0.500 | 0.886 | 0.507 | 0.500 | 0.884 | 0.681 |
| Q_ils vs N_0.95 | 200 | 226 | 0.518 | 0.526 | 0.494 | 0.518 | 0.492 | 0.472 | 0.533 | 0.604 | 0.596 | 0.602 | 0.428 | 0.530 | 0.500 | 0.533 | 0.500 | 0.479 | 0.482 | 0.500 | 0.528 | 0.500 |
| decoy: Q_deep vs decoy | 312 | 199 | 0.486 | 0.482 | 0.498 | 0.537 | 0.513 | 0.492 | 0.510 | 0.526 | 0.523 | 0.525 | 0.503 | 0.504 | 0.517 | 0.481 | 0.500 | 0.495 | 0.509 | 0.500 | 0.496 | 0.517 |
| allele genomes 0 | 573 | 498 | 0.761 | 0.686 | 0.797 | 0.522 | 0.519 | 0.844 | 0.164 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.665 | 0.500 | 0.846 | 0.508 | 0.500 | 0.844 | 0.500 |
| allele genomes 2 | 576 | 497 | 0.787 | 0.686 | 0.809 | 0.585 | 0.499 | 0.837 | 0.178 | 0.794 | 0.668 | 0.641 | 0.281 | 0.511 | 0.853 | 0.675 | 0.500 | 0.838 | 0.523 | 0.500 | 0.898 | 0.853 |
| allele genomes 8 | 576 | 502 | 0.786 | 0.681 | 0.807 | 0.593 | 0.495 | 0.843 | 0.169 | 0.768 | 0.601 | 0.581 | 0.250 | 0.504 | 0.860 | 0.648 | 0.500 | 0.841 | 0.483 | 0.500 | 0.893 | 0.860 |
| congeners 2 | 864 | 749 | 0.791 | 0.683 | 0.807 | 0.574 | 0.509 | 0.841 | 0.175 | 0.639 | 0.576 | 0.539 | 0.360 | 0.509 | 0.659 | 0.650 | 0.500 | 0.842 | 0.504 | 0.500 | 0.847 | 0.658 |
| congeners 6 | 861 | 748 | 0.771 | 0.686 | 0.802 | 0.559 | 0.503 | 0.836 | 0.170 | 0.615 | 0.544 | 0.555 | 0.382 | 0.497 | 0.664 | 0.687 | 0.500 | 0.839 | 0.505 | 0.500 | 0.841 | 0.663 |
| twin 0 | 864 | 755 | 0.883 | 0.777 | 0.876 | 0.676 | 0.489 | 0.875 | 0.131 | 0.635 | 0.547 | 0.511 | 0.319 | 0.494 | 0.690 | 0.663 | 0.500 | 0.879 | 0.502 | 0.500 | 0.867 | 0.690 |
| twin 1 | 861 | 742 | 0.778 | 0.601 | 0.721 | 0.548 | 0.537 | 0.798 | 0.217 | 0.621 | 0.572 | 0.580 | 0.410 | 0.513 | 0.640 | 0.677 | 0.500 | 0.798 | 0.505 | 0.500 | 0.820 | 0.639 |
| regime low | 574 | 500 | 0.758 | 0.669 | 0.801 | 0.562 | 0.514 | 0.828 | 0.184 | 0.634 | 0.560 | 0.558 | 0.374 | 0.508 | 0.662 | 0.665 | 0.500 | 0.833 | 0.496 | 0.500 | 0.802 | 0.662 |
| regime high | 575 | 495 | 0.798 | 0.717 | 0.843 | 0.573 | 0.500 | 0.861 | 0.158 | 0.660 | 0.562 | 0.513 | 0.358 | 0.484 | 0.667 | 0.654 | 0.500 | 0.864 | 0.502 | 0.500 | 0.858 | 0.666 |
| regime none | 576 | 502 | 0.778 | 0.670 | 0.774 | 0.565 | 0.506 | 0.839 | 0.167 | 0.589 | 0.554 | 0.565 | 0.379 | 0.519 | 0.657 | 0.668 | 0.500 | 0.836 | 0.513 | 0.500 | 0.880 | 0.656 |
| depth 3 | 440 | 359 | 0.750 | 0.758 | 0.775 | 0.568 | 0.486 | 0.788 | 0.228 | 0.612 | 0.561 | 0.549 | 0.418 | 0.512 | 0.641 | 0.607 | 0.500 | 0.791 | 0.514 | 0.500 | 0.809 | 0.639 |
| depth 10 | 435 | 380 | 0.771 | 0.714 | 0.800 | 0.550 | 0.513 | 0.850 | 0.164 | 0.619 | 0.560 | 0.547 | 0.376 | 0.511 | 0.661 | 0.667 | 0.500 | 0.850 | 0.480 | 0.500 | 0.851 | 0.661 |
| depth 30 | 424 | 381 | 0.799 | 0.675 | 0.825 | 0.569 | 0.504 | 0.864 | 0.146 | 0.633 | 0.557 | 0.543 | 0.325 | 0.494 | 0.672 | 0.687 | 0.500 | 0.866 | 0.516 | 0.500 | 0.861 | 0.672 |
| depth 100 | 426 | 377 | 0.794 | 0.637 | 0.830 | 0.576 | 0.521 | 0.858 | 0.146 | 0.644 | 0.561 | 0.552 | 0.319 | 0.497 | 0.690 | 0.711 | 0.500 | 0.861 | 0.517 | 0.500 | 0.862 | 0.690 |
| all, within identity bands of 0.005 | 1725 | 1497 | 0.564 | 0.548 | 0.656 | 0.680 | 0.698 | 0.731 | 0.287 | 0.685 | 0.681 | 0.637 | 0.462 | 0.588 | 0.635 | 0.550 | 0.500 | 0.730 | 0.503 | 0.500 | 0.767 | 0.635 |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 159 | 98 | 0.799 | 0.898 | 0.976 | 0.743 | 0.675 | 0.787 | 0.971 | 0.786 |
| Q_deep vs N_0 | congeners 6 | 153 | 97 | 0.806 | 0.900 | 0.959 | 0.671 | 0.671 | 0.731 | 0.966 | 0.726 |
| Q_deep vs N_0 | no twin | 162 | 97 | 0.990 | 0.978 | 0.984 | 0.755 | 0.679 | 0.802 | 0.988 | 0.799 |
| Q_deep vs N_0 | twin | 150 | 98 | 0.910 | 0.864 | 0.951 | 0.674 | 0.666 | 0.720 | 0.955 | 0.719 |
| Q_deep vs N_0 | no allele genomes | 83 | 66 | 0.820 | 0.882 | 0.975 | 0.500 | 0.500 | 0.500 | 0.975 | 0.500 |
| Q_deep vs N_0 | allele genomes | 229 | 129 | 0.799 | 0.907 | 0.967 | 0.835 | 0.773 | 0.944 | 0.981 | 0.945 |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 22 | 0.753 | 0.871 | 0.978 | 0.769 | 0.840 | 0.950 | 0.996 | 0.950 |
| Q_deep vs N_0 | recombination none | 107 | 65 | 0.820 | 0.907 | 0.986 | 0.692 | 0.709 | 0.764 | 0.979 | 0.760 |
| Q_deep vs N_0 | recombination high | 102 | 64 | 0.796 | 0.909 | 0.945 | 0.747 | 0.669 | 0.788 | 0.979 | 0.784 |
| Q_deep vs N_0.5 | congeners 2 | 159 | 102 | 0.687 | 0.793 | 0.851 | 0.733 | 0.687 | 0.751 | 0.898 | 0.750 |
| Q_deep vs N_0.5 | congeners 6 | 153 | 104 | 0.714 | 0.834 | 0.878 | 0.637 | 0.639 | 0.696 | 0.874 | 0.696 |
| Q_deep vs N_0.5 | no twin | 162 | 102 | 0.946 | 0.938 | 0.913 | 0.725 | 0.666 | 0.801 | 0.925 | 0.800 |
| Q_deep vs N_0.5 | twin | 150 | 104 | 0.752 | 0.729 | 0.827 | 0.660 | 0.655 | 0.659 | 0.831 | 0.659 |
| Q_deep vs N_0.5 | no allele genomes | 83 | 70 | 0.724 | 0.789 | 0.862 | 0.500 | 0.500 | 0.500 | 0.862 | 0.500 |
| Q_deep vs N_0.5 | allele genomes | 229 | 136 | 0.687 | 0.831 | 0.873 | 0.800 | 0.747 | 0.881 | 0.929 | 0.880 |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.633 | 0.736 | 0.863 | 0.721 | 0.830 | 0.920 | 0.963 | 0.916 |
| Q_deep vs N_0.5 | recombination none | 107 | 68 | 0.711 | 0.801 | 0.886 | 0.684 | 0.713 | 0.748 | 0.923 | 0.747 |
| Q_deep vs N_0.5 | recombination high | 102 | 72 | 0.694 | 0.850 | 0.896 | 0.722 | 0.650 | 0.738 | 0.953 | 0.738 |
| Q_deep vs N_0.8 | congeners 2 | 159 | 113 | 0.604 | 0.666 | 0.674 | 0.702 | 0.644 | 0.699 | 0.761 | 0.699 |
| Q_deep vs N_0.8 | congeners 6 | 153 | 106 | 0.630 | 0.700 | 0.685 | 0.596 | 0.622 | 0.610 | 0.676 | 0.611 |
| Q_deep vs N_0.8 | no twin | 162 | 112 | 0.728 | 0.756 | 0.692 | 0.678 | 0.640 | 0.721 | 0.766 | 0.722 |
| Q_deep vs N_0.8 | twin | 150 | 107 | 0.663 | 0.613 | 0.667 | 0.633 | 0.627 | 0.612 | 0.667 | 0.611 |
| Q_deep vs N_0.8 | no allele genomes | 83 | 69 | 0.635 | 0.671 | 0.689 | 0.500 | 0.500 | 0.500 | 0.689 | 0.500 |
| Q_deep vs N_0.8 | allele genomes | 229 | 150 | 0.610 | 0.695 | 0.683 | 0.756 | 0.718 | 0.756 | 0.777 | 0.756 |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.594 | 0.655 | 0.589 | 0.681 | 0.841 | 0.666 | 0.733 | 0.670 |
| Q_deep vs N_0.8 | recombination none | 107 | 76 | 0.600 | 0.617 | 0.667 | 0.662 | 0.696 | 0.655 | 0.738 | 0.655 |
| Q_deep vs N_0.8 | recombination high | 102 | 75 | 0.625 | 0.771 | 0.721 | 0.677 | 0.625 | 0.676 | 0.753 | 0.676 |
| Q_deep vs N_0.95 | congeners 2 | 159 | 112 | 0.553 | 0.519 | 0.604 | 0.710 | 0.673 | 0.639 | 0.684 | 0.639 |
| Q_deep vs N_0.95 | congeners 6 | 153 | 114 | 0.569 | 0.585 | 0.580 | 0.593 | 0.630 | 0.585 | 0.616 | 0.584 |
| Q_deep vs N_0.95 | no twin | 162 | 116 | 0.579 | 0.590 | 0.593 | 0.656 | 0.652 | 0.664 | 0.671 | 0.664 |
| Q_deep vs N_0.95 | twin | 150 | 110 | 0.599 | 0.506 | 0.588 | 0.647 | 0.646 | 0.587 | 0.638 | 0.585 |
| Q_deep vs N_0.95 | no allele genomes | 83 | 76 | 0.579 | 0.526 | 0.603 | 0.500 | 0.500 | 0.500 | 0.603 | 0.500 |
| Q_deep vs N_0.95 | allele genomes | 229 | 150 | 0.562 | 0.570 | 0.594 | 0.740 | 0.736 | 0.660 | 0.693 | 0.661 |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.570 | 0.661 | 0.581 | 0.645 | 0.864 | 0.613 | 0.627 | 0.612 |
| Q_deep vs N_0.95 | recombination none | 107 | 75 | 0.583 | 0.547 | 0.616 | 0.659 | 0.712 | 0.634 | 0.680 | 0.634 |
| Q_deep vs N_0.95 | recombination high | 102 | 73 | 0.564 | 0.618 | 0.642 | 0.685 | 0.649 | 0.632 | 0.721 | 0.632 |

## L2c: separation, se (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1717 | 1491 | 0.758 | 0.680 | 0.774 | 0.562 | 0.518 | 0.801 | 0.215 | 0.625 | 0.564 | 0.547 | 0.394 | 0.497 | 0.643 | 0.622 | 0.500 | 0.802 | 0.506 | 0.500 | 0.809 | - |
| Q_deep vs N_0 | 310 | 201 | 0.753 | 0.755 | 0.845 | 0.585 | 0.613 | 0.919 | 0.106 | 0.711 | 0.672 | 0.667 | 0.423 | 0.638 | 0.722 | 0.558 | 0.500 | 0.927 | 0.494 | 0.500 | 0.907 | - |
| Q_deep vs N_0.5 | 310 | 204 | 0.664 | 0.662 | 0.752 | 0.567 | 0.554 | 0.786 | 0.237 | 0.699 | 0.662 | 0.654 | 0.459 | 0.612 | 0.692 | 0.563 | 0.500 | 0.795 | 0.485 | 0.500 | 0.810 | - |
| Q_deep vs N_0.8 | 310 | 214 | 0.600 | 0.607 | 0.643 | 0.545 | 0.504 | 0.636 | 0.371 | 0.653 | 0.629 | 0.629 | 0.439 | 0.556 | 0.622 | 0.559 | 0.500 | 0.644 | 0.512 | 0.500 | 0.679 | - |
| Q_deep vs N_0.95 | 310 | 227 | 0.512 | 0.551 | 0.520 | 0.533 | 0.514 | 0.532 | 0.464 | 0.658 | 0.641 | 0.650 | 0.463 | 0.556 | 0.592 | 0.507 | 0.500 | 0.536 | 0.500 | 0.500 | 0.602 | - |
| Q_deep vs N_1 | 310 | 225 | 0.503 | 0.520 | 0.501 | 0.535 | 0.513 | 0.486 | 0.515 | 0.640 | 0.627 | 0.631 | 0.467 | 0.540 | 0.541 | 0.528 | 0.500 | 0.484 | 0.486 | 0.500 | 0.550 | - |
| Q_near vs N_0.95 | 331 | 227 | 0.722 | 0.652 | 0.696 | 0.582 | 0.501 | 0.736 | 0.278 | 0.886 | 0.772 | 0.742 | 0.411 | 0.606 | 0.890 | 0.599 | 0.500 | 0.734 | 0.516 | 0.500 | 0.896 | - |
| Q_lone vs N_0.95 | 276 | 227 | 0.733 | 0.655 | 0.751 | 0.510 | 0.476 | 0.763 | 0.239 | 0.350 | 0.329 | 0.320 | 0.344 | 0.298 | 0.345 | 0.610 | 0.500 | 0.762 | 0.514 | 0.500 | 0.595 | - |
| Q_rep vs N_0 | 201 | 201 | 0.972 | 0.889 | 0.976 | 0.624 | 0.594 | 0.982 | 0.033 | 0.482 | 0.393 | 0.346 | 0.349 | 0.382 | 0.681 | 0.808 | 0.500 | 0.983 | 0.472 | 0.500 | 0.977 | - |
| Q_imp vs N (b <= 0.95) | 401 | 846 | 0.790 | 0.706 | 0.809 | 0.553 | 0.527 | 0.838 | 0.185 | 0.663 | 0.597 | 0.573 | 0.366 | 0.512 | 0.662 | 0.630 | 0.500 | 0.838 | 0.524 | 0.500 | 0.847 | - |
| Q (no import) vs N_imp | 1316 | 420 | 0.779 | 0.676 | 0.794 | 0.573 | 0.512 | 0.843 | 0.170 | 0.613 | 0.547 | 0.539 | 0.407 | 0.501 | 0.665 | 0.625 | 0.500 | 0.843 | 0.505 | 0.500 | 0.852 | - |
| Q_ils vs N_0.95 | 198 | 227 | 0.498 | 0.524 | 0.472 | 0.506 | 0.486 | 0.446 | 0.561 | 0.603 | 0.593 | 0.600 | 0.450 | 0.530 | 0.501 | 0.481 | 0.500 | 0.450 | 0.494 | 0.500 | 0.494 | - |
| decoy: Q_deep vs decoy | 310 | 198 | 0.473 | 0.494 | 0.493 | 0.535 | 0.507 | 0.485 | 0.505 | 0.521 | 0.527 | 0.525 | 0.526 | 0.502 | 0.510 | 0.531 | 0.500 | 0.484 | 0.507 | 0.500 | 0.497 | - |
| allele genomes 0 | 571 | 496 | 0.752 | 0.686 | 0.782 | 0.527 | 0.525 | 0.815 | 0.190 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.616 | 0.500 | 0.817 | 0.519 | 0.500 | 0.815 | - |
| allele genomes 2 | 574 | 498 | 0.762 | 0.677 | 0.774 | 0.571 | 0.525 | 0.790 | 0.235 | 0.781 | 0.663 | 0.641 | 0.326 | 0.508 | 0.813 | 0.637 | 0.500 | 0.791 | 0.525 | 0.500 | 0.847 | - |
| allele genomes 8 | 572 | 497 | 0.761 | 0.676 | 0.765 | 0.591 | 0.506 | 0.808 | 0.209 | 0.768 | 0.615 | 0.585 | 0.316 | 0.492 | 0.819 | 0.616 | 0.500 | 0.805 | 0.474 | 0.500 | 0.860 | - |
| congeners 2 | 860 | 743 | 0.767 | 0.679 | 0.772 | 0.561 | 0.527 | 0.802 | 0.216 | 0.638 | 0.576 | 0.542 | 0.378 | 0.493 | 0.637 | 0.623 | 0.500 | 0.802 | 0.504 | 0.500 | 0.807 | - |
| congeners 6 | 857 | 748 | 0.755 | 0.681 | 0.775 | 0.562 | 0.511 | 0.802 | 0.211 | 0.613 | 0.553 | 0.552 | 0.411 | 0.500 | 0.650 | 0.631 | 0.500 | 0.803 | 0.507 | 0.500 | 0.813 | - |
| twin 0 | 862 | 750 | 0.862 | 0.764 | 0.843 | 0.634 | 0.507 | 0.850 | 0.158 | 0.629 | 0.555 | 0.513 | 0.343 | 0.490 | 0.673 | 0.622 | 0.500 | 0.853 | 0.516 | 0.500 | 0.845 | - |
| twin 1 | 855 | 741 | 0.744 | 0.603 | 0.697 | 0.564 | 0.553 | 0.752 | 0.270 | 0.624 | 0.577 | 0.577 | 0.433 | 0.504 | 0.622 | 0.633 | 0.500 | 0.750 | 0.497 | 0.500 | 0.776 | - |
| regime low | 574 | 501 | 0.735 | 0.669 | 0.774 | 0.546 | 0.528 | 0.795 | 0.220 | 0.634 | 0.566 | 0.561 | 0.400 | 0.508 | 0.645 | 0.626 | 0.500 | 0.798 | 0.497 | 0.500 | 0.783 | - |
| regime high | 570 | 489 | 0.775 | 0.695 | 0.801 | 0.567 | 0.502 | 0.817 | 0.203 | 0.654 | 0.568 | 0.519 | 0.380 | 0.474 | 0.636 | 0.615 | 0.500 | 0.818 | 0.508 | 0.500 | 0.814 | - |
| regime none | 573 | 501 | 0.766 | 0.677 | 0.752 | 0.576 | 0.526 | 0.800 | 0.215 | 0.590 | 0.558 | 0.558 | 0.402 | 0.511 | 0.649 | 0.628 | 0.500 | 0.798 | 0.513 | 0.500 | 0.835 | - |
| depth 3 | 432 | 353 | 0.717 | 0.737 | 0.716 | 0.542 | 0.506 | 0.720 | 0.304 | 0.617 | 0.583 | 0.549 | 0.439 | 0.516 | 0.608 | 0.540 | 0.500 | 0.720 | 0.494 | 0.500 | 0.738 | - |
| depth 10 | 435 | 380 | 0.754 | 0.740 | 0.771 | 0.553 | 0.523 | 0.814 | 0.206 | 0.619 | 0.564 | 0.547 | 0.415 | 0.499 | 0.647 | 0.634 | 0.500 | 0.814 | 0.488 | 0.500 | 0.816 | - |
| depth 30 | 424 | 381 | 0.791 | 0.686 | 0.819 | 0.567 | 0.518 | 0.848 | 0.163 | 0.633 | 0.553 | 0.542 | 0.345 | 0.484 | 0.656 | 0.658 | 0.500 | 0.847 | 0.544 | 0.500 | 0.848 | - |
| depth 100 | 426 | 377 | 0.779 | 0.628 | 0.814 | 0.595 | 0.529 | 0.846 | 0.161 | 0.647 | 0.561 | 0.554 | 0.341 | 0.490 | 0.686 | 0.697 | 0.500 | 0.850 | 0.510 | 0.500 | 0.852 | - |
| all, within identity bands of 0.005 | 1717 | 1491 | 0.558 | 0.554 | 0.641 | 0.663 | 0.691 | 0.698 | 0.323 | 0.687 | 0.684 | 0.640 | 0.474 | 0.581 | 0.620 | 0.536 | 0.500 | 0.699 | 0.498 | 0.500 | 0.732 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 158 | 99 | 0.720 | 0.814 | 0.909 | 0.743 | 0.686 | 0.724 | 0.885 | - |
| Q_deep vs N_0 | congeners 6 | 152 | 102 | 0.783 | 0.871 | 0.930 | 0.676 | 0.650 | 0.721 | 0.934 | - |
| Q_deep vs N_0 | no twin | 161 | 100 | 0.950 | 0.909 | 0.947 | 0.761 | 0.686 | 0.787 | 0.953 | - |
| Q_deep vs N_0 | twin | 149 | 101 | 0.830 | 0.825 | 0.893 | 0.671 | 0.648 | 0.671 | 0.866 | - |
| Q_deep vs N_0 | no allele genomes | 82 | 66 | 0.779 | 0.847 | 0.928 | 0.500 | 0.500 | 0.500 | 0.928 | - |
| Q_deep vs N_0 | allele genomes | 228 | 135 | 0.744 | 0.849 | 0.919 | 0.831 | 0.757 | 0.866 | 0.915 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 22 | 0.724 | 0.758 | 0.912 | 0.759 | 0.799 | 0.928 | 0.960 | - |
| Q_deep vs N_0 | recombination none | 107 | 67 | 0.772 | 0.846 | 0.956 | 0.697 | 0.701 | 0.745 | 0.919 | - |
| Q_deep vs N_0 | recombination high | 101 | 65 | 0.739 | 0.845 | 0.877 | 0.756 | 0.668 | 0.737 | 0.898 | - |
| Q_deep vs N_0.5 | congeners 2 | 158 | 100 | 0.648 | 0.728 | 0.764 | 0.745 | 0.687 | 0.722 | 0.824 | - |
| Q_deep vs N_0.5 | congeners 6 | 152 | 104 | 0.685 | 0.774 | 0.810 | 0.649 | 0.625 | 0.660 | 0.801 | - |
| Q_deep vs N_0.5 | no twin | 161 | 100 | 0.867 | 0.871 | 0.858 | 0.732 | 0.667 | 0.765 | 0.883 | - |
| Q_deep vs N_0.5 | twin | 149 | 104 | 0.701 | 0.670 | 0.723 | 0.675 | 0.638 | 0.639 | 0.739 | - |
| Q_deep vs N_0.5 | no allele genomes | 82 | 70 | 0.690 | 0.767 | 0.804 | 0.500 | 0.500 | 0.500 | 0.804 | - |
| Q_deep vs N_0.5 | allele genomes | 228 | 134 | 0.660 | 0.753 | 0.786 | 0.813 | 0.727 | 0.819 | 0.837 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.593 | 0.615 | 0.773 | 0.743 | 0.774 | 0.831 | 0.812 | - |
| Q_deep vs N_0.5 | recombination none | 107 | 67 | 0.676 | 0.739 | 0.830 | 0.699 | 0.699 | 0.719 | 0.857 | - |
| Q_deep vs N_0.5 | recombination high | 101 | 71 | 0.664 | 0.776 | 0.795 | 0.746 | 0.647 | 0.713 | 0.856 | - |
| Q_deep vs N_0.8 | congeners 2 | 158 | 109 | 0.593 | 0.634 | 0.624 | 0.695 | 0.657 | 0.659 | 0.700 | - |
| Q_deep vs N_0.8 | congeners 6 | 152 | 105 | 0.616 | 0.650 | 0.648 | 0.606 | 0.601 | 0.587 | 0.664 | - |
| Q_deep vs N_0.8 | no twin | 161 | 111 | 0.686 | 0.688 | 0.663 | 0.686 | 0.645 | 0.671 | 0.711 | - |
| Q_deep vs N_0.8 | twin | 149 | 103 | 0.629 | 0.603 | 0.609 | 0.626 | 0.611 | 0.594 | 0.651 | - |
| Q_deep vs N_0.8 | no allele genomes | 82 | 67 | 0.665 | 0.706 | 0.682 | 0.500 | 0.500 | 0.500 | 0.682 | - |
| Q_deep vs N_0.8 | allele genomes | 228 | 147 | 0.578 | 0.620 | 0.623 | 0.754 | 0.711 | 0.696 | 0.710 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.586 | 0.576 | 0.592 | 0.710 | 0.746 | 0.598 | 0.685 | - |
| Q_deep vs N_0.8 | recombination none | 107 | 73 | 0.607 | 0.612 | 0.667 | 0.652 | 0.663 | 0.623 | 0.717 | - |
| Q_deep vs N_0.8 | recombination high | 101 | 73 | 0.584 | 0.649 | 0.625 | 0.677 | 0.624 | 0.610 | 0.684 | - |
| Q_deep vs N_0.95 | congeners 2 | 158 | 113 | 0.512 | 0.492 | 0.548 | 0.701 | 0.670 | 0.614 | 0.627 | - |
| Q_deep vs N_0.95 | congeners 6 | 152 | 114 | 0.511 | 0.549 | 0.518 | 0.618 | 0.628 | 0.569 | 0.576 | - |
| Q_deep vs N_0.95 | no twin | 161 | 116 | 0.542 | 0.545 | 0.553 | 0.663 | 0.652 | 0.647 | 0.641 | - |
| Q_deep vs N_0.95 | twin | 149 | 111 | 0.503 | 0.484 | 0.512 | 0.657 | 0.646 | 0.558 | 0.570 | - |
| Q_deep vs N_0.95 | no allele genomes | 82 | 78 | 0.534 | 0.544 | 0.540 | 0.500 | 0.500 | 0.500 | 0.540 | - |
| Q_deep vs N_0.95 | allele genomes | 228 | 149 | 0.519 | 0.518 | 0.543 | 0.740 | 0.727 | 0.614 | 0.627 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.530 | 0.583 | 0.515 | 0.680 | 0.800 | 0.589 | 0.573 | - |
| Q_deep vs N_0.95 | recombination none | 107 | 77 | 0.524 | 0.526 | 0.546 | 0.673 | 0.698 | 0.626 | 0.629 | - |
| Q_deep vs N_0.95 | recombination high | 101 | 73 | 0.528 | 0.553 | 0.549 | 0.686 | 0.648 | 0.592 | 0.642 | - |

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
| pe | identity, depth +ancestry | 0.780 | 0.836 | 0.056 | 0.034 | 0.073 |
| pe | identity, depth +alleles | 0.780 | 0.853 | 0.073 | 0.035 | 0.119 |
| pe | identity, depth +polymorphic | 0.780 | 0.813 | 0.033 | 0.002 | 0.063 |
| pe | identity, depth +weights | 0.780 | 0.840 | 0.060 | 0.044 | 0.073 |
| pe | identity, depth +all three | 0.780 | 0.914 | 0.134 | 0.099 | 0.171 |
| pe | identity, depth +ancestry, polymorphic | 0.780 | 0.894 | 0.114 | 0.084 | 0.147 |
| pe | identity, depth +all four | 0.780 | 0.915 | 0.135 | 0.098 | 0.173 |
| pe | identity, depth +oracle fixed agreement | 0.780 | 0.781 | 0.001 | -0.006 | 0.009 |
| pe, genera with alleles | identity, depth +ancestry | 0.787 | 0.840 | 0.053 | 0.028 | 0.076 |
| pe, genera with alleles | identity, depth +alleles | 0.787 | 0.922 | 0.134 | 0.075 | 0.201 |
| pe, genera with alleles | identity, depth +polymorphic | 0.787 | 0.903 | 0.116 | 0.071 | 0.167 |
| pe, genera with alleles | identity, depth +weights | 0.787 | 0.844 | 0.057 | 0.034 | 0.077 |
| pe, genera with alleles | identity, depth +all three | 0.787 | 0.947 | 0.160 | 0.100 | 0.218 |
| pe, genera with alleles | identity, depth +ancestry, polymorphic | 0.787 | 0.925 | 0.137 | 0.089 | 0.189 |
| pe, genera with alleles | identity, depth +all four | 0.787 | 0.948 | 0.160 | 0.101 | 0.220 |
| pe, genera with alleles | identity, depth +oracle fixed agreement | 0.787 | 0.867 | 0.079 | 0.054 | 0.104 |
| se | identity, depth +ancestry | 0.761 | 0.810 | 0.050 | 0.030 | 0.064 |
| se | identity, depth +alleles | 0.761 | 0.843 | 0.082 | 0.045 | 0.130 |
| se | identity, depth +polymorphic | 0.761 | 0.791 | 0.030 | 0.005 | 0.058 |
| se | identity, depth +weights | 0.761 | 0.806 | 0.046 | 0.031 | 0.056 |
| se | identity, depth +all three | 0.761 | 0.900 | 0.139 | 0.106 | 0.173 |
| se | identity, depth +ancestry, polymorphic | 0.761 | 0.866 | 0.105 | 0.075 | 0.135 |
| se | identity, depth +all four | 0.761 | 0.900 | 0.140 | 0.106 | 0.173 |
| se | identity, depth +oracle fixed agreement | 0.761 | 0.761 | 0.000 | 0.000 | 0.000 |
| se, genera with alleles | identity, depth +ancestry | 0.763 | 0.806 | 0.043 | 0.020 | 0.064 |
| se, genera with alleles | identity, depth +alleles | 0.763 | 0.910 | 0.148 | 0.087 | 0.206 |
| se, genera with alleles | identity, depth +polymorphic | 0.763 | 0.866 | 0.103 | 0.062 | 0.148 |
| se, genera with alleles | identity, depth +weights | 0.763 | 0.805 | 0.042 | 0.025 | 0.058 |
| se, genera with alleles | identity, depth +all three | 0.763 | 0.933 | 0.170 | 0.113 | 0.225 |
| se, genera with alleles | identity, depth +ancestry, polymorphic | 0.763 | 0.893 | 0.130 | 0.086 | 0.174 |
| se, genera with alleles | identity, depth +all four | 0.763 | 0.933 | 0.170 | 0.113 | 0.224 |
| se, genera with alleles | identity, depth +oracle fixed agreement | 0.763 | 0.763 | 0.000 | 0.000 | 0.000 |
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
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 279 | 0.265 | 0.355 | 0.452 | 0.341 | 0.495 | 0.000 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.078 | 0.109 | 0.062 | 0.078 | 0.109 | 0.000 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 573 | 0.377 | 0.529 | 0.000 | 0.000 | 0.154 | 0.000 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 809 | 0.396 | 0.533 | 0.743 | 0.477 | 0.791 | 0.000 |
| pe | Q | depth | 10 | 435 | 0.361 | 0.478 | 0.421 | 0.278 | 0.480 | 0.000 |
| pe | Q | depth | 100 | 426 | 0.338 | 0.467 | 0.423 | 0.303 | 0.521 | 0.000 |
| pe | Q | depth | 3 | 440 | 0.398 | 0.523 | 0.423 | 0.282 | 0.536 | 0.000 |
| pe | Q | depth | 30 | 424 | 0.328 | 0.479 | 0.429 | 0.264 | 0.486 | 0.000 |
| pe | Q | distance to rep | 0.005-0.015 | 773 | 0.362 | 0.503 | 0.545 | 0.310 | 0.512 | 0.000 |
| pe | Q | distance to rep | 0.015-0.03 | 473 | 0.047 | 0.279 | 0.493 | 0.370 | 0.362 | 0.000 |
| pe | Q | distance to rep | <0.005 | 321 | 0.975 | 0.991 | 0.084 | 0.075 | 0.860 | 0.000 |
| pe | Q | distance to rep | >=0.03 | 158 | 0.000 | 0.006 | 0.316 | 0.297 | 0.190 | 0.000 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 567 | 0.354 | 0.517 | 0.471 | 0.291 | 0.492 | 0.000 |
| pe | Q | genes with a congener's segment | <=0.1 | 464 | 0.317 | 0.448 | 0.384 | 0.216 | 0.483 | 0.000 |
| pe | Q | genes with a congener's segment | >0.3 | 257 | 0.327 | 0.514 | 0.502 | 0.296 | 0.444 | 0.000 |
| pe | Q | genes with a congener's segment | none | 437 | 0.419 | 0.474 | 0.359 | 0.332 | 0.586 | 0.000 |
| pe | Q | regime | high | 575 | 0.360 | 0.539 | 0.480 | 0.256 | 0.438 | 0.000 |
| pe | Q | regime | low | 574 | 0.300 | 0.449 | 0.413 | 0.267 | 0.493 | 0.000 |
| pe | Q | regime | none | 576 | 0.410 | 0.472 | 0.378 | 0.323 | 0.587 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 949 | 0.348 | 0.474 | 0.749 | 0.479 | 0.765 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 582 | 0.375 | 0.521 | 0.000 | 0.005 | 0.139 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 183 | 0.306 | 0.415 | 0.109 | 0.153 | 0.301 | 0.000 |
| pe | Q | role | Q_deep | 312 | 0.112 | 0.186 | 0.429 | 0.372 | 0.314 | 0.000 |
| pe | Q | role | Q_ils | 200 | 0.075 | 0.070 | 0.405 | 0.325 | 0.110 | 0.000 |
| pe | Q | role | Q_imp | 401 | 0.322 | 0.504 | 0.506 | 0.314 | 0.511 | 0.000 |
| pe | Q | role | Q_lone | 277 | 0.437 | 0.635 | 0.126 | 0.094 | 0.329 | 0.000 |
| pe | Q | role | Q_near | 334 | 0.419 | 0.581 | 0.746 | 0.458 | 0.886 | 0.000 |
| pe | Q | role | Q_rep | 201 | 0.871 | 0.975 | 0.144 | 0.000 | 0.801 | 0.000 |
| pe | Q | twin genus | 0 | 864 | 0.251 | 0.519 | 0.416 | 0.229 | 0.486 | 0.000 |
| pe | Q | twin genus | 1 | 861 | 0.462 | 0.455 | 0.432 | 0.334 | 0.526 | 0.000 |
| pe | N | depth | 10 | 380 | 0.050 | 0.058 | 0.058 | 0.058 | 0.053 | 0.000 |
| pe | N | depth | 100 | 377 | 0.011 | 0.032 | 0.000 | 0.008 | 0.019 | 0.000 |
| pe | N | depth | 3 | 359 | 0.123 | 0.081 | 0.123 | 0.092 | 0.106 | 0.000 |
| pe | N | depth | 30 | 381 | 0.021 | 0.031 | 0.024 | 0.029 | 0.026 | 0.000 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 875 | 0.048 | 0.054 | 0.049 | 0.050 | 0.059 | 0.000 |
| pe | N | genes with a segment of S or a congener | >0.3 | 255 | 0.059 | 0.016 | 0.059 | 0.059 | 0.016 | 0.000 |
| pe | N | genes with a segment of S or a congener | none | 367 | 0.049 | 0.065 | 0.046 | 0.027 | 0.052 | 0.000 |
| pe | N | regime | high | 495 | 0.055 | 0.048 | 0.046 | 0.053 | 0.048 | 0.000 |
| pe | N | regime | low | 500 | 0.040 | 0.050 | 0.064 | 0.054 | 0.058 | 0.000 |
| pe | N | regime | none | 502 | 0.056 | 0.052 | 0.040 | 0.032 | 0.044 | 0.000 |
| pe | N | role | N_0 | 195 | 0.005 | 0.005 | 0.031 | 0.026 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 206 | 0.049 | 0.015 | 0.044 | 0.044 | 0.010 | 0.000 |
| pe | N | role | N_0.8 | 219 | 0.027 | 0.032 | 0.064 | 0.064 | 0.059 | 0.000 |
| pe | N | role | N_0.95 | 226 | 0.053 | 0.088 | 0.066 | 0.044 | 0.062 | 0.000 |
| pe | N | role | N_1 | 227 | 0.115 | 0.163 | 0.057 | 0.031 | 0.163 | 0.000 |
| pe | N | role | N_imp0.2 | 213 | 0.028 | 0.019 | 0.028 | 0.056 | 0.033 | 0.000 |
| pe | N | role | N_imp0.4 | 211 | 0.066 | 0.014 | 0.057 | 0.057 | 0.009 | 0.000 |
| pe | N | twin genus | 0 | 755 | 0.000 | 0.038 | 0.024 | 0.037 | 0.023 | 0.000 |
| pe | N | twin genus | 1 | 742 | 0.101 | 0.062 | 0.077 | 0.055 | 0.078 | 0.000 |
| se | Q | MRCA with nearest allele genome | 0.002-0.01 | 277 | 0.264 | 0.274 | 0.448 | 0.303 | 0.051 |  |
| se | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.078 | 0.078 | 0.062 | 0.078 | 0.000 |  |
| se | Q | MRCA with nearest allele genome | no allele genome | 571 | 0.352 | 0.417 | 0.000 | 0.000 | 0.000 |  |
| se | Q | MRCA with nearest allele genome | recent (<0.002) | 805 | 0.352 | 0.431 | 0.750 | 0.445 | 0.130 |  |
| se | Q | depth | 10 | 435 | 0.340 | 0.375 | 0.425 | 0.253 | 0.124 |  |
| se | Q | depth | 100 | 426 | 0.300 | 0.376 | 0.434 | 0.272 | 0.031 |  |
| se | Q | depth | 3 | 432 | 0.354 | 0.424 | 0.403 | 0.266 | 0.032 |  |
| se | Q | depth | 30 | 424 | 0.314 | 0.377 | 0.443 | 0.250 | 0.090 |  |
| se | Q | distance to rep | 0.005-0.015 | 766 | 0.305 | 0.376 | 0.539 | 0.291 | 0.093 |  |
| se | Q | distance to rep | 0.015-0.03 | 473 | 0.057 | 0.152 | 0.488 | 0.328 | 0.032 |  |
| se | Q | distance to rep | <0.005 | 320 | 0.941 | 0.953 | 0.119 | 0.081 | 0.097 |  |
| se | Q | distance to rep | >=0.03 | 158 | 0.000 | 0.006 | 0.316 | 0.272 | 0.013 |  |
| se | Q | genes with a congener's segment | 0.1-0.3 | 562 | 0.315 | 0.384 | 0.473 | 0.283 | 0.091 |  |
| se | Q | genes with a congener's segment | <=0.1 | 463 | 0.289 | 0.343 | 0.397 | 0.192 | 0.060 |  |
| se | Q | genes with a congener's segment | >0.3 | 257 | 0.300 | 0.428 | 0.486 | 0.261 | 0.054 |  |
| se | Q | genes with a congener's segment | none | 435 | 0.400 | 0.416 | 0.361 | 0.303 | 0.060 |  |
| se | Q | regime | high | 570 | 0.314 | 0.411 | 0.467 | 0.239 | 0.086 |  |
| se | Q | regime | low | 574 | 0.279 | 0.345 | 0.434 | 0.235 | 0.066 |  |
| se | Q | regime | none | 573 | 0.389 | 0.408 | 0.379 | 0.307 | 0.056 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 945 | 0.312 | 0.385 | 0.747 | 0.439 | 0.114 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 579 | 0.351 | 0.402 | 0.000 | 0.005 | 0.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | partly | 182 | 0.291 | 0.319 | 0.143 | 0.159 | 0.060 |  |
| se | Q | role | Q_deep | 310 | 0.113 | 0.081 | 0.423 | 0.339 | 0.042 |  |
| se | Q | role | Q_ils | 198 | 0.066 | 0.030 | 0.399 | 0.288 | 0.030 |  |
| se | Q | role | Q_imp | 401 | 0.304 | 0.411 | 0.501 | 0.292 | 0.072 |  |
| se | Q | role | Q_lone | 276 | 0.377 | 0.478 | 0.127 | 0.098 | 0.025 |  |
| se | Q | role | Q_near | 331 | 0.366 | 0.471 | 0.773 | 0.417 | 0.148 |  |
| se | Q | role | Q_rep | 201 | 0.831 | 0.905 | 0.149 | 0.015 | 0.075 |  |
| se | Q | twin genus | 0 | 862 | 0.225 | 0.393 | 0.412 | 0.208 | 0.044 |  |
| se | Q | twin genus | 1 | 855 | 0.430 | 0.382 | 0.441 | 0.313 | 0.095 |  |
| se | N | depth | 10 | 380 | 0.063 | 0.047 | 0.055 | 0.074 | 0.029 |  |
| se | N | depth | 100 | 377 | 0.003 | 0.003 | 0.005 | 0.008 | 0.000 |  |
| se | N | depth | 3 | 353 | 0.119 | 0.130 | 0.102 | 0.099 | 0.025 |  |
| se | N | depth | 30 | 381 | 0.018 | 0.010 | 0.039 | 0.024 | 0.008 |  |
| se | N | genes with a segment of S or a congener | <=0.3 | 872 | 0.054 | 0.055 | 0.044 | 0.056 | 0.018 |  |
| se | N | genes with a segment of S or a congener | >0.3 | 251 | 0.040 | 0.016 | 0.060 | 0.056 | 0.012 |  |
| se | N | genes with a segment of S or a congener | none | 368 | 0.046 | 0.046 | 0.057 | 0.033 | 0.011 |  |
| se | N | regime | high | 489 | 0.055 | 0.059 | 0.037 | 0.063 | 0.020 |  |
| se | N | regime | low | 501 | 0.052 | 0.040 | 0.062 | 0.052 | 0.014 |  |
| se | N | regime | none | 501 | 0.042 | 0.040 | 0.050 | 0.036 | 0.012 |  |
| se | N | role | N_0 | 201 | 0.045 | 0.025 | 0.025 | 0.030 | 0.000 |  |
| se | N | role | N_0.5 | 204 | 0.049 | 0.044 | 0.025 | 0.049 | 0.010 |  |
| se | N | role | N_0.8 | 214 | 0.014 | 0.042 | 0.047 | 0.061 | 0.023 |  |
| se | N | role | N_0.95 | 227 | 0.093 | 0.079 | 0.079 | 0.048 | 0.035 |  |
| se | N | role | N_1 | 225 | 0.084 | 0.093 | 0.080 | 0.058 | 0.018 |  |
| se | N | role | N_imp0.2 | 212 | 0.014 | 0.019 | 0.019 | 0.052 | 0.009 |  |
| se | N | role | N_imp0.4 | 208 | 0.043 | 0.014 | 0.067 | 0.053 | 0.010 |  |
| se | N | twin genus | 0 | 750 | 0.004 | 0.020 | 0.028 | 0.041 | 0.001 |  |
| se | N | twin genus | 1 | 741 | 0.096 | 0.073 | 0.072 | 0.059 | 0.030 |  |
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
| pe_alleles | target rows | 1725 | 1497 | 0.923 | 0.730 | 0.830 | 0.139 |  |
| pe_alleles | hard: Q_deep vs N_0.8/0.95 | 312 | 445 | 0.791 | 0.449 | 0.657 | 0.207 |  |
| pe_alleles | twin genera | 861 | 742 | 0.926 | 0.726 | 0.851 | 0.168 |  |
| pe_alleles | no allele genomes | 573 | 498 | 0.822 | 0.431 | 0.696 | 0.211 |  |
| pe_alleles | recombination high | 575 | 495 | 0.939 | 0.744 | 0.856 | 0.156 |  |
| pe_alleles | all test rows | 2685 | 3651 | 0.976 | - | 0.891 | 0.062 | 0.902 |
| pe_ancestry+alleles | target rows | 1725 | 1497 | 0.922 | 0.734 | 0.823 | 0.139 |  |
| pe_ancestry+alleles | hard: Q_deep vs N_0.8/0.95 | 312 | 445 | 0.795 | 0.465 | 0.670 | 0.204 |  |
| pe_ancestry+alleles | twin genera | 861 | 742 | 0.923 | 0.739 | 0.846 | 0.171 |  |
| pe_ancestry+alleles | no allele genomes | 573 | 498 | 0.828 | 0.473 | 0.667 | 0.187 |  |
| pe_ancestry+alleles | recombination high | 575 | 495 | 0.940 | 0.758 | 0.850 | 0.145 |  |
| pe_ancestry+alleles | all test rows | 2685 | 3651 | 0.977 | - | 0.886 | 0.058 | 0.902 |
| pe_ancestry | target rows | 1725 | 1497 | 0.858 | 0.506 | 0.718 | 0.182 |  |
| pe_ancestry | hard: Q_deep vs N_0.8/0.95 | 312 | 445 | 0.667 | 0.212 | 0.526 | 0.301 |  |
| pe_ancestry | twin genera | 861 | 742 | 0.824 | 0.436 | 0.656 | 0.198 |  |
| pe_ancestry | no allele genomes | 573 | 498 | 0.829 | 0.473 | 0.691 | 0.185 |  |
| pe_ancestry | recombination high | 575 | 495 | 0.884 | 0.553 | 0.786 | 0.206 |  |
| pe_ancestry | all test rows | 2685 | 3651 | 0.960 | - | 0.819 | 0.076 | 0.852 |
| pe_base | target rows | 1725 | 1497 | 0.856 | 0.504 | 0.711 | 0.176 |  |
| pe_base | hard: Q_deep vs N_0.8/0.95 | 312 | 445 | 0.669 | 0.212 | 0.519 | 0.281 |  |
| pe_base | twin genera | 861 | 742 | 0.821 | 0.439 | 0.649 | 0.189 |  |
| pe_base | no allele genomes | 573 | 498 | 0.828 | 0.462 | 0.670 | 0.167 |  |
| pe_base | recombination high | 575 | 495 | 0.887 | 0.555 | 0.774 | 0.188 |  |
| pe_base | all test rows | 2685 | 3651 | 0.957 | - | 0.815 | 0.078 | 0.848 |
| pe_default | target rows | 1725 | 1497 | 0.925 | 0.747 | 0.823 | 0.126 |  |
| pe_default | hard: Q_deep vs N_0.8/0.95 | 312 | 445 | 0.798 | 0.484 | 0.654 | 0.184 |  |
| pe_default | twin genera | 861 | 742 | 0.931 | 0.748 | 0.864 | 0.166 |  |
| pe_default | no allele genomes | 573 | 498 | 0.824 | 0.468 | 0.682 | 0.199 |  |
| pe_default | recombination high | 575 | 495 | 0.942 | 0.770 | 0.838 | 0.135 |  |
| pe_default | all test rows | 2685 | 3651 | 0.978 | - | 0.886 | 0.053 | 0.905 |
| pe_default_on_noallele | target rows | 1719 | 1476 | 0.917 | 0.721 | 0.818 | 0.136 |  |
| pe_default_on_noallele | hard: Q_deep vs N_0.8/0.95 | 311 | 439 | 0.784 | 0.437 | 0.633 | 0.196 |  |
| pe_default_on_noallele | twin genera | 855 | 724 | 0.917 | 0.720 | 0.854 | 0.185 |  |
| pe_default_on_noallele | no allele genomes | 573 | 498 | 0.823 | 0.469 | 0.677 | 0.211 |  |
| pe_default_on_noallele | recombination high | 574 | 487 | 0.933 | 0.730 | 0.840 | 0.148 |  |
| pe_default_on_noallele | all test rows | 2679 | 3858 | 0.977 | - | 0.883 | 0.054 | 0.901 |
| pe_default_on_shuffled | target rows | 1720 | 1482 | 0.837 | 0.456 | 0.583 | 0.119 |  |
| pe_default_on_shuffled | hard: Q_deep vs N_0.8/0.95 | 311 | 441 | 0.620 | 0.129 | 0.315 | 0.177 |  |
| pe_default_on_shuffled | twin genera | 856 | 730 | 0.810 | 0.436 | 0.611 | 0.168 |  |
| pe_default_on_shuffled | no allele genomes | 573 | 498 | 0.823 | 0.469 | 0.679 | 0.211 |  |
| pe_default_on_shuffled | recombination high | 575 | 487 | 0.865 | 0.497 | 0.638 | 0.121 |  |
| pe_default_on_shuffled | all test rows | 2680 | 3838 | 0.954 | - | 0.732 | 0.047 | 0.813 |
| pe_weights | target rows | 1725 | 1497 | 0.857 | 0.503 | 0.710 | 0.181 |  |
| pe_weights | hard: Q_deep vs N_0.8/0.95 | 312 | 445 | 0.664 | 0.224 | 0.510 | 0.297 |  |
| pe_weights | twin genera | 861 | 742 | 0.822 | 0.439 | 0.653 | 0.199 |  |
| pe_weights | no allele genomes | 573 | 498 | 0.829 | 0.475 | 0.675 | 0.187 |  |
| pe_weights | recombination high | 575 | 495 | 0.884 | 0.543 | 0.769 | 0.210 |  |
| pe_weights | all test rows | 2685 | 3651 | 0.958 | - | 0.813 | 0.079 | 0.847 |
| pe_without weights | target rows | 1725 | 1497 | 0.925 | 0.749 | 0.833 | 0.141 |  |
| pe_without weights | hard: Q_deep vs N_0.8/0.95 | 312 | 445 | 0.802 | 0.478 | 0.676 | 0.202 |  |
| pe_without weights | twin genera | 861 | 742 | 0.927 | 0.743 | 0.861 | 0.177 |  |
| pe_without weights | no allele genomes | 573 | 498 | 0.832 | 0.476 | 0.702 | 0.213 |  |
| pe_without weights | recombination high | 575 | 495 | 0.942 | 0.767 | 0.875 | 0.162 |  |
| pe_without weights | all test rows | 2685 | 3651 | 0.978 | - | 0.893 | 0.059 | 0.905 |

The models' calls on the target rows by bin: 'called' the share called at the model's knob (Q: sensitivity; N: false-positive rate), 'at 5% FPR' the share above the score 5% of the model's N rows exceed.
| set | class | by | level | rows | base called | ancestry called | alleles called | default called | base at 5% FPR | ancestry at 5% FPR | alleles at 5% FPR | default at 5% FPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 279 | 0.638 | 0.624 | 0.774 | 0.778 | 0.444 | 0.444 | 0.649 | 0.670 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.391 | 0.391 | 0.312 | 0.266 | 0.141 | 0.094 | 0.141 | 0.141 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 573 | 0.670 | 0.691 | 0.696 | 0.682 | 0.440 | 0.449 | 0.513 | 0.546 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 809 | 0.791 | 0.796 | 0.984 | 0.983 | 0.598 | 0.600 | 0.959 | 0.964 |
| pe | Q | depth | 10 | 435 | 0.708 | 0.715 | 0.830 | 0.832 | 0.501 | 0.494 | 0.745 | 0.756 |
| pe | Q | depth | 100 | 426 | 0.746 | 0.751 | 0.845 | 0.843 | 0.563 | 0.554 | 0.763 | 0.791 |
| pe | Q | depth | 3 | 440 | 0.661 | 0.668 | 0.782 | 0.770 | 0.400 | 0.423 | 0.655 | 0.675 |
| pe | Q | depth | 30 | 424 | 0.731 | 0.741 | 0.863 | 0.849 | 0.554 | 0.554 | 0.762 | 0.769 |
| pe | Q | distance to rep | 0.005-0.015 | 773 | 0.717 | 0.719 | 0.894 | 0.900 | 0.473 | 0.462 | 0.777 | 0.805 |
| pe | Q | distance to rep | 0.015-0.03 | 473 | 0.641 | 0.660 | 0.784 | 0.753 | 0.359 | 0.368 | 0.628 | 0.643 |
| pe | Q | distance to rep | <0.005 | 321 | 0.997 | 0.997 | 0.994 | 0.994 | 0.972 | 0.981 | 0.984 | 0.988 |
| pe | Q | distance to rep | >=0.03 | 158 | 0.316 | 0.323 | 0.316 | 0.310 | 0.133 | 0.165 | 0.291 | 0.291 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 567 | 0.765 | 0.771 | 0.877 | 0.866 | 0.533 | 0.547 | 0.760 | 0.792 |
| pe | Q | genes with a congener's segment | <=0.1 | 464 | 0.634 | 0.638 | 0.705 | 0.703 | 0.496 | 0.481 | 0.625 | 0.642 |
| pe | Q | genes with a congener's segment | >0.3 | 257 | 0.732 | 0.759 | 0.868 | 0.864 | 0.490 | 0.502 | 0.774 | 0.786 |
| pe | Q | genes with a congener's segment | none | 437 | 0.712 | 0.712 | 0.879 | 0.872 | 0.483 | 0.481 | 0.778 | 0.778 |
| pe | Q | regime | high | 575 | 0.774 | 0.786 | 0.856 | 0.838 | 0.537 | 0.537 | 0.744 | 0.770 |
| pe | Q | regime | low | 574 | 0.652 | 0.662 | 0.744 | 0.746 | 0.502 | 0.503 | 0.666 | 0.679 |
| pe | Q | regime | none | 576 | 0.708 | 0.707 | 0.889 | 0.885 | 0.472 | 0.476 | 0.781 | 0.792 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 949 | 0.750 | 0.752 | 0.969 | 0.969 | 0.555 | 0.556 | 0.936 | 0.944 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 11 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 582 | 0.674 | 0.691 | 0.687 | 0.675 | 0.431 | 0.436 | 0.500 | 0.531 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 183 | 0.612 | 0.612 | 0.546 | 0.525 | 0.437 | 0.432 | 0.383 | 0.399 |
| pe | Q | role | Q_deep | 312 | 0.519 | 0.526 | 0.657 | 0.654 | 0.247 | 0.244 | 0.490 | 0.510 |
| pe | Q | role | Q_ils | 200 | 0.460 | 0.445 | 0.615 | 0.615 | 0.120 | 0.130 | 0.420 | 0.470 |
| pe | Q | role | Q_imp | 401 | 0.731 | 0.751 | 0.883 | 0.883 | 0.506 | 0.529 | 0.786 | 0.808 |
| pe | Q | role | Q_lone | 277 | 0.747 | 0.758 | 0.805 | 0.773 | 0.560 | 0.538 | 0.693 | 0.693 |
| pe | Q | role | Q_near | 334 | 0.823 | 0.832 | 0.991 | 0.985 | 0.662 | 0.656 | 0.973 | 0.979 |
| pe | Q | role | Q_rep | 201 | 0.985 | 0.980 | 0.970 | 0.975 | 0.940 | 0.945 | 0.950 | 0.960 |
| pe | Q | twin genus | 0 | 864 | 0.773 | 0.780 | 0.808 | 0.782 | 0.561 | 0.567 | 0.704 | 0.709 |
| pe | Q | twin genus | 1 | 861 | 0.649 | 0.656 | 0.851 | 0.864 | 0.446 | 0.444 | 0.757 | 0.785 |
| pe | N | depth | 10 | 380 | 0.192 | 0.200 | 0.161 | 0.145 | 0.068 | 0.063 | 0.063 | 0.068 |
| pe | N | depth | 100 | 377 | 0.186 | 0.175 | 0.117 | 0.119 | 0.045 | 0.056 | 0.042 | 0.037 |
| pe | N | depth | 3 | 359 | 0.167 | 0.178 | 0.142 | 0.123 | 0.042 | 0.036 | 0.053 | 0.050 |
| pe | N | depth | 30 | 381 | 0.160 | 0.173 | 0.136 | 0.118 | 0.045 | 0.045 | 0.042 | 0.045 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 875 | 0.187 | 0.202 | 0.144 | 0.136 | 0.064 | 0.066 | 0.047 | 0.058 |
| pe | N | genes with a segment of S or a congener | >0.3 | 255 | 0.071 | 0.059 | 0.075 | 0.059 | 0.008 | 0.008 | 0.027 | 0.008 |
| pe | N | genes with a segment of S or a congener | none | 367 | 0.223 | 0.218 | 0.172 | 0.150 | 0.046 | 0.041 | 0.074 | 0.060 |
| pe | N | regime | high | 495 | 0.188 | 0.206 | 0.156 | 0.135 | 0.040 | 0.046 | 0.048 | 0.053 |
| pe | N | regime | low | 500 | 0.146 | 0.152 | 0.108 | 0.110 | 0.068 | 0.064 | 0.038 | 0.048 |
| pe | N | regime | none | 502 | 0.195 | 0.187 | 0.153 | 0.133 | 0.042 | 0.040 | 0.064 | 0.050 |
| pe | N | role | N_0 | 195 | 0.000 | 0.000 | 0.005 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 206 | 0.029 | 0.034 | 0.029 | 0.034 | 0.000 | 0.005 | 0.005 | 0.005 |
| pe | N | role | N_0.8 | 219 | 0.178 | 0.237 | 0.142 | 0.132 | 0.009 | 0.014 | 0.023 | 0.037 |
| pe | N | role | N_0.95 | 226 | 0.381 | 0.363 | 0.270 | 0.235 | 0.115 | 0.111 | 0.106 | 0.097 |
| pe | N | role | N_1 | 227 | 0.454 | 0.458 | 0.348 | 0.335 | 0.189 | 0.181 | 0.159 | 0.176 |
| pe | N | role | N_imp0.2 | 213 | 0.056 | 0.061 | 0.056 | 0.047 | 0.009 | 0.014 | 0.014 | 0.014 |
| pe | N | role | N_imp0.4 | 211 | 0.085 | 0.066 | 0.085 | 0.066 | 0.009 | 0.009 | 0.028 | 0.005 |
| pe | N | twin genus | 0 | 755 | 0.164 | 0.166 | 0.110 | 0.087 | 0.046 | 0.048 | 0.040 | 0.032 |
| pe | N | twin genus | 1 | 742 | 0.189 | 0.198 | 0.168 | 0.166 | 0.054 | 0.053 | 0.061 | 0.069 |
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
