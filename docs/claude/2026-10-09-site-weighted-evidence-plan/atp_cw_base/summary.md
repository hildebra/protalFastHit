# Ancestry true-positive test: /home/falk/atp_cw_base

Worlds: test (seed 11), train (seed 23); 36 test genera. See docs/claude/2026-10-09-ancestry-true-positive-test/README.md for the design and the predictions.

## Verdicts

| hypothesis | prediction | outcome | verdict |
|---|---|---|---|
| H1a table | stored alleles are the allele genomes' differences (>= 0.95 exact) | min share exact 1.000 | PASS |
| H3 sites | consensus sites = stem + representative's lineage (no recombination, >= 3 congeners): >= 0.98 of them, recall >= 0.95 | share 0.930, recall 0.866 | FAIL |
| H3 fallback | with 2 congeners the nearest congener's own derived states enter (share in history clearly below 1) | mean share 0.502 | PASS |
| H7 novel reads | the allele scores draw some of a novel species' reads onto the target (it carries the alleles' ancestral base at the representative's own sites) | 6386 of 116107 reads (0.055) | INFO |
| H7 settling | the strains' reads the allele scores move go to their species (>= 0.9) | 1.000 of 3285 moved | PASS |
| H1b/H2b oracle | protal's ancestry, allele and polymorphic features equal the oracle's (error-free reads, every row within 1e-6) | worst share within 1e-6: 1.000 | PASS |
| H1 leak | an allele genome's own reads: allele_explained_share >= 0.95, polymorphic_novel_share <= 0.02 | 0.956, 0.000 | PASS |
| H2 calibration | planted A_a: ancestry_agreement = a (+-0.03) | worst |diff| 0.015 | PASS |
| H1 calibration | planted P_f: allele_explained_share rises with f (a step near 0.5: a read is explained once it shares more than half of the allele's edits in its span) | f 0: 0.000, f 0.25: 0.051, f 0.5: 0.286, f 0.75: 0.542, f 1: 0.705 | PASS |
| H4 ancestry, clean strata | the same with >= 3 congeners and no recombination (the consensus applies) | b 0: 0.999, b 0.5: 0.857 | INFO |
| H5 fixed sites | deep strains with their representative-lineage sites covered (>= 0.8) against late novel species (b 0.8, 0.95), genera with alleles: ancestry_fixed_gain adds to ancestry_agreement (CI above 0) | AUC 0.653 -> 0.885, gain 0.232 [0.154, 0.312]; ancestry_fixed_agreement alone 0.839, the true fixed-site agreement alone 0.839 (160 Q, 297 N rows) | PASS |
| H4 ancestry | ancestry_agreement separates deep strains from novel species that left the stem early (b <= 0.5): AUC >= 0.9 | b 0: 0.906, b 0.5: 0.815 | FAIL |
| H6 confound | allele_explained_share favours strains with a stored relative (Q_near) over late twins, but not lone strains | Q_near 0.909, Q_lone 0.319 | PASS |
| H9 decoy | the decoy (a tip of the deep strains' clade labelled a species) against Q_deep: no feature separates (AUC 0.5 +- 0.1) | largest |AUC - 0.5|: top_identity 0.472 (310 Q, 200 decoy rows) | PASS |
| H9b crown lineage | N_1, a novel lineage leaving at S's crown, against Q_deep: the ancestry features cannot tell them apart; the allele features can where Q_deep shares its side of the crown with allele genomes | ancestry_agreement 0.529, allele_explained_share 0.640 | INFO |
| H4/H5 gain | the three groups add separation over identity and depth (pe, CI above 0) | 0.116 [0.069, 0.164] | PASS |
| H5 ceiling | the true fixed-site agreement alone (the oracle's) adds separation; compare with '+all three' for the gap | 0.000 [-0.006, 0.006] | INFO |
| H5 ceiling, genera with alleles | the true fixed-site agreement (the oracle's) over identity and depth, where the fixed sites exist; compare with '+all three' | 0.074 [0.047, 0.096] | INFO |
| H8 calls pe target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.063 [0.035, 0.093] | PASS |
| H8 calls pe hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.132 [0.066, 0.197] | PASS |
| H8 calls hifi target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.056 [0.034, 0.081] | PASS |
| H8 calls hifi hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.148 [0.070, 0.219] | PASS |
| H9 shuffled | shuffled alleles cost the default model separation (the alleles carry information) | AUC 0.921 -> 0.830 | PASS |

## L0: the database's tables against the trees

Per test genus: copies of the target's genes with stored alleles, alleles per copy, and the share of stored alleles whose edits are exactly an allele genome's differences in the allele's range.
| genus | congeners | twin | allele_genomes | regime | copies | copies_with_alleles | alleles_per_copy | share_exact |
|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 0 | 0 | none | 120 | 0 | 0.000 | - |
| 2 | 2 | 0 | 0 | low | 120 | 0 | 0.000 | - |
| 3 | 2 | 0 | 0 | high | 120 | 0 | 0.000 | - |
| 4 | 2 | 0 | 2 | none | 120 | 120 | 1.992 | 1.000 |
| 5 | 2 | 0 | 2 | low | 120 | 119 | 1.942 | 1.000 |
| 6 | 2 | 0 | 2 | high | 120 | 116 | 1.858 | 1.000 |
| 7 | 2 | 0 | 8 | none | 120 | 117 | 3.258 | 1.000 |
| 8 | 2 | 0 | 8 | low | 120 | 119 | 3.950 | 1.000 |
| 9 | 2 | 0 | 8 | high | 120 | 117 | 3.825 | 1.000 |
| 10 | 2 | 1 | 0 | none | 120 | 0 | 0.000 | - |
| 11 | 2 | 1 | 0 | low | 120 | 0 | 0.000 | - |
| 12 | 2 | 1 | 0 | high | 120 | 0 | 0.000 | - |
| 13 | 2 | 1 | 2 | none | 120 | 117 | 1.892 | 1.000 |
| 14 | 2 | 1 | 2 | low | 120 | 120 | 1.933 | 1.000 |
| 15 | 2 | 1 | 2 | high | 120 | 116 | 1.825 | 1.000 |
| 16 | 2 | 1 | 8 | none | 120 | 118 | 3.800 | 1.000 |
| 17 | 2 | 1 | 8 | low | 120 | 119 | 3.775 | 1.000 |
| 18 | 2 | 1 | 8 | high | 120 | 119 | 3.850 | 1.000 |
| 19 | 6 | 0 | 0 | none | 120 | 0 | 0.000 | - |
| 20 | 6 | 0 | 0 | low | 120 | 0 | 0.000 | - |
| 21 | 6 | 0 | 0 | high | 120 | 0 | 0.000 | - |
| 22 | 6 | 0 | 2 | none | 120 | 120 | 1.792 | 1.000 |
| 23 | 6 | 0 | 2 | low | 120 | 116 | 1.875 | 1.000 |
| 24 | 6 | 0 | 2 | high | 120 | 119 | 1.858 | 1.000 |
| 25 | 6 | 0 | 8 | none | 120 | 120 | 4.000 | 1.000 |
| 26 | 6 | 0 | 8 | low | 120 | 118 | 3.900 | 1.000 |
| 27 | 6 | 0 | 8 | high | 120 | 115 | 3.742 | 1.000 |
| 28 | 6 | 1 | 0 | none | 120 | 0 | 0.000 | - |
| 29 | 6 | 1 | 0 | low | 120 | 0 | 0.000 | - |
| 30 | 6 | 1 | 0 | high | 120 | 0 | 0.000 | - |
| 31 | 6 | 1 | 2 | none | 120 | 118 | 1.925 | 1.000 |
| 32 | 6 | 1 | 2 | low | 120 | 103 | 1.233 | 1.000 |
| 33 | 6 | 1 | 2 | high | 120 | 112 | 1.608 | 1.000 |
| 34 | 6 | 1 | 8 | none | 120 | 84 | 1.300 | 1.000 |
| 35 | 6 | 1 | 8 | low | 120 | 119 | 3.725 | 1.000 |
| 36 | 6 | 1 | 8 | high | 120 | 115 | 2.808 | 1.000 |

Per strain with an allele genome in its species: the share of genes in which its nearest allele genome is among the stored alleles (K = 4, chosen for the sample's coverage since 2026-10-09 evening, farthest first before; the build rejects alleles as far as the nearest congener's copy).
| genome | genus | role | allele_genomes | twin | mrca_nearest_allele | share_genes_nearest_stored |
|---|---|---|---|---|---|---|
| GCA_998004004.1 | 4 | Q_deep | 2 | 0 | 0.000919 | 1.000 |
| GCA_998004005.1 | 4 | Q_deep | 2 | 0 | 0.000919 | 1.000 |
| GCA_998004017.1 | 4 | Q_ils | 2 | 0 | 0.000919 | 1.000 |
| GCA_998004013.1 | 4 | Q_imp0.2 | 2 | 0 | 0.000465 | 1.000 |
| GCA_998004015.1 | 4 | Q_imp0.4 | 2 | 0 | 0.000465 | 1.000 |
| GCA_998004003.1 | 4 | Q_lone | 2 | 0 | 0.003716 | 0.992 |
| GCA_998004001.1 | 4 | Q_near | 2 | 0 | 0.000465 | 1.000 |
| GCA_998004002.1 | 4 | Q_near | 2 | 0 | 0.000678 | 0.992 |
| GCA_998004006.1 | 4 | Q_rep | 2 | 0 | 0.003716 | 0.992 |
| GCA_998005004.1 | 5 | Q_deep | 2 | 0 | 0.004292 | 0.983 |
| GCA_998005005.1 | 5 | Q_deep | 2 | 0 | 0.004292 | 0.983 |
| GCA_998005017.1 | 5 | Q_ils | 2 | 0 | 0.004292 | 0.983 |
| GCA_998005013.1 | 5 | Q_imp0.2 | 2 | 0 | 0.000542 | 0.983 |
| GCA_998005015.1 | 5 | Q_imp0.4 | 2 | 0 | 0.000542 | 0.983 |
| GCA_998005003.1 | 5 | Q_lone | 2 | 0 | 0.007043 | 0.958 |
| GCA_998005001.1 | 5 | Q_near | 2 | 0 | 0.000542 | 0.983 |
| GCA_998005002.1 | 5 | Q_near | 2 | 0 | 9.2e-05 | 0.958 |
| GCA_998005006.1 | 5 | Q_rep | 2 | 0 | 0.000373 | 0.958 |
| GCA_998006004.1 | 6 | Q_deep | 2 | 0 | 0.008763 | 0.917 |
| GCA_998006005.1 | 6 | Q_deep | 2 | 0 | 0.008763 | 0.917 |
| GCA_998006017.1 | 6 | Q_ils | 2 | 0 | 0.008763 | 0.917 |
| GCA_998006013.1 | 6 | Q_imp0.2 | 2 | 0 | 0.001551 | 0.942 |
| GCA_998006015.1 | 6 | Q_imp0.4 | 2 | 0 | 0.001551 | 0.942 |
| GCA_998006003.1 | 6 | Q_lone | 2 | 0 | 0.004133 | 0.942 |
| GCA_998006001.1 | 6 | Q_near | 2 | 0 | 0.001551 | 0.942 |
| GCA_998006002.1 | 6 | Q_near | 2 | 0 | 4.3e-05 | 0.917 |
| GCA_998006006.1 | 6 | Q_rep | 2 | 0 | 0.004133 | 0.942 |
| GCA_998007004.1 | 7 | Q_deep | 8 | 0 | 0.011217 | 0.692 |
| GCA_998007005.1 | 7 | Q_deep | 8 | 0 | 0.011217 | 0.692 |
| GCA_998007017.1 | 7 | Q_ils | 8 | 0 | 0.011217 | 0.692 |
| GCA_998007013.1 | 7 | Q_imp0.2 | 8 | 0 | 0.000103 | 0.692 |
| GCA_998007015.1 | 7 | Q_imp0.4 | 8 | 0 | 0.000103 | 0.692 |
| GCA_998007001.1 | 7 | Q_near | 8 | 0 | 0.000103 | 0.692 |
| GCA_998007002.1 | 7 | Q_near | 8 | 0 | 2.3e-05 | 0.500 |
| GCA_998007003.1 | 7 | Q_near | 8 | 0 | 6.1e-05 | 0.550 |
| GCA_998007006.1 | 7 | Q_rep | 8 | 0 | 0.000784 | 0.550 |
| GCA_998008004.1 | 8 | Q_deep | 8 | 0 | 0.008927 | 0.625 |
| GCA_998008005.1 | 8 | Q_deep | 8 | 0 | 0.008927 | 0.625 |
| GCA_998008017.1 | 8 | Q_ils | 8 | 0 | 0.008927 | 0.625 |
| GCA_998008013.1 | 8 | Q_imp0.2 | 8 | 0 | 8e-06 | 0.625 |
(171 more rows in l0_nearest_allele_stored.tsv)

The build's line: Strain alleles: 7640 alleles (106298 edits) of 2776 gene copies of 24 species, up to 4 each (of 38880 full-reference copies of the database's species and genes: 0 of genomes outside --allele_genome_share 1, 26017 identical to the representative's, 2245 repeated, 45 covering less than 0.5 of it or more than 0.1 apart, 361 as far as the nearest congener's copy or farther; 10212 alleles offered): /home/falk/atp_cw_base/test/db/strain_alleles.tsv

The ancestry sites (the oracle's port of protal's rule, from the database's references; the oracle equals protal, L2a) against the trees: the share on S's stem or the representative's lineage, the recall of each, and why history sites are missed or other sites taken (by the congeners compared at the position).
| genus | congeners | twin | regime | sites | share_in_history | stem_recall | rep_lineage_recall | missed: < 3 congeners, the nearest agrees | missed: a congener carries another base | missed: not compared | other: consensus | other: the nearest's own (< 3 congeners) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 0 | none | 6409 | 0.504 | 0.961 | 0.965 | 51 |  | 67 |  | 3179 |
| 2 | 2 | 0 | low | 8537 | 0.510 | 0.954 | 0.950 | 102 |  | 107 |  | 4186 |
| 3 | 2 | 0 | high | 6939 | 0.477 | 0.912 | 0.856 | 272 |  | 60 |  | 3629 |
| 4 | 2 | 0 | none | 8094 | 0.499 | 0.951 | 0.959 | 73 |  | 114 |  | 4056 |
| 5 | 2 | 0 | low | 6749 | 0.482 | 0.928 | 0.947 | 123 |  | 79 |  | 3493 |
| 6 | 2 | 0 | high | 8316 | 0.501 | 0.888 | 0.874 | 364 |  | 121 |  | 4149 |
| 7 | 2 | 0 | none | 7005 | 0.499 | 0.965 | 0.962 | 46 |  | 77 |  | 3513 |
| 8 | 2 | 0 | low | 8609 | 0.501 | 0.926 | 0.899 | 242 |  | 114 |  | 4297 |
| 9 | 2 | 0 | high | 8919 | 0.496 | 0.874 | 0.871 | 426 |  | 124 |  | 4494 |
| 10 | 2 | 1 | none | 3189 | 0.497 | 0.974 | 0.968 | 12 |  | 35 |  | 1605 |
| 11 | 2 | 1 | low | 3117 | 0.508 | 0.967 | 0.907 | 72 |  | 26 |  | 1535 |
| 12 | 2 | 1 | high | 3119 | 0.500 | 0.910 | 0.816 | 188 |  | 39 |  | 1561 |
| 13 | 2 | 1 | none | 2812 | 0.511 | 0.963 | 0.968 | 16 |  | 32 |  | 1376 |
| 14 | 2 | 1 | low | 3462 | 0.510 | 0.975 | 0.958 | 17 |  | 35 |  | 1698 |
| 15 | 2 | 1 | high | 3248 | 0.518 | 0.901 | 0.939 | 83 |  | 32 |  | 1565 |
| 16 | 2 | 1 | none | 3149 | 0.506 | 0.976 | 0.975 | 10 |  | 29 |  | 1555 |
| 17 | 2 | 1 | low | 2844 | 0.500 | 0.957 | 0.942 | 45 |  | 24 |  | 1423 |
| 18 | 2 | 1 | high | 3092 | 0.495 | 0.919 | 0.960 | 70 |  | 15 |  | 1563 |
| 19 | 6 | 0 | none | 4272 | 0.958 | 0.895 | 0.891 | 10 | 376 | 85 | 161 | 20 |
| 20 | 6 | 0 | low | 4799 | 0.922 | 0.868 | 0.871 | 21 | 503 | 110 | 347 | 29 |
| 21 | 6 | 0 | high | 3717 | 0.959 | 0.763 | 0.717 | 26 | 971 | 98 | 141 | 13 |
| 22 | 6 | 0 | none | 3183 | 0.941 | 0.901 | 0.909 |  | 239 | 67 | 168 | 20 |
| 23 | 6 | 0 | low | 4569 | 0.928 | 0.875 | 0.862 | 7 | 485 | 115 | 304 | 23 |
| 24 | 6 | 0 | high | 3726 | 0.932 | 0.782 | 0.711 | 18 | 929 | 79 | 237 | 15 |
| 25 | 6 | 0 | none | 4912 | 0.956 | 0.916 | 0.907 | 22 | 292 | 121 | 199 | 16 |
| 26 | 6 | 0 | low | 4717 | 0.941 | 0.869 | 0.841 | 11 | 555 | 109 | 262 | 17 |
| 27 | 6 | 0 | high | 4048 | 0.913 | 0.766 | 0.690 | 23 | 1067 | 107 | 346 | 8 |
| 28 | 6 | 1 | none | 1617 | 0.932 | 0.883 | 0.893 |  | 162 | 26 | 94 | 16 |
| 29 | 6 | 1 | low | 1602 | 0.918 | 0.889 | 0.878 | 1 | 168 | 23 | 121 | 11 |
| 30 | 6 | 1 | high | 1568 | 0.844 | 0.810 | 0.536 | 1 | 557 | 23 | 235 | 10 |
| 31 | 6 | 1 | none | 1436 | 0.958 | 0.870 | 0.866 | 1 | 179 | 28 | 58 | 2 |
| 32 | 6 | 1 | low | 1520 | 0.916 | 0.907 | 0.887 | 1 | 139 | 23 | 117 | 11 |
| 33 | 6 | 1 | high | 1830 | 0.808 | 0.816 | 0.723 | 8 | 367 | 37 | 331 | 20 |
| 34 | 6 | 1 | none | 1534 | 0.930 | 0.866 | 0.868 | 1 | 175 | 40 | 96 | 11 |
| 35 | 6 | 1 | low | 1750 | 0.874 | 0.856 | 0.882 |  | 203 | 35 | 202 | 19 |
| 36 | 6 | 1 | high | 1356 | 0.902 | 0.785 | 0.498 | 4 | 622 | 22 | 121 | 12 |

## L1: the focus genomes' reads, with and without the allele scores

The paired-end focus reads (best record of each mate): the share on the target without and with the allele scores, and the reads the scores moved (Q: the target is right; N: no species is).
| class | twin | reads | on_target_noallele | on_target_default | moved | moved_to_target | share_moves_to_target |
|---|---|---|---|---|---|---|---|
| Q | 0 | 66822 | 0.944 | 0.952 | 537 | 537 | 1.000 |
| Q | 1 | 65779 | 0.888 | 0.930 | 2748 | 2748 | 1.000 |
| N | 0 | 57409 | 0.846 | 0.876 | 1763 | 1763 | 1.000 |
| N | 1 | 58698 | 0.757 | 0.836 | 4623 | 4623 | 1.000 |
| class | twin | from | to | reads |
|---|---|---|---|---|
| D | 0 | congener | target | 8 |
| D | 0 | sister | target | 101 |
| D | 1 | congener | target | 10 |
| D | 1 | sister | target | 458 |
| N | 0 | congener | target | 168 |
| N | 0 | sister | target | 1595 |
| N | 1 | congener | target | 90 |
| N | 1 | sister | target | 4533 |
| Q | 0 | congener | target | 84 |
| Q | 0 | sister | target | 453 |
| Q | 1 | congener | target | 110 |
| Q | 1 | sister | target | 2638 |
The protal log's 'strain alleles:' lines (every sample of the default run): 329760 unsure reads took the shifts, 27463 of them to another species.

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
| default | pe | ancestry_sites_per_record | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_agreement | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_congener_share | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_indel_sites_per_record | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_indel_congener_share | 3418 | 0.000 | 1.000 |
| default | pe | allele_explained_share | 3418 | 0.000 | 1.000 |
| default | pe | allele_identity_gain | 3418 | 0.000 | 1.000 |
| default | pe | polymorphic_known_share | 3418 | 0.000 | 1.000 |
| default | pe | polymorphic_novel_share | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_gain | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_agreement | 3418 | 0.000 | 1.000 |
| default | pe | conserved_mismatch_ratio | 3418 | 0.000 | 1.000 |
| default | pe | conserved_mismatch_rate | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_agreement_weighted | 3418 | 0.000 | 1.000 |
| default | pe | nonsynonymous_share | 3418 | 0.000 | 1.000 |
| default | pe | nonsynonymous_conserved_rate | 3418 | 0.000 | 1.000 |
| default | pe | allele_copy_share | 3418 | 0.000 | 1.000 |
| default | pe | allele_sites_per_kb | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_share | 3418 | 0.000 | 1.000 |
| default | pe | column_weight_coverage | 3418 | 0.000 | 1.000 |
(40 more rows in l2_oracle.tsv)

## L2b: the planted genomes and the leaky control (error-free reads)

Mean features of the planted genomes (P_f: a share f of the farthest allele genome's differences; A_a: the congeners' base at a share 1 - a of the consensus sites) and of Q_leak (an allele genome itself). truth_agreement: the genome's agreement at the consensus sites by position.
| role | rows | truth_agreement | ancestry_agreement | allele_explained_share | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | identity |
|---|---|---|---|---|---|---|---|---|
| A_a0 | 7 | 0.001 | 0.002 | 0.147 | 0.310 | 0.026 | 0.001 | 0.976 |
| A_a0.25 | 7 | 0.252 | 0.265 | 0.103 | 0.251 | 0.017 | 0.001 | 0.978 |
| A_a0.5 | 6 | 0.500 | 0.501 | 0.057 | 0.169 | 0.015 | 0.000 | 0.980 |
| A_a0.75 | 7 | 0.751 | 0.748 | 0.022 | 0.084 | 0.008 | -0.001 | 0.980 |
| A_a1 | 6 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0 | 7 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0.25 | 6 | 0.865 | 0.925 | 0.051 | 0.188 | 0.000 | 0.072 | 0.981 |
| P_f0.5 | 7 | 0.943 | 0.839 | 0.286 | 0.432 | 0.000 | 0.161 | 0.980 |
| P_f0.75 | 6 | 0.828 | 0.784 | 0.542 | 0.595 | 0.000 | 0.215 | 0.979 |
| P_f1 | 6 | 0.479 | 0.722 | 0.705 | 0.805 | 0.000 | 0.278 | 0.978 |
| Q_leak | 110 | 0.412 | 0.804 | 0.956 | 0.523 | 0.000 | 0.190 | 0.989 |

## L2c: separation, pe (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1724 | 1494 | 0.780 | 0.686 | 0.799 | 0.566 | 0.498 | 0.815 | 0.186 | 0.622 | 0.549 | 0.537 | 0.426 | 0.486 | 0.653 | 0.610 | 0.500 | 0.815 | 0.501 | 0.500 | 0.822 | 0.651 |
| Q_deep vs N_0 | 310 | 196 | 0.796 | 0.782 | 0.873 | 0.559 | 0.548 | 0.906 | 0.092 | 0.707 | 0.644 | 0.640 | 0.466 | 0.683 | 0.754 | 0.521 | 0.500 | 0.906 | 0.485 | 0.500 | 0.928 | 0.742 |
| Q_deep vs N_0.5 | 310 | 203 | 0.715 | 0.678 | 0.804 | 0.566 | 0.508 | 0.815 | 0.186 | 0.687 | 0.632 | 0.638 | 0.444 | 0.593 | 0.709 | 0.548 | 0.500 | 0.815 | 0.465 | 0.500 | 0.835 | 0.709 |
| Q_deep vs N_0.8 | 310 | 215 | 0.613 | 0.614 | 0.674 | 0.554 | 0.507 | 0.683 | 0.323 | 0.658 | 0.624 | 0.626 | 0.486 | 0.552 | 0.656 | 0.518 | 0.500 | 0.683 | 0.488 | 0.500 | 0.719 | 0.657 |
| Q_deep vs N_0.95 | 310 | 227 | 0.529 | 0.544 | 0.539 | 0.551 | 0.517 | 0.577 | 0.423 | 0.655 | 0.635 | 0.646 | 0.483 | 0.542 | 0.596 | 0.479 | 0.500 | 0.577 | 0.453 | 0.500 | 0.625 | 0.596 |
| Q_deep vs N_1 | 310 | 228 | 0.487 | 0.511 | 0.484 | 0.520 | 0.512 | 0.529 | 0.473 | 0.640 | 0.629 | 0.636 | 0.482 | 0.528 | 0.549 | 0.508 | 0.500 | 0.529 | 0.468 | 0.500 | 0.563 | 0.549 |
| Q_near vs N_0.95 | 334 | 227 | 0.743 | 0.646 | 0.726 | 0.600 | 0.493 | 0.763 | 0.238 | 0.909 | 0.783 | 0.752 | 0.462 | 0.599 | 0.890 | 0.620 | 0.500 | 0.763 | 0.508 | 0.500 | 0.929 | 0.890 |
| Q_lone vs N_0.95 | 277 | 227 | 0.762 | 0.664 | 0.774 | 0.544 | 0.491 | 0.775 | 0.226 | 0.319 | 0.295 | 0.299 | 0.358 | 0.281 | 0.346 | 0.635 | 0.500 | 0.775 | 0.473 | 0.500 | 0.591 | 0.345 |
| Q_rep vs N_0 | 201 | 196 | 0.993 | 0.891 | 0.987 | 0.600 | 0.510 | 1.000 | 0.000 | 0.458 | 0.350 | 0.315 | 0.404 | 0.385 | 0.707 | 0.851 | 0.500 | 1.000 | 0.494 | 0.500 | 1.000 | 0.691 |
| Q_imp vs N (b <= 0.95) | 402 | 841 | 0.818 | 0.739 | 0.842 | 0.572 | 0.500 | 0.857 | 0.144 | 0.666 | 0.586 | 0.566 | 0.400 | 0.498 | 0.675 | 0.585 | 0.500 | 0.857 | 0.510 | 0.500 | 0.860 | 0.671 |
| Q (no import) vs N_imp | 1322 | 425 | 0.801 | 0.661 | 0.823 | 0.574 | 0.499 | 0.844 | 0.156 | 0.612 | 0.531 | 0.524 | 0.439 | 0.492 | 0.675 | 0.619 | 0.500 | 0.844 | 0.502 | 0.500 | 0.871 | 0.675 |
| Q_ils vs N_0.95 | 200 | 227 | 0.529 | 0.537 | 0.510 | 0.526 | 0.500 | 0.489 | 0.511 | 0.585 | 0.570 | 0.592 | 0.474 | 0.487 | 0.465 | 0.399 | 0.500 | 0.489 | 0.492 | 0.500 | 0.480 | 0.465 |
| decoy: Q_deep vs decoy | 310 | 200 | 0.478 | 0.472 | 0.485 | 0.492 | 0.519 | 0.494 | 0.503 | 0.525 | 0.524 | 0.513 | 0.520 | 0.512 | 0.526 | 0.478 | 0.500 | 0.494 | 0.493 | 0.500 | 0.520 | 0.527 |
| allele genomes 0 | 573 | 496 | 0.767 | 0.684 | 0.799 | 0.535 | 0.509 | 0.807 | 0.192 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.606 | 0.500 | 0.807 | 0.493 | 0.500 | 0.807 | 0.500 |
| allele genomes 2 | 576 | 497 | 0.803 | 0.681 | 0.809 | 0.589 | 0.503 | 0.827 | 0.176 | 0.787 | 0.646 | 0.616 | 0.397 | 0.469 | 0.840 | 0.610 | 0.500 | 0.827 | 0.488 | 0.500 | 0.866 | 0.840 |
| allele genomes 8 | 575 | 501 | 0.772 | 0.693 | 0.794 | 0.581 | 0.485 | 0.815 | 0.185 | 0.769 | 0.581 | 0.573 | 0.336 | 0.481 | 0.836 | 0.616 | 0.500 | 0.815 | 0.523 | 0.500 | 0.895 | 0.836 |
| congeners 2 | 862 | 746 | 0.796 | 0.703 | 0.795 | 0.568 | 0.478 | 0.838 | 0.164 | 0.629 | 0.562 | 0.528 | 0.422 | 0.470 | 0.646 | 0.636 | 0.500 | 0.838 | 0.495 | 0.500 | 0.820 | 0.646 |
| congeners 6 | 862 | 748 | 0.770 | 0.670 | 0.802 | 0.569 | 0.523 | 0.850 | 0.152 | 0.617 | 0.539 | 0.544 | 0.431 | 0.503 | 0.666 | 0.612 | 0.500 | 0.850 | 0.507 | 0.500 | 0.846 | 0.662 |
| twin 0 | 864 | 759 | 0.875 | 0.773 | 0.869 | 0.701 | 0.493 | 0.828 | 0.173 | 0.633 | 0.543 | 0.509 | 0.375 | 0.483 | 0.680 | 0.627 | 0.500 | 0.828 | 0.499 | 0.500 | 0.842 | 0.679 |
| twin 1 | 860 | 735 | 0.788 | 0.607 | 0.713 | 0.525 | 0.516 | 0.803 | 0.197 | 0.614 | 0.558 | 0.562 | 0.469 | 0.487 | 0.632 | 0.597 | 0.500 | 0.803 | 0.501 | 0.500 | 0.805 | 0.630 |
| regime low | 574 | 496 | 0.771 | 0.695 | 0.797 | 0.560 | 0.497 | 0.801 | 0.200 | 0.637 | 0.558 | 0.562 | 0.410 | 0.496 | 0.655 | 0.597 | 0.500 | 0.801 | 0.486 | 0.500 | 0.784 | 0.654 |
| regime high | 575 | 492 | 0.796 | 0.695 | 0.833 | 0.572 | 0.501 | 0.855 | 0.143 | 0.647 | 0.539 | 0.492 | 0.438 | 0.460 | 0.654 | 0.609 | 0.500 | 0.855 | 0.512 | 0.500 | 0.840 | 0.652 |
| regime none | 575 | 506 | 0.774 | 0.670 | 0.770 | 0.569 | 0.497 | 0.798 | 0.204 | 0.587 | 0.549 | 0.555 | 0.431 | 0.502 | 0.652 | 0.628 | 0.500 | 0.798 | 0.505 | 0.500 | 0.853 | 0.649 |
| depth 3 | 439 | 354 | 0.759 | 0.754 | 0.755 | 0.556 | 0.474 | 0.787 | 0.216 | 0.614 | 0.547 | 0.528 | 0.463 | 0.481 | 0.633 | 0.574 | 0.500 | 0.787 | 0.477 | 0.500 | 0.785 | 0.630 |
| depth 10 | 435 | 382 | 0.779 | 0.717 | 0.805 | 0.575 | 0.520 | 0.836 | 0.163 | 0.613 | 0.551 | 0.539 | 0.444 | 0.492 | 0.655 | 0.594 | 0.500 | 0.836 | 0.518 | 0.500 | 0.833 | 0.652 |
| depth 30 | 424 | 381 | 0.793 | 0.683 | 0.828 | 0.562 | 0.499 | 0.819 | 0.182 | 0.625 | 0.545 | 0.533 | 0.392 | 0.484 | 0.660 | 0.636 | 0.500 | 0.819 | 0.494 | 0.500 | 0.834 | 0.659 |
| depth 100 | 426 | 377 | 0.790 | 0.632 | 0.815 | 0.573 | 0.501 | 0.820 | 0.179 | 0.640 | 0.554 | 0.548 | 0.359 | 0.492 | 0.677 | 0.654 | 0.500 | 0.820 | 0.522 | 0.500 | 0.844 | 0.677 |
| all, within identity bands of 0.005 | 1724 | 1494 | 0.560 | 0.546 | 0.644 | 0.686 | 0.690 | 0.731 | 0.270 | 0.678 | 0.676 | 0.626 | 0.503 | 0.554 | 0.622 | 0.466 | 0.500 | 0.731 | 0.481 | 0.500 | 0.762 | 0.621 |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 158 | 97 | 0.794 | 0.858 | 0.968 | 0.743 | 0.644 | 0.776 | 0.949 | 0.776 |
| Q_deep vs N_0 | congeners 6 | 152 | 99 | 0.801 | 0.889 | 0.990 | 0.674 | 0.634 | 0.766 | 0.995 | 0.741 |
| Q_deep vs N_0 | no twin | 162 | 100 | 0.988 | 0.979 | 0.907 | 0.765 | 0.663 | 0.806 | 0.952 | 0.796 |
| Q_deep vs N_0 | twin | 148 | 96 | 0.923 | 0.787 | 0.909 | 0.658 | 0.610 | 0.711 | 0.906 | 0.698 |
| Q_deep vs N_0 | no allele genomes | 82 | 65 | 0.821 | 0.855 | 0.866 | 0.500 | 0.500 | 0.500 | 0.866 | 0.500 |
| Q_deep vs N_0 | allele genomes | 228 | 131 | 0.785 | 0.885 | 0.922 | 0.836 | 0.707 | 0.922 | 0.965 | 0.922 |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 22 | 0.735 | 0.843 | 0.999 | 0.673 | 0.711 | 1.000 | 1.000 | 1.000 |
| Q_deep vs N_0 | recombination none | 106 | 69 | 0.773 | 0.865 | 0.868 | 0.669 | 0.658 | 0.770 | 0.940 | 0.754 |
| Q_deep vs N_0 | recombination high | 102 | 60 | 0.792 | 0.875 | 0.964 | 0.762 | 0.622 | 0.776 | 0.980 | 0.766 |
| Q_deep vs N_0.5 | congeners 2 | 158 | 102 | 0.706 | 0.799 | 0.863 | 0.729 | 0.656 | 0.735 | 0.840 | 0.735 |
| Q_deep vs N_0.5 | congeners 6 | 152 | 101 | 0.732 | 0.810 | 0.907 | 0.644 | 0.620 | 0.709 | 0.895 | 0.709 |
| Q_deep vs N_0.5 | no twin | 162 | 102 | 0.920 | 0.904 | 0.815 | 0.728 | 0.660 | 0.776 | 0.888 | 0.776 |
| Q_deep vs N_0.5 | twin | 148 | 101 | 0.817 | 0.717 | 0.809 | 0.654 | 0.613 | 0.656 | 0.782 | 0.656 |
| Q_deep vs N_0.5 | no allele genomes | 82 | 69 | 0.741 | 0.804 | 0.768 | 0.500 | 0.500 | 0.500 | 0.768 | 0.500 |
| Q_deep vs N_0.5 | allele genomes | 228 | 134 | 0.708 | 0.811 | 0.842 | 0.804 | 0.708 | 0.853 | 0.890 | 0.853 |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.583 | 0.705 | 0.857 | 0.633 | 0.723 | 0.866 | 0.899 | 0.868 |
| Q_deep vs N_0.5 | recombination none | 106 | 68 | 0.688 | 0.758 | 0.798 | 0.655 | 0.669 | 0.709 | 0.855 | 0.711 |
| Q_deep vs N_0.5 | recombination high | 102 | 73 | 0.704 | 0.824 | 0.854 | 0.739 | 0.616 | 0.733 | 0.895 | 0.733 |
| Q_deep vs N_0.8 | congeners 2 | 158 | 109 | 0.607 | 0.641 | 0.673 | 0.706 | 0.649 | 0.679 | 0.720 | 0.679 |
| Q_deep vs N_0.8 | congeners 6 | 152 | 106 | 0.630 | 0.703 | 0.755 | 0.610 | 0.600 | 0.643 | 0.741 | 0.645 |
| Q_deep vs N_0.8 | no twin | 162 | 112 | 0.713 | 0.718 | 0.665 | 0.694 | 0.650 | 0.703 | 0.750 | 0.704 |
| Q_deep vs N_0.8 | twin | 148 | 103 | 0.670 | 0.630 | 0.710 | 0.628 | 0.600 | 0.624 | 0.684 | 0.624 |
| Q_deep vs N_0.8 | no allele genomes | 82 | 67 | 0.637 | 0.664 | 0.679 | 0.500 | 0.500 | 0.500 | 0.679 | 0.500 |
| Q_deep vs N_0.8 | allele genomes | 228 | 148 | 0.606 | 0.684 | 0.686 | 0.768 | 0.706 | 0.763 | 0.789 | 0.763 |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 24 | 0.596 | 0.639 | 0.695 | 0.611 | 0.675 | 0.746 | 0.774 | 0.755 |
| Q_deep vs N_0.8 | recombination none | 106 | 74 | 0.603 | 0.588 | 0.659 | 0.637 | 0.652 | 0.657 | 0.737 | 0.658 |
| Q_deep vs N_0.8 | recombination high | 102 | 74 | 0.614 | 0.758 | 0.742 | 0.696 | 0.603 | 0.682 | 0.753 | 0.682 |
| Q_deep vs N_0.95 | congeners 2 | 158 | 112 | 0.522 | 0.534 | 0.557 | 0.695 | 0.655 | 0.608 | 0.647 | 0.608 |
| Q_deep vs N_0.95 | congeners 6 | 152 | 115 | 0.537 | 0.545 | 0.588 | 0.616 | 0.636 | 0.576 | 0.609 | 0.576 |
| Q_deep vs N_0.95 | no twin | 162 | 116 | 0.553 | 0.578 | 0.571 | 0.674 | 0.660 | 0.635 | 0.637 | 0.635 |
| Q_deep vs N_0.95 | twin | 148 | 111 | 0.553 | 0.498 | 0.584 | 0.643 | 0.632 | 0.571 | 0.610 | 0.571 |
| Q_deep vs N_0.95 | no allele genomes | 82 | 78 | 0.523 | 0.511 | 0.526 | 0.500 | 0.500 | 0.500 | 0.526 | 0.500 |
| Q_deep vs N_0.95 | allele genomes | 228 | 149 | 0.543 | 0.558 | 0.604 | 0.738 | 0.721 | 0.618 | 0.668 | 0.618 |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.504 | 0.624 | 0.562 | 0.601 | 0.720 | 0.572 | 0.590 | 0.569 |
| Q_deep vs N_0.95 | recombination none | 106 | 76 | 0.507 | 0.524 | 0.559 | 0.648 | 0.673 | 0.607 | 0.635 | 0.606 |
| Q_deep vs N_0.95 | recombination high | 102 | 73 | 0.535 | 0.546 | 0.601 | 0.684 | 0.639 | 0.621 | 0.662 | 0.622 |

## L2c: separation, se (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1714 | 1491 | 0.758 | 0.676 | 0.769 | 0.577 | 0.504 | 0.784 | 0.217 | 0.620 | 0.553 | 0.539 | 0.439 | 0.487 | 0.645 | 0.579 | 0.500 | 0.784 | 0.495 | 0.500 | 0.804 | - |
| Q_deep vs N_0 | 309 | 199 | 0.762 | 0.752 | 0.842 | 0.578 | 0.561 | 0.874 | 0.140 | 0.688 | 0.634 | 0.631 | 0.490 | 0.652 | 0.746 | 0.470 | 0.500 | 0.874 | 0.523 | 0.500 | 0.907 | - |
| Q_deep vs N_0.5 | 309 | 202 | 0.665 | 0.655 | 0.726 | 0.587 | 0.524 | 0.764 | 0.227 | 0.677 | 0.635 | 0.629 | 0.482 | 0.603 | 0.711 | 0.531 | 0.500 | 0.764 | 0.489 | 0.500 | 0.839 | - |
| Q_deep vs N_0.8 | 309 | 213 | 0.602 | 0.584 | 0.649 | 0.572 | 0.509 | 0.654 | 0.344 | 0.651 | 0.617 | 0.624 | 0.491 | 0.546 | 0.631 | 0.503 | 0.500 | 0.654 | 0.500 | 0.500 | 0.711 | - |
| Q_deep vs N_0.95 | 309 | 226 | 0.524 | 0.532 | 0.530 | 0.559 | 0.512 | 0.561 | 0.435 | 0.638 | 0.621 | 0.623 | 0.492 | 0.532 | 0.584 | 0.484 | 0.500 | 0.561 | 0.512 | 0.500 | 0.609 | - |
| Q_deep vs N_1 | 309 | 227 | 0.464 | 0.486 | 0.465 | 0.539 | 0.520 | 0.502 | 0.490 | 0.642 | 0.632 | 0.638 | 0.485 | 0.550 | 0.559 | 0.484 | 0.500 | 0.502 | 0.474 | 0.500 | 0.564 | - |
| Q_near vs N_0.95 | 333 | 226 | 0.733 | 0.644 | 0.701 | 0.612 | 0.484 | 0.732 | 0.268 | 0.892 | 0.772 | 0.732 | 0.457 | 0.580 | 0.866 | 0.595 | 0.500 | 0.732 | 0.518 | 0.500 | 0.887 | - |
| Q_lone vs N_0.95 | 276 | 226 | 0.760 | 0.663 | 0.748 | 0.548 | 0.487 | 0.737 | 0.259 | 0.320 | 0.295 | 0.301 | 0.390 | 0.284 | 0.328 | 0.591 | 0.500 | 0.737 | 0.484 | 0.500 | 0.555 | - |
| Q_rep vs N_0 | 201 | 199 | 0.984 | 0.883 | 0.983 | 0.599 | 0.533 | 0.995 | 0.020 | 0.462 | 0.363 | 0.320 | 0.445 | 0.373 | 0.700 | 0.755 | 0.500 | 0.995 | 0.522 | 0.500 | 0.988 | - |
| Q_imp vs N (b <= 0.95) | 398 | 840 | 0.793 | 0.721 | 0.814 | 0.582 | 0.507 | 0.835 | 0.169 | 0.658 | 0.585 | 0.565 | 0.409 | 0.493 | 0.666 | 0.582 | 0.500 | 0.835 | 0.500 | 0.500 | 0.849 | - |
| Q (no import) vs N_imp | 1316 | 424 | 0.779 | 0.658 | 0.785 | 0.581 | 0.502 | 0.800 | 0.200 | 0.613 | 0.539 | 0.527 | 0.446 | 0.493 | 0.663 | 0.581 | 0.500 | 0.800 | 0.486 | 0.500 | 0.837 | - |
| Q_ils vs N_0.95 | 197 | 226 | 0.529 | 0.550 | 0.514 | 0.520 | 0.489 | 0.490 | 0.508 | 0.559 | 0.553 | 0.565 | 0.492 | 0.473 | 0.463 | 0.434 | 0.500 | 0.490 | 0.534 | 0.500 | 0.463 | - |
| decoy: Q_deep vs decoy | 309 | 198 | 0.460 | 0.461 | 0.454 | 0.510 | 0.512 | 0.492 | 0.503 | 0.520 | 0.527 | 0.515 | 0.513 | 0.522 | 0.520 | 0.483 | 0.500 | 0.492 | 0.531 | 0.500 | 0.536 | - |
| allele genomes 0 | 570 | 497 | 0.751 | 0.680 | 0.775 | 0.557 | 0.518 | 0.774 | 0.226 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.571 | 0.500 | 0.774 | 0.489 | 0.500 | 0.774 | - |
| allele genomes 2 | 572 | 496 | 0.773 | 0.675 | 0.774 | 0.591 | 0.510 | 0.801 | 0.200 | 0.791 | 0.665 | 0.627 | 0.410 | 0.481 | 0.811 | 0.589 | 0.500 | 0.801 | 0.482 | 0.500 | 0.856 | - |
| allele genomes 8 | 572 | 498 | 0.756 | 0.671 | 0.764 | 0.589 | 0.487 | 0.778 | 0.221 | 0.745 | 0.575 | 0.562 | 0.388 | 0.474 | 0.828 | 0.580 | 0.500 | 0.778 | 0.513 | 0.500 | 0.866 | - |
| congeners 2 | 857 | 747 | 0.771 | 0.685 | 0.767 | 0.587 | 0.479 | 0.802 | 0.201 | 0.628 | 0.569 | 0.533 | 0.441 | 0.476 | 0.638 | 0.603 | 0.500 | 0.802 | 0.497 | 0.500 | 0.804 | - |
| congeners 6 | 857 | 744 | 0.750 | 0.667 | 0.772 | 0.567 | 0.536 | 0.818 | 0.186 | 0.613 | 0.541 | 0.543 | 0.438 | 0.500 | 0.658 | 0.572 | 0.500 | 0.818 | 0.493 | 0.500 | 0.825 | - |
| twin 0 | 858 | 756 | 0.854 | 0.753 | 0.841 | 0.677 | 0.506 | 0.808 | 0.192 | 0.631 | 0.547 | 0.515 | 0.390 | 0.486 | 0.676 | 0.596 | 0.500 | 0.808 | 0.498 | 0.500 | 0.834 | - |
| twin 1 | 856 | 735 | 0.742 | 0.602 | 0.682 | 0.561 | 0.513 | 0.760 | 0.242 | 0.611 | 0.563 | 0.559 | 0.484 | 0.488 | 0.621 | 0.568 | 0.500 | 0.760 | 0.491 | 0.500 | 0.776 | - |
| regime low | 573 | 498 | 0.743 | 0.672 | 0.760 | 0.570 | 0.507 | 0.772 | 0.231 | 0.638 | 0.567 | 0.566 | 0.422 | 0.496 | 0.642 | 0.557 | 0.500 | 0.772 | 0.475 | 0.500 | 0.768 | - |
| regime high | 571 | 489 | 0.776 | 0.678 | 0.802 | 0.580 | 0.503 | 0.821 | 0.177 | 0.635 | 0.540 | 0.496 | 0.449 | 0.461 | 0.650 | 0.590 | 0.500 | 0.821 | 0.523 | 0.500 | 0.823 | - |
| regime none | 570 | 504 | 0.758 | 0.677 | 0.750 | 0.580 | 0.503 | 0.767 | 0.234 | 0.590 | 0.551 | 0.551 | 0.447 | 0.501 | 0.645 | 0.594 | 0.500 | 0.767 | 0.487 | 0.500 | 0.828 | - |
| depth 3 | 429 | 351 | 0.718 | 0.725 | 0.712 | 0.560 | 0.485 | 0.722 | 0.281 | 0.601 | 0.562 | 0.536 | 0.483 | 0.495 | 0.621 | 0.518 | 0.500 | 0.722 | 0.486 | 0.500 | 0.743 | - |
| depth 10 | 435 | 382 | 0.764 | 0.727 | 0.774 | 0.579 | 0.522 | 0.817 | 0.184 | 0.610 | 0.547 | 0.533 | 0.457 | 0.486 | 0.650 | 0.590 | 0.500 | 0.817 | 0.510 | 0.500 | 0.827 | - |
| depth 30 | 424 | 381 | 0.777 | 0.685 | 0.802 | 0.575 | 0.504 | 0.801 | 0.200 | 0.638 | 0.551 | 0.540 | 0.404 | 0.481 | 0.649 | 0.608 | 0.500 | 0.800 | 0.496 | 0.500 | 0.822 | - |
| depth 100 | 426 | 377 | 0.778 | 0.627 | 0.800 | 0.591 | 0.506 | 0.804 | 0.195 | 0.645 | 0.557 | 0.549 | 0.385 | 0.490 | 0.677 | 0.624 | 0.500 | 0.804 | 0.491 | 0.500 | 0.843 | - |
| all, within identity bands of 0.005 | 1714 | 1491 | 0.552 | 0.550 | 0.633 | 0.685 | 0.659 | 0.700 | 0.304 | 0.679 | 0.677 | 0.633 | 0.508 | 0.561 | 0.619 | 0.481 | 0.500 | 0.700 | 0.489 | 0.500 | 0.742 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 158 | 99 | 0.749 | 0.824 | 0.913 | 0.729 | 0.644 | 0.758 | 0.915 | - |
| Q_deep vs N_0 | congeners 6 | 151 | 100 | 0.775 | 0.861 | 0.976 | 0.647 | 0.616 | 0.760 | 0.980 | - |
| Q_deep vs N_0 | no twin | 161 | 99 | 0.974 | 0.939 | 0.890 | 0.745 | 0.664 | 0.803 | 0.948 | - |
| Q_deep vs N_0 | twin | 148 | 100 | 0.844 | 0.779 | 0.861 | 0.639 | 0.594 | 0.701 | 0.870 | - |
| Q_deep vs N_0 | no allele genomes | 82 | 68 | 0.783 | 0.848 | 0.847 | 0.500 | 0.500 | 0.500 | 0.847 | - |
| Q_deep vs N_0 | allele genomes | 227 | 131 | 0.757 | 0.856 | 0.893 | 0.806 | 0.690 | 0.910 | 0.938 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 22 | 0.729 | 0.864 | 1.000 | 0.650 | 0.652 | 1.000 | 1.000 | - |
| Q_deep vs N_0 | recombination none | 105 | 70 | 0.763 | 0.836 | 0.846 | 0.666 | 0.653 | 0.776 | 0.938 | - |
| Q_deep vs N_0 | recombination high | 101 | 61 | 0.735 | 0.826 | 0.908 | 0.720 | 0.618 | 0.764 | 0.937 | - |
| Q_deep vs N_0.5 | congeners 2 | 158 | 100 | 0.646 | 0.703 | 0.754 | 0.725 | 0.652 | 0.726 | 0.839 | - |
| Q_deep vs N_0.5 | congeners 6 | 151 | 102 | 0.691 | 0.747 | 0.871 | 0.629 | 0.603 | 0.719 | 0.894 | - |
| Q_deep vs N_0.5 | no twin | 161 | 101 | 0.858 | 0.829 | 0.779 | 0.724 | 0.668 | 0.777 | 0.886 | - |
| Q_deep vs N_0.5 | twin | 148 | 101 | 0.702 | 0.645 | 0.746 | 0.636 | 0.587 | 0.660 | 0.791 | - |
| Q_deep vs N_0.5 | no allele genomes | 82 | 68 | 0.669 | 0.758 | 0.714 | 0.500 | 0.500 | 0.500 | 0.714 | - |
| Q_deep vs N_0.5 | allele genomes | 227 | 134 | 0.674 | 0.723 | 0.798 | 0.790 | 0.693 | 0.858 | 0.897 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.605 | 0.625 | 0.843 | 0.646 | 0.696 | 0.911 | 0.949 | - |
| Q_deep vs N_0.5 | recombination none | 105 | 66 | 0.653 | 0.688 | 0.758 | 0.651 | 0.650 | 0.729 | 0.860 | - |
| Q_deep vs N_0.5 | recombination high | 101 | 72 | 0.656 | 0.756 | 0.772 | 0.713 | 0.613 | 0.733 | 0.891 | - |
| Q_deep vs N_0.8 | congeners 2 | 158 | 111 | 0.586 | 0.629 | 0.604 | 0.706 | 0.661 | 0.646 | 0.716 | - |
| Q_deep vs N_0.8 | congeners 6 | 151 | 102 | 0.622 | 0.673 | 0.743 | 0.592 | 0.584 | 0.627 | 0.733 | - |
| Q_deep vs N_0.8 | no twin | 161 | 112 | 0.693 | 0.690 | 0.654 | 0.690 | 0.645 | 0.692 | 0.747 | - |
| Q_deep vs N_0.8 | twin | 148 | 101 | 0.613 | 0.607 | 0.657 | 0.618 | 0.601 | 0.587 | 0.672 | - |
| Q_deep vs N_0.8 | no allele genomes | 82 | 66 | 0.628 | 0.694 | 0.657 | 0.500 | 0.500 | 0.500 | 0.657 | - |
| Q_deep vs N_0.8 | allele genomes | 227 | 147 | 0.598 | 0.637 | 0.656 | 0.756 | 0.702 | 0.717 | 0.768 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.601 | 0.585 | 0.687 | 0.596 | 0.633 | 0.683 | 0.733 | - |
| Q_deep vs N_0.8 | recombination none | 105 | 72 | 0.583 | 0.584 | 0.638 | 0.640 | 0.648 | 0.625 | 0.731 | - |
| Q_deep vs N_0.8 | recombination high | 101 | 74 | 0.613 | 0.711 | 0.704 | 0.673 | 0.603 | 0.647 | 0.744 | - |
| Q_deep vs N_0.95 | congeners 2 | 158 | 111 | 0.528 | 0.513 | 0.529 | 0.680 | 0.639 | 0.594 | 0.611 | - |
| Q_deep vs N_0.95 | congeners 6 | 151 | 115 | 0.523 | 0.550 | 0.577 | 0.597 | 0.601 | 0.568 | 0.605 | - |
| Q_deep vs N_0.95 | no twin | 161 | 116 | 0.546 | 0.564 | 0.566 | 0.677 | 0.654 | 0.632 | 0.650 | - |
| Q_deep vs N_0.95 | twin | 148 | 110 | 0.533 | 0.500 | 0.556 | 0.605 | 0.591 | 0.552 | 0.568 | - |
| Q_deep vs N_0.95 | no allele genomes | 82 | 76 | 0.548 | 0.568 | 0.529 | 0.500 | 0.500 | 0.500 | 0.529 | - |
| Q_deep vs N_0.95 | allele genomes | 227 | 150 | 0.524 | 0.522 | 0.579 | 0.725 | 0.686 | 0.601 | 0.640 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.474 | 0.561 | 0.562 | 0.599 | 0.662 | 0.560 | 0.548 | - |
| Q_deep vs N_0.95 | recombination none | 105 | 77 | 0.501 | 0.515 | 0.546 | 0.630 | 0.650 | 0.603 | 0.622 | - |
| Q_deep vs N_0.95 | recombination high | 101 | 71 | 0.518 | 0.521 | 0.572 | 0.671 | 0.618 | 0.590 | 0.627 | - |

## L2c: separation, hifi (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | conserved_mismatch_ratio | conserved_mismatch_rate | ancestry_agreement_weighted | nonsynonymous_share | nonsynonymous_conserved_rate | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 576 | 513 | 0.814 | 0.758 | 0.877 | 0.604 | 0.500 | 0.840 | 0.162 | 0.662 | 0.600 | 0.540 | 0.296 | 0.491 | 0.669 | 0.354 | 0.500 | 0.840 | 0.527 | 0.500 | 0.829 | - |
| Q_deep vs N_0 | 101 | 75 | 0.805 | 0.789 | 0.924 | 0.666 | 0.489 | 0.939 | 0.057 | 0.745 | 0.712 | 0.709 | 0.458 | 0.709 | 0.747 | 0.372 | 0.500 | 0.939 | 0.477 | 0.500 | 0.946 | - |
| Q_deep vs N_0.5 | 101 | 74 | 0.777 | 0.756 | 0.967 | 0.606 | 0.453 | 0.849 | 0.153 | 0.682 | 0.670 | 0.669 | 0.362 | 0.608 | 0.716 | 0.491 | 0.500 | 0.849 | 0.481 | 0.500 | 0.826 | - |
| Q_deep vs N_0.8 | 101 | 65 | 0.657 | 0.715 | 0.791 | 0.599 | 0.488 | 0.708 | 0.297 | 0.673 | 0.671 | 0.686 | 0.387 | 0.572 | 0.672 | 0.524 | 0.500 | 0.708 | 0.509 | 0.500 | 0.724 | - |
| Q_deep vs N_0.95 | 101 | 73 | 0.584 | 0.637 | 0.629 | 0.562 | 0.456 | 0.569 | 0.430 | 0.643 | 0.652 | 0.668 | 0.416 | 0.525 | 0.608 | 0.531 | 0.500 | 0.569 | 0.477 | 0.500 | 0.614 | - |
| Q_deep vs N_1 | 101 | 72 | 0.538 | 0.617 | 0.579 | 0.521 | 0.476 | 0.531 | 0.469 | 0.625 | 0.643 | 0.660 | 0.367 | 0.469 | 0.549 | 0.507 | 0.500 | 0.531 | 0.482 | 0.500 | 0.524 | - |
| Q_near vs N_0.95 | 100 | 73 | 0.812 | 0.833 | 0.805 | 0.635 | 0.484 | 0.785 | 0.217 | 1.000 | 0.883 | 0.717 | 0.280 | 0.592 | 0.976 | 0.465 | 0.500 | 0.785 | 0.536 | 0.500 | 0.982 | - |
| Q_lone vs N_0.95 | 92 | 73 | 0.748 | 0.770 | 0.769 | 0.532 | 0.456 | 0.758 | 0.244 | 0.338 | 0.305 | 0.305 | 0.217 | 0.274 | 0.374 | 0.446 | 0.500 | 0.758 | 0.506 | 0.500 | 0.601 | - |
| Q_rep vs N_0 | 73 | 75 | 0.991 | 0.960 | 0.999 | 0.715 | 0.515 | 1.000 | 0.000 | 0.652 | 0.500 | 0.380 | 0.380 | 0.477 | 0.783 | 0.288 | 0.500 | 1.000 | 0.599 | 0.500 | 1.000 | - |
| Q_imp vs N (b <= 0.95) | 136 | 287 | 0.825 | 0.822 | 0.908 | 0.621 | 0.548 | 0.877 | 0.124 | 0.660 | 0.618 | 0.534 | 0.288 | 0.478 | 0.625 | 0.211 | 0.500 | 0.877 | 0.522 | 0.500 | 0.833 | - |
| Q (no import) vs N_imp | 440 | 154 | 0.857 | 0.657 | 0.931 | 0.615 | 0.482 | 0.877 | 0.125 | 0.665 | 0.587 | 0.539 | 0.285 | 0.511 | 0.714 | 0.388 | 0.500 | 0.877 | 0.533 | 0.500 | 0.893 | - |
| Q_ils vs N_0.95 | 74 | 73 | 0.547 | 0.587 | 0.535 | 0.488 | 0.483 | 0.502 | 0.501 | 0.630 | 0.631 | 0.635 | 0.393 | 0.511 | 0.485 | 0.372 | 0.500 | 0.502 | 0.481 | 0.500 | 0.489 | - |
| decoy: Q_deep vs decoy | 101 | 63 | 0.476 | 0.468 | 0.487 | 0.531 | 0.491 | 0.448 | 0.551 | 0.528 | 0.534 | 0.534 | 0.533 | 0.513 | 0.521 | 0.519 | 0.500 | 0.448 | 0.503 | 0.500 | 0.485 | - |
| allele genomes 0 | 192 | 172 | 0.806 | 0.740 | 0.862 | 0.585 | 0.495 | 0.829 | 0.171 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.336 | 0.500 | 0.829 | 0.496 | 0.500 | 0.829 | - |
| allele genomes 2 | 192 | 171 | 0.824 | 0.753 | 0.887 | 0.596 | 0.504 | 0.848 | 0.153 | 0.869 | 0.775 | 0.660 | 0.139 | 0.496 | 0.870 | 0.368 | 0.500 | 0.848 | 0.541 | 0.500 | 0.874 | - |
| allele genomes 8 | 192 | 170 | 0.818 | 0.780 | 0.883 | 0.645 | 0.500 | 0.850 | 0.151 | 0.876 | 0.676 | 0.530 | 0.064 | 0.470 | 0.886 | 0.357 | 0.500 | 0.850 | 0.555 | 0.500 | 0.906 | - |
| congeners 2 | 288 | 257 | 0.834 | 0.792 | 0.900 | 0.618 | 0.501 | 0.898 | 0.109 | 0.672 | 0.622 | 0.528 | 0.316 | 0.480 | 0.673 | 0.304 | 0.500 | 0.898 | 0.535 | 0.500 | 0.828 | - |
| congeners 6 | 288 | 256 | 0.795 | 0.732 | 0.858 | 0.599 | 0.509 | 0.863 | 0.139 | 0.655 | 0.583 | 0.555 | 0.278 | 0.508 | 0.677 | 0.306 | 0.500 | 0.863 | 0.519 | 0.500 | 0.851 | - |
| twin 0 | 288 | 259 | 0.916 | 0.851 | 0.919 | 0.754 | 0.498 | 0.859 | 0.141 | 0.672 | 0.605 | 0.509 | 0.276 | 0.477 | 0.683 | 0.360 | 0.500 | 0.859 | 0.493 | 0.500 | 0.856 | - |
| twin 1 | 288 | 254 | 0.871 | 0.743 | 0.828 | 0.545 | 0.517 | 0.824 | 0.176 | 0.654 | 0.599 | 0.569 | 0.303 | 0.502 | 0.657 | 0.347 | 0.500 | 0.824 | 0.554 | 0.500 | 0.804 | - |
| regime low | 192 | 168 | 0.787 | 0.718 | 0.865 | 0.593 | 0.506 | 0.814 | 0.186 | 0.675 | 0.613 | 0.596 | 0.310 | 0.520 | 0.666 | 0.393 | 0.500 | 0.814 | 0.505 | 0.500 | 0.778 | - |
| regime high | 192 | 166 | 0.854 | 0.836 | 0.952 | 0.622 | 0.486 | 0.905 | 0.097 | 0.703 | 0.624 | 0.491 | 0.301 | 0.470 | 0.687 | 0.250 | 0.500 | 0.905 | 0.525 | 0.500 | 0.880 | - |
| regime none | 192 | 179 | 0.806 | 0.730 | 0.833 | 0.595 | 0.507 | 0.809 | 0.193 | 0.617 | 0.557 | 0.539 | 0.278 | 0.497 | 0.675 | 0.377 | 0.500 | 0.809 | 0.568 | 0.500 | 0.851 | - |
| depth 1 | 159 | 135 | 0.806 | 0.757 | 0.885 | 0.619 | 0.521 | 0.852 | 0.148 | 0.675 | 0.600 | 0.547 | 0.301 | 0.494 | 0.686 | 0.353 | 0.500 | 0.852 | 0.497 | 0.500 | 0.830 | - |
| depth 2 | 136 | 122 | 0.815 | 0.754 | 0.874 | 0.598 | 0.494 | 0.828 | 0.175 | 0.624 | 0.555 | 0.497 | 0.246 | 0.450 | 0.632 | 0.362 | 0.500 | 0.828 | 0.559 | 0.500 | 0.811 | - |
| depth 4 | 138 | 127 | 0.829 | 0.791 | 0.889 | 0.574 | 0.480 | 0.836 | 0.164 | 0.672 | 0.629 | 0.563 | 0.308 | 0.511 | 0.676 | 0.356 | 0.500 | 0.836 | 0.522 | 0.500 | 0.839 | - |
| depth 0.5 | 143 | 129 | 0.809 | 0.734 | 0.860 | 0.607 | 0.506 | 0.843 | 0.158 | 0.675 | 0.613 | 0.550 | 0.323 | 0.509 | 0.678 | 0.342 | 0.500 | 0.843 | 0.530 | 0.500 | 0.840 | - |
| all, within identity bands of 0.005 | 576 | 513 | 0.621 | 0.580 | 0.723 | 0.660 | 0.688 | 0.718 | 0.284 | 0.635 | 0.637 | 0.583 | 0.331 | 0.494 | 0.577 | 0.379 | 0.500 | 0.718 | 0.530 | 0.500 | 0.698 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 44 | 34 | 0.779 | 0.899 | 1.000 | 0.800 | 0.752 | 0.803 | 0.989 | - |
| Q_deep vs N_0 | congeners 6 | 57 | 41 | 0.819 | 0.940 | 1.000 | 0.706 | 0.683 | 0.730 | 0.998 | - |
| Q_deep vs N_0 | no twin | 48 | 38 | 0.997 | 0.993 | 0.985 | 0.795 | 0.735 | 0.820 | 0.988 | - |
| Q_deep vs N_0 | twin | 53 | 37 | 0.916 | 0.829 | 0.904 | 0.701 | 0.673 | 0.680 | 0.903 | - |
| Q_deep vs N_0 | no allele genomes | 27 | 31 | 0.827 | 0.892 | 0.916 | 0.500 | 0.500 | 0.500 | 0.916 | - |
| Q_deep vs N_0 | allele genomes | 74 | 44 | 0.807 | 0.939 | 0.948 | 0.887 | 0.815 | 0.904 | 0.972 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 11 | 7 | 0.870 | 0.987 | 1.000 | 0.688 | 0.636 | 1.000 | 1.000 | - |
| Q_deep vs N_0 | recombination none | 35 | 25 | 0.835 | 0.915 | 0.930 | 0.706 | 0.680 | 0.727 | 0.953 | - |
| Q_deep vs N_0 | recombination high | 29 | 23 | 0.838 | 0.970 | 1.000 | 0.783 | 0.699 | 0.790 | 0.999 | - |
| Q_deep vs N_0.5 | congeners 2 | 44 | 41 | 0.761 | 0.975 | 0.978 | 0.778 | 0.713 | 0.794 | 0.898 | - |
| Q_deep vs N_0.5 | congeners 6 | 57 | 33 | 0.792 | 0.960 | 0.952 | 0.602 | 0.628 | 0.706 | 0.890 | - |
| Q_deep vs N_0.5 | no twin | 48 | 38 | 0.986 | 0.997 | 0.873 | 0.731 | 0.667 | 0.774 | 0.899 | - |
| Q_deep vs N_0.5 | twin | 53 | 36 | 0.975 | 0.941 | 0.823 | 0.644 | 0.673 | 0.664 | 0.751 | - |
| Q_deep vs N_0.5 | no allele genomes | 27 | 24 | 0.843 | 0.946 | 0.735 | 0.500 | 0.500 | 0.500 | 0.735 | - |
| Q_deep vs N_0.5 | allele genomes | 74 | 50 | 0.759 | 0.974 | 0.888 | 0.810 | 0.784 | 0.878 | 0.900 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 11 | 6 | 0.773 | 0.773 | 0.909 | 0.727 | 0.758 | 0.985 | 1.000 | - |
| Q_deep vs N_0.5 | recombination none | 35 | 22 | 0.816 | 0.945 | 0.817 | 0.716 | 0.718 | 0.761 | 0.908 | - |
| Q_deep vs N_0.5 | recombination high | 29 | 25 | 0.761 | 1.000 | 0.968 | 0.770 | 0.680 | 0.778 | 0.932 | - |
| Q_deep vs N_0.8 | congeners 2 | 44 | 33 | 0.663 | 0.815 | 0.796 | 0.745 | 0.690 | 0.733 | 0.785 | - |
| Q_deep vs N_0.8 | congeners 6 | 57 | 32 | 0.639 | 0.769 | 0.758 | 0.661 | 0.690 | 0.677 | 0.731 | - |
| Q_deep vs N_0.8 | no twin | 48 | 30 | 0.796 | 0.823 | 0.683 | 0.713 | 0.689 | 0.719 | 0.783 | - |
| Q_deep vs N_0.8 | twin | 53 | 35 | 0.805 | 0.767 | 0.740 | 0.646 | 0.684 | 0.643 | 0.681 | - |
| Q_deep vs N_0.8 | no allele genomes | 27 | 23 | 0.713 | 0.760 | 0.681 | 0.500 | 0.500 | 0.500 | 0.681 | - |
| Q_deep vs N_0.8 | allele genomes | 74 | 42 | 0.645 | 0.822 | 0.715 | 0.774 | 0.801 | 0.771 | 0.798 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 11 | 9 | 0.596 | 0.758 | 0.737 | 0.707 | 0.747 | 0.808 | 0.838 | - |
| Q_deep vs N_0.8 | recombination none | 35 | 23 | 0.629 | 0.684 | 0.719 | 0.635 | 0.648 | 0.658 | 0.735 | - |
| Q_deep vs N_0.8 | recombination high | 29 | 17 | 0.625 | 0.982 | 0.880 | 0.765 | 0.722 | 0.759 | 0.824 | - |
| Q_deep vs N_0.95 | congeners 2 | 44 | 38 | 0.566 | 0.678 | 0.623 | 0.706 | 0.671 | 0.689 | 0.699 | - |
| Q_deep vs N_0.95 | congeners 6 | 57 | 35 | 0.595 | 0.607 | 0.591 | 0.618 | 0.665 | 0.601 | 0.593 | - |
| Q_deep vs N_0.95 | no twin | 48 | 37 | 0.651 | 0.680 | 0.617 | 0.681 | 0.683 | 0.634 | 0.655 | - |
| Q_deep vs N_0.95 | twin | 53 | 36 | 0.593 | 0.593 | 0.534 | 0.630 | 0.660 | 0.592 | 0.587 | - |
| Q_deep vs N_0.95 | no allele genomes | 27 | 23 | 0.652 | 0.651 | 0.578 | 0.500 | 0.500 | 0.500 | 0.578 | - |
| Q_deep vs N_0.95 | allele genomes | 74 | 50 | 0.574 | 0.623 | 0.560 | 0.737 | 0.786 | 0.667 | 0.668 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 11 | 7 | 0.623 | 0.688 | 0.636 | 0.727 | 0.727 | 0.831 | 0.831 | - |
| Q_deep vs N_0.95 | recombination none | 35 | 25 | 0.609 | 0.618 | 0.581 | 0.640 | 0.648 | 0.654 | 0.679 | - |
| Q_deep vs N_0.95 | recombination high | 29 | 26 | 0.576 | 0.768 | 0.672 | 0.714 | 0.718 | 0.677 | 0.703 | - |

## L2d: the gain of each group over identity and depth (logistic models, genera held out)

AUC of grouped-CV logistic models on the target rows; gain over identity, top identity and fragments with a 95% interval over resampled genera.
| set | model | auc_base | auc | gain | gain_lo | gain_hi |
|---|---|---|---|---|---|---|
| pe | identity, depth +ancestry | 0.780 | 0.822 | 0.043 | 0.022 | 0.064 |
| pe | identity, depth +alleles | 0.780 | 0.856 | 0.076 | 0.036 | 0.123 |
| pe | identity, depth +polymorphic | 0.780 | 0.798 | 0.019 | -0.014 | 0.047 |
| pe | identity, depth +weights | 0.780 | 0.823 | 0.044 | 0.023 | 0.065 |
| pe | identity, depth +all three | 0.780 | 0.896 | 0.116 | 0.069 | 0.164 |
| pe | identity, depth +ancestry, polymorphic | 0.780 | 0.870 | 0.090 | 0.060 | 0.119 |
| pe | identity, depth +all four | 0.780 | 0.898 | 0.119 | 0.074 | 0.165 |
| pe | identity, depth +oracle fixed agreement | 0.780 | 0.780 | 0.000 | -0.006 | 0.006 |
| pe, genera with alleles | identity, depth +ancestry | 0.784 | 0.828 | 0.044 | 0.015 | 0.071 |
| pe, genera with alleles | identity, depth +alleles | 0.784 | 0.927 | 0.144 | 0.084 | 0.201 |
| pe, genera with alleles | identity, depth +polymorphic | 0.784 | 0.892 | 0.108 | 0.072 | 0.145 |
| pe, genera with alleles | identity, depth +weights | 0.784 | 0.827 | 0.044 | 0.015 | 0.069 |
| pe, genera with alleles | identity, depth +all three | 0.784 | 0.951 | 0.167 | 0.108 | 0.219 |
| pe, genera with alleles | identity, depth +ancestry, polymorphic | 0.784 | 0.915 | 0.132 | 0.092 | 0.171 |
| pe, genera with alleles | identity, depth +all four | 0.784 | 0.950 | 0.166 | 0.106 | 0.219 |
| pe, genera with alleles | identity, depth +oracle fixed agreement | 0.784 | 0.858 | 0.074 | 0.047 | 0.096 |
| se | identity, depth +ancestry | 0.759 | 0.797 | 0.038 | 0.018 | 0.057 |
| se | identity, depth +alleles | 0.759 | 0.847 | 0.087 | 0.046 | 0.135 |
| se | identity, depth +polymorphic | 0.759 | 0.780 | 0.021 | -0.002 | 0.043 |
| se | identity, depth +weights | 0.759 | 0.795 | 0.036 | 0.021 | 0.054 |
| se | identity, depth +all three | 0.759 | 0.885 | 0.125 | 0.084 | 0.171 |
| se | identity, depth +ancestry, polymorphic | 0.759 | 0.848 | 0.088 | 0.061 | 0.115 |
| se | identity, depth +all four | 0.759 | 0.885 | 0.125 | 0.085 | 0.171 |
| se | identity, depth +oracle fixed agreement | 0.759 | 0.759 | 0.000 | 0.000 | 0.000 |
| se, genera with alleles | identity, depth +ancestry | 0.760 | 0.800 | 0.040 | 0.015 | 0.061 |
| se, genera with alleles | identity, depth +alleles | 0.760 | 0.915 | 0.154 | 0.093 | 0.211 |
| se, genera with alleles | identity, depth +polymorphic | 0.760 | 0.864 | 0.104 | 0.073 | 0.138 |
| se, genera with alleles | identity, depth +weights | 0.760 | 0.799 | 0.038 | 0.015 | 0.058 |
| se, genera with alleles | identity, depth +all three | 0.760 | 0.938 | 0.178 | 0.117 | 0.233 |
| se, genera with alleles | identity, depth +ancestry, polymorphic | 0.760 | 0.892 | 0.131 | 0.094 | 0.173 |
| se, genera with alleles | identity, depth +all four | 0.760 | 0.938 | 0.177 | 0.117 | 0.232 |
| se, genera with alleles | identity, depth +oracle fixed agreement | 0.760 | 0.760 | 0.000 | 0.000 | 0.000 |
| hifi | identity, depth +ancestry | 0.805 | 0.853 | 0.047 | 0.016 | 0.074 |
| hifi | identity, depth +alleles | 0.805 | 0.884 | 0.079 | 0.039 | 0.127 |
| hifi | identity, depth +polymorphic | 0.805 | 0.885 | 0.080 | 0.030 | 0.121 |
| hifi | identity, depth +weights | 0.805 | 0.847 | 0.042 | 0.017 | 0.065 |
| hifi | identity, depth +all three | 0.805 | 0.942 | 0.137 | 0.089 | 0.179 |
| hifi | identity, depth +ancestry, polymorphic | 0.805 | 0.916 | 0.110 | 0.070 | 0.149 |
| hifi | identity, depth +all four | 0.805 | 0.943 | 0.138 | 0.092 | 0.181 |
| hifi | identity, depth +oracle fixed agreement | 0.805 | 0.805 | 0.000 | 0.000 | 0.000 |
(8 more rows in l2_increment.tsv)

## Sensitivity: by divergence of the lineage and by recombination with congeners

Per feature, the threshold 5% of all N rows exceed (the feature oriented by its AUC); rate_above: the share of Q rows above it (sensitivity) or of N rows (false-positive rate) in each bin. 'MRCA with nearest allele genome' is how recently the strain shared an ancestor with a strain the database knows; 'genes with a congener's segment' how much of it came from a congener.
| set | class | by | level | rows | identity | ancestry_agreement | allele_explained_share | polymorphic_known_share | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 278 | 0.295 | 0.388 | 0.403 | 0.284 | 0.464 | 0.000 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.078 | 0.141 | 0.078 | 0.078 | 0.078 | 0.000 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 573 | 0.353 | 0.527 | 0.000 | 0.000 | 0.094 | 0.000 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 809 | 0.386 | 0.525 | 0.739 | 0.468 | 0.695 | 0.000 |
| pe | Q | depth | 10 | 435 | 0.352 | 0.501 | 0.405 | 0.246 | 0.414 | 0.000 |
| pe | Q | depth | 100 | 426 | 0.326 | 0.469 | 0.423 | 0.289 | 0.406 | 0.000 |
| pe | Q | depth | 3 | 439 | 0.405 | 0.526 | 0.424 | 0.264 | 0.497 | 0.000 |
| pe | Q | depth | 30 | 424 | 0.309 | 0.460 | 0.408 | 0.276 | 0.422 | 0.000 |
| pe | Q | distance to rep | 0.005-0.015 | 790 | 0.365 | 0.516 | 0.473 | 0.266 | 0.422 | 0.000 |
| pe | Q | distance to rep | 0.015-0.03 | 433 | 0.048 | 0.333 | 0.527 | 0.363 | 0.370 | 0.000 |
| pe | Q | distance to rep | <0.005 | 298 | 0.977 | 0.966 | 0.074 | 0.070 | 0.738 | 0.000 |
| pe | Q | distance to rep | >=0.03 | 203 | 0.005 | 0.020 | 0.448 | 0.369 | 0.182 | 0.000 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 570 | 0.316 | 0.488 | 0.428 | 0.260 | 0.358 | 0.000 |
| pe | Q | genes with a congener's segment | <=0.1 | 454 | 0.333 | 0.476 | 0.407 | 0.225 | 0.467 | 0.000 |
| pe | Q | genes with a congener's segment | >0.3 | 251 | 0.299 | 0.514 | 0.502 | 0.307 | 0.331 | 0.000 |
| pe | Q | genes with a congener's segment | none | 449 | 0.434 | 0.492 | 0.356 | 0.303 | 0.559 | 0.000 |
| pe | Q | regime | high | 575 | 0.350 | 0.520 | 0.445 | 0.228 | 0.350 | 0.000 |
| pe | Q | regime | low | 574 | 0.291 | 0.455 | 0.423 | 0.261 | 0.430 | 0.000 |
| pe | Q | regime | none | 575 | 0.405 | 0.494 | 0.376 | 0.317 | 0.525 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 933 | 0.338 | 0.477 | 0.738 | 0.459 | 0.685 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 589 | 0.368 | 0.518 | 0.000 | 0.002 | 0.092 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 198 | 0.328 | 0.455 | 0.131 | 0.172 | 0.268 | 0.000 |
| pe | Q | role | Q_deep | 310 | 0.126 | 0.168 | 0.432 | 0.323 | 0.313 | 0.000 |
| pe | Q | role | Q_ils | 200 | 0.090 | 0.105 | 0.370 | 0.275 | 0.075 | 0.000 |
| pe | Q | role | Q_imp | 402 | 0.306 | 0.515 | 0.507 | 0.328 | 0.376 | 0.000 |
| pe | Q | role | Q_lone | 277 | 0.397 | 0.632 | 0.101 | 0.079 | 0.220 | 0.000 |
| pe | Q | role | Q_near | 334 | 0.395 | 0.569 | 0.760 | 0.461 | 0.841 | 0.000 |
| pe | Q | role | Q_rep | 201 | 0.891 | 0.990 | 0.104 | 0.000 | 0.721 | 0.000 |
| pe | Q | twin genus | 0 | 864 | 0.243 | 0.543 | 0.407 | 0.222 | 0.443 | 0.000 |
| pe | Q | twin genus | 1 | 860 | 0.455 | 0.436 | 0.422 | 0.315 | 0.427 | 0.000 |
| pe | N | depth | 10 | 382 | 0.055 | 0.026 | 0.050 | 0.065 | 0.047 | 0.000 |
| pe | N | depth | 100 | 377 | 0.011 | 0.029 | 0.011 | 0.005 | 0.019 | 0.000 |
| pe | N | depth | 3 | 354 | 0.121 | 0.119 | 0.102 | 0.110 | 0.110 | 0.000 |
| pe | N | depth | 30 | 381 | 0.018 | 0.031 | 0.042 | 0.024 | 0.029 | 0.000 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 861 | 0.043 | 0.053 | 0.042 | 0.055 | 0.046 | 0.000 |
| pe | N | genes with a segment of S or a congener | >0.3 | 266 | 0.045 | 0.015 | 0.049 | 0.056 | 0.019 | 0.000 |
| pe | N | genes with a segment of S or a congener | none | 367 | 0.071 | 0.068 | 0.071 | 0.035 | 0.082 | 0.000 |
| pe | N | regime | high | 492 | 0.061 | 0.055 | 0.051 | 0.057 | 0.033 | 0.000 |
| pe | N | regime | low | 496 | 0.026 | 0.044 | 0.038 | 0.058 | 0.056 | 0.000 |
| pe | N | regime | none | 506 | 0.063 | 0.051 | 0.061 | 0.036 | 0.061 | 0.000 |
| pe | N | role | N_0 | 196 | 0.010 | 0.000 | 0.031 | 0.061 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 203 | 0.025 | 0.010 | 0.025 | 0.044 | 0.015 | 0.000 |
| pe | N | role | N_0.8 | 215 | 0.042 | 0.028 | 0.047 | 0.047 | 0.028 | 0.000 |
| pe | N | role | N_0.95 | 227 | 0.079 | 0.123 | 0.075 | 0.031 | 0.075 | 0.000 |
| pe | N | role | N_1 | 228 | 0.110 | 0.149 | 0.096 | 0.053 | 0.193 | 0.000 |
| pe | N | role | N_imp0.2 | 213 | 0.033 | 0.009 | 0.033 | 0.070 | 0.009 | 0.000 |
| pe | N | role | N_imp0.4 | 212 | 0.042 | 0.014 | 0.038 | 0.047 | 0.014 | 0.000 |
| pe | N | twin genus | 0 | 759 | 0.003 | 0.063 | 0.030 | 0.028 | 0.042 | 0.000 |
| pe | N | twin genus | 1 | 735 | 0.099 | 0.037 | 0.071 | 0.073 | 0.059 | 0.000 |
| se | Q | MRCA with nearest allele genome | 0.002-0.01 | 275 | 0.280 | 0.364 | 0.400 | 0.247 | 0.076 |  |
| se | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.078 | 0.109 | 0.156 | 0.062 | 0.000 |  |
| se | Q | MRCA with nearest allele genome | no allele genome | 570 | 0.332 | 0.453 | 0.000 | 0.000 | 0.000 |  |
| se | Q | MRCA with nearest allele genome | recent (<0.002) | 805 | 0.352 | 0.477 | 0.716 | 0.415 | 0.130 |  |
| se | Q | depth | 10 | 435 | 0.340 | 0.444 | 0.393 | 0.221 | 0.120 |  |
| se | Q | depth | 100 | 426 | 0.298 | 0.401 | 0.420 | 0.246 | 0.035 |  |
| se | Q | depth | 3 | 429 | 0.364 | 0.485 | 0.385 | 0.238 | 0.051 |  |
| se | Q | depth | 30 | 424 | 0.290 | 0.417 | 0.427 | 0.243 | 0.087 |  |
| se | Q | distance to rep | 0.005-0.015 | 785 | 0.322 | 0.460 | 0.464 | 0.237 | 0.073 |  |
| se | Q | distance to rep | 0.015-0.03 | 429 | 0.051 | 0.249 | 0.503 | 0.322 | 0.051 |  |
| se | Q | distance to rep | <0.005 | 298 | 0.933 | 0.933 | 0.077 | 0.060 | 0.141 |  |
| se | Q | distance to rep | >=0.03 | 202 | 0.005 | 0.015 | 0.460 | 0.317 | 0.025 |  |
| se | Q | genes with a congener's segment | 0.1-0.3 | 567 | 0.323 | 0.429 | 0.422 | 0.217 | 0.065 |  |
| se | Q | genes with a congener's segment | <=0.1 | 453 | 0.272 | 0.406 | 0.395 | 0.188 | 0.071 |  |
| se | Q | genes with a congener's segment | >0.3 | 248 | 0.278 | 0.476 | 0.476 | 0.290 | 0.048 |  |
| se | Q | genes with a congener's segment | none | 446 | 0.401 | 0.457 | 0.359 | 0.283 | 0.101 |  |
| se | Q | regime | high | 571 | 0.336 | 0.473 | 0.431 | 0.191 | 0.060 |  |
| se | Q | regime | low | 573 | 0.258 | 0.389 | 0.408 | 0.220 | 0.072 |  |
| se | Q | regime | none | 570 | 0.375 | 0.449 | 0.379 | 0.300 | 0.089 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 928 | 0.311 | 0.430 | 0.717 | 0.408 | 0.128 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 586 | 0.340 | 0.447 | 0.000 | 0.000 | 0.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | partly | 196 | 0.316 | 0.429 | 0.158 | 0.138 | 0.036 |  |
| se | Q | role | Q_deep | 309 | 0.091 | 0.142 | 0.430 | 0.272 | 0.052 |  |
| se | Q | role | Q_ils | 197 | 0.071 | 0.091 | 0.335 | 0.259 | 0.015 |  |
| se | Q | role | Q_imp | 398 | 0.304 | 0.462 | 0.480 | 0.302 | 0.065 |  |
| se | Q | role | Q_lone | 276 | 0.366 | 0.496 | 0.109 | 0.072 | 0.040 |  |
| se | Q | role | Q_near | 333 | 0.351 | 0.523 | 0.751 | 0.390 | 0.159 |  |
| se | Q | role | Q_rep | 201 | 0.861 | 0.955 | 0.129 | 0.005 | 0.085 |  |
| se | Q | twin genus | 0 | 858 | 0.226 | 0.456 | 0.407 | 0.188 | 0.056 |  |
| se | Q | twin genus | 1 | 856 | 0.421 | 0.418 | 0.405 | 0.286 | 0.091 |  |
| se | N | depth | 10 | 382 | 0.050 | 0.045 | 0.068 | 0.076 | 0.016 |  |
| se | N | depth | 100 | 377 | 0.005 | 0.008 | 0.005 | 0.008 | 0.000 |  |
| se | N | depth | 3 | 351 | 0.131 | 0.145 | 0.100 | 0.105 | 0.023 |  |
| se | N | depth | 30 | 381 | 0.021 | 0.010 | 0.031 | 0.016 | 0.008 |  |
| se | N | genes with a segment of S or a congener | <=0.3 | 863 | 0.051 | 0.053 | 0.044 | 0.052 | 0.009 |  |
| se | N | genes with a segment of S or a congener | >0.3 | 264 | 0.042 | 0.030 | 0.034 | 0.061 | 0.008 |  |
| se | N | genes with a segment of S or a congener | none | 364 | 0.055 | 0.058 | 0.077 | 0.038 | 0.019 |  |
| se | N | regime | high | 489 | 0.063 | 0.053 | 0.053 | 0.063 | 0.008 |  |
| se | N | regime | low | 498 | 0.034 | 0.046 | 0.034 | 0.048 | 0.012 |  |
| se | N | regime | none | 504 | 0.054 | 0.052 | 0.063 | 0.040 | 0.014 |  |
| se | N | role | N_0 | 199 | 0.025 | 0.005 | 0.030 | 0.055 | 0.005 |  |
| se | N | role | N_0.5 | 202 | 0.030 | 0.010 | 0.030 | 0.040 | 0.005 |  |
| se | N | role | N_0.8 | 213 | 0.033 | 0.028 | 0.052 | 0.047 | 0.009 |  |
| se | N | role | N_0.95 | 226 | 0.066 | 0.106 | 0.075 | 0.053 | 0.018 |  |
| se | N | role | N_1 | 227 | 0.106 | 0.128 | 0.097 | 0.040 | 0.031 |  |
| se | N | role | N_imp0.2 | 212 | 0.038 | 0.038 | 0.033 | 0.057 | 0.005 |  |
| se | N | role | N_imp0.4 | 212 | 0.047 | 0.024 | 0.028 | 0.061 | 0.005 |  |
| se | N | twin genus | 0 | 756 | 0.003 | 0.038 | 0.025 | 0.033 | 0.003 |  |
| se | N | twin genus | 1 | 735 | 0.099 | 0.063 | 0.076 | 0.068 | 0.020 |  |
| hifi | Q | MRCA with nearest allele genome | 0.002-0.01 | 97 | 0.361 | 0.381 | 0.515 | 0.381 | 0.454 |  |
| hifi | Q | MRCA with nearest allele genome | long ago (>=0.01) | 20 | 0.200 | 0.300 | 0.000 | 0.150 | 0.200 |  |
| hifi | Q | MRCA with nearest allele genome | no allele genome | 192 | 0.417 | 0.516 | 0.000 | 0.000 | 0.073 |  |
| hifi | Q | MRCA with nearest allele genome | recent (<0.002) | 267 | 0.536 | 0.588 | 0.940 | 0.528 | 0.566 |  |
| hifi | Q | depth | 0.5 | 143 | 0.462 | 0.524 | 0.524 | 0.329 | 0.364 |  |
| hifi | Q | depth | 1 | 159 | 0.453 | 0.484 | 0.535 | 0.302 | 0.384 |  |
| hifi | Q | depth | 2 | 136 | 0.500 | 0.544 | 0.507 | 0.294 | 0.338 |  |
| hifi | Q | depth | 4 | 138 | 0.406 | 0.529 | 0.522 | 0.333 | 0.391 |  |
| hifi | Q | distance to rep | 0.005-0.015 | 274 | 0.595 | 0.518 | 0.573 | 0.314 | 0.350 |  |
| hifi | Q | distance to rep | 0.015-0.03 | 148 | 0.034 | 0.426 | 0.514 | 0.385 | 0.216 |  |
| hifi | Q | distance to rep | <0.005 | 94 | 1.000 | 1.000 | 0.404 | 0.064 | 0.830 |  |
| hifi | Q | distance to rep | >=0.03 | 60 | 0.000 | 0.000 | 0.500 | 0.533 | 0.117 |  |
| hifi | Q | genes with a congener's segment | 0.1-0.3 | 190 | 0.463 | 0.505 | 0.553 | 0.353 | 0.258 |  |
| hifi | Q | genes with a congener's segment | <=0.1 | 152 | 0.441 | 0.454 | 0.513 | 0.289 | 0.467 |  |
| hifi | Q | genes with a congener's segment | >0.3 | 82 | 0.366 | 0.622 | 0.524 | 0.232 | 0.024 |  |
| hifi | Q | genes with a congener's segment | none | 152 | 0.507 | 0.546 | 0.493 | 0.336 | 0.599 |  |
| hifi | Q | regime | high | 192 | 0.495 | 0.583 | 0.573 | 0.286 | 0.260 |  |
| hifi | Q | regime | low | 192 | 0.375 | 0.438 | 0.521 | 0.323 | 0.380 |  |
| hifi | Q | regime | none | 192 | 0.495 | 0.536 | 0.474 | 0.333 | 0.469 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 314 | 0.490 | 0.525 | 0.908 | 0.529 | 0.605 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | no such site | 2 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 198 | 0.434 | 0.505 | 0.000 | 0.000 | 0.071 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | partly | 62 | 0.323 | 0.516 | 0.258 | 0.242 | 0.113 |  |
| hifi | Q | role | Q_deep | 101 | 0.218 | 0.198 | 0.475 | 0.436 | 0.307 |  |
| hifi | Q | role | Q_ils | 74 | 0.149 | 0.122 | 0.473 | 0.432 | 0.027 |  |
| hifi | Q | role | Q_imp | 136 | 0.382 | 0.559 | 0.551 | 0.309 | 0.059 |  |
| hifi | Q | role | Q_lone | 92 | 0.489 | 0.620 | 0.185 | 0.141 | 0.283 |  |
| hifi | Q | role | Q_near | 100 | 0.620 | 0.640 | 1.000 | 0.500 | 0.880 |  |
| hifi | Q | role | Q_rep | 73 | 0.959 | 1.000 | 0.356 | 0.000 | 0.795 |  |
| hifi | Q | twin genus | 0 | 288 | 0.319 | 0.611 | 0.524 | 0.267 | 0.372 |  |
| hifi | Q | twin genus | 1 | 288 | 0.590 | 0.427 | 0.521 | 0.361 | 0.368 |  |
| hifi | N | depth | 0.5 | 129 | 0.047 | 0.054 | 0.054 | 0.070 | 0.047 |  |
| hifi | N | depth | 1 | 135 | 0.044 | 0.030 | 0.052 | 0.067 | 0.037 |  |
| hifi | N | depth | 2 | 122 | 0.074 | 0.066 | 0.057 | 0.033 | 0.049 |  |
| hifi | N | depth | 4 | 127 | 0.039 | 0.055 | 0.039 | 0.031 | 0.071 |  |
| hifi | N | genes with a segment of S or a congener | <=0.3 | 289 | 0.045 | 0.038 | 0.031 | 0.035 | 0.038 |  |
| hifi | N | genes with a segment of S or a congener | >0.3 | 98 | 0.000 | 0.000 | 0.020 | 0.082 | 0.000 |  |
| hifi | N | genes with a segment of S or a congener | none | 126 | 0.103 | 0.119 | 0.119 | 0.063 | 0.119 |  |
| hifi | N | regime | high | 166 | 0.030 | 0.012 | 0.030 | 0.072 | 0.012 |  |
| hifi | N | regime | low | 168 | 0.048 | 0.054 | 0.024 | 0.012 | 0.054 |  |
| hifi | N | regime | none | 179 | 0.073 | 0.084 | 0.095 | 0.067 | 0.084 |  |
| hifi | N | role | N_0 | 75 | 0.027 | 0.000 | 0.013 | 0.013 | 0.000 |  |
| hifi | N | role | N_0.5 | 74 | 0.000 | 0.000 | 0.000 | 0.068 | 0.000 |  |
| hifi | N | role | N_0.8 | 65 | 0.046 | 0.046 | 0.000 | 0.031 | 0.000 |  |
| hifi | N | role | N_0.95 | 73 | 0.110 | 0.137 | 0.096 | 0.068 | 0.055 |  |
| hifi | N | role | N_1 | 72 | 0.181 | 0.181 | 0.222 | 0.042 | 0.306 |  |
| hifi | N | role | N_imp0.2 | 77 | 0.000 | 0.000 | 0.000 | 0.065 | 0.000 |  |
| hifi | N | role | N_imp0.4 | 77 | 0.000 | 0.000 | 0.026 | 0.065 | 0.000 |  |
| hifi | N | twin genus | 0 | 259 | 0.004 | 0.073 | 0.008 | 0.073 | 0.039 |  |
| hifi | N | twin genus | 1 | 254 | 0.098 | 0.028 | 0.094 | 0.028 | 0.063 |  |

## L3: the presence models (trained on the train world, called on the test world)

Models by read type and feature set (base: the default set without ancestry, alleles and polymorphic; _on_noallele / _on_shuffled: the default model on the test world's runs without allele scores and with the shuffled table). Q and N of 'all test rows': present and absent taxa of every row.
| model | stratum | Q | N | auc | tpr_at_5pct_fpr | sensitivity_at_knob | fpr_at_knob | f1_at_knob |
|---|---|---|---|---|---|---|---|---|
| hifi_base | target rows | 576 | 513 | 0.885 | 0.582 | 0.700 | 0.135 |  |
| hifi_base | hard: Q_deep vs N_0.8/0.95 | 101 | 138 | 0.699 | 0.248 | 0.485 | 0.261 |  |
| hifi_base | twin genera | 288 | 254 | 0.850 | 0.514 | 0.628 | 0.118 |  |
| hifi_base | no allele genomes | 192 | 172 | 0.851 | 0.536 | 0.625 | 0.134 |  |
| hifi_base | recombination high | 192 | 166 | 0.947 | 0.729 | 0.891 | 0.211 |  |
| hifi_base | all test rows | 896 | 938 | 0.955 | - | 0.807 | 0.074 | 0.857 |
| hifi_default | target rows | 576 | 513 | 0.941 | 0.802 | 0.847 | 0.080 |  |
| hifi_default | hard: Q_deep vs N_0.8/0.95 | 101 | 138 | 0.846 | 0.535 | 0.723 | 0.159 |  |
| hifi_default | twin genera | 288 | 254 | 0.931 | 0.781 | 0.882 | 0.130 |  |
| hifi_default | no allele genomes | 192 | 172 | 0.858 | 0.516 | 0.714 | 0.157 |  |
| hifi_default | recombination high | 192 | 166 | 0.974 | 0.865 | 0.906 | 0.090 |  |
| hifi_default | all test rows | 896 | 938 | 0.975 | - | 0.902 | 0.044 | 0.926 |
| pe_alleles | target rows | 1724 | 1494 | 0.926 | 0.752 | 0.819 | 0.109 |  |
| pe_alleles | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.807 | 0.542 | 0.652 | 0.163 |  |
| pe_alleles | twin genera | 860 | 735 | 0.934 | 0.745 | 0.840 | 0.129 |  |
| pe_alleles | no allele genomes | 573 | 496 | 0.821 | 0.480 | 0.644 | 0.165 |  |
| pe_alleles | recombination high | 575 | 492 | 0.938 | 0.755 | 0.833 | 0.118 |  |
| pe_alleles | all test rows | 2684 | 3524 | 0.976 | - | 0.883 | 0.049 | 0.907 |
| pe_ancestry+alleles | target rows | 1724 | 1494 | 0.915 | 0.758 | 0.797 | 0.093 |  |
| pe_ancestry+alleles | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.784 | 0.477 | 0.581 | 0.133 |  |
| pe_ancestry+alleles | twin genera | 860 | 735 | 0.924 | 0.757 | 0.786 | 0.075 |  |
| pe_ancestry+alleles | no allele genomes | 573 | 496 | 0.813 | 0.497 | 0.593 | 0.103 |  |
| pe_ancestry+alleles | recombination high | 575 | 492 | 0.934 | 0.757 | 0.830 | 0.106 |  |
| pe_ancestry+alleles | all test rows | 2684 | 3524 | 0.972 | - | 0.869 | 0.042 | 0.903 |
| pe_ancestry | target rows | 1724 | 1494 | 0.858 | 0.516 | 0.719 | 0.169 |  |
| pe_ancestry | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.671 | 0.235 | 0.545 | 0.278 |  |
| pe_ancestry | twin genera | 860 | 735 | 0.835 | 0.464 | 0.644 | 0.155 |  |
| pe_ancestry | no allele genomes | 573 | 496 | 0.823 | 0.480 | 0.637 | 0.139 |  |
| pe_ancestry | recombination high | 575 | 492 | 0.896 | 0.515 | 0.810 | 0.193 |  |
| pe_ancestry | all test rows | 2684 | 3524 | 0.956 | - | 0.820 | 0.074 | 0.855 |
| pe_base | target rows | 1724 | 1494 | 0.858 | 0.502 | 0.715 | 0.169 |  |
| pe_base | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.671 | 0.206 | 0.535 | 0.269 |  |
| pe_base | twin genera | 860 | 735 | 0.829 | 0.459 | 0.636 | 0.165 |  |
| pe_base | no allele genomes | 573 | 496 | 0.826 | 0.497 | 0.634 | 0.135 |  |
| pe_base | recombination high | 575 | 492 | 0.892 | 0.517 | 0.790 | 0.199 |  |
| pe_base | all test rows | 2684 | 3524 | 0.956 | - | 0.816 | 0.075 | 0.853 |
| pe_default | target rows | 1724 | 1494 | 0.921 | 0.740 | 0.807 | 0.095 |  |
| pe_default | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.803 | 0.494 | 0.626 | 0.143 |  |
| pe_default | twin genera | 860 | 735 | 0.930 | 0.751 | 0.807 | 0.091 |  |
| pe_default | no allele genomes | 573 | 496 | 0.818 | 0.482 | 0.614 | 0.135 |  |
| pe_default | recombination high | 575 | 492 | 0.937 | 0.748 | 0.835 | 0.112 |  |
| pe_default | all test rows | 2684 | 3524 | 0.974 | - | 0.876 | 0.043 | 0.907 |
| pe_default_on_noallele | target rows | 1714 | 1467 | 0.917 | 0.730 | 0.785 | 0.093 |  |
| pe_default_on_noallele | hard: Q_deep vs N_0.8/0.95 | 307 | 436 | 0.798 | 0.440 | 0.593 | 0.135 |  |
| pe_default_on_noallele | twin genera | 850 | 710 | 0.922 | 0.725 | 0.767 | 0.093 |  |
| pe_default_on_noallele | no allele genomes | 573 | 496 | 0.816 | 0.478 | 0.600 | 0.125 |  |
| pe_default_on_noallele | recombination high | 572 | 486 | 0.933 | 0.722 | 0.806 | 0.109 |  |
| pe_default_on_noallele | all test rows | 2674 | 3757 | 0.974 | - | 0.862 | 0.039 | 0.900 |
| pe_default_on_shuffled | target rows | 1716 | 1469 | 0.830 | 0.486 | 0.565 | 0.094 |  |
| pe_default_on_shuffled | hard: Q_deep vs N_0.8/0.95 | 307 | 437 | 0.630 | 0.121 | 0.306 | 0.137 |  |
| pe_default_on_shuffled | twin genera | 852 | 712 | 0.823 | 0.464 | 0.536 | 0.091 |  |
| pe_default_on_shuffled | no allele genomes | 573 | 496 | 0.816 | 0.475 | 0.600 | 0.127 |  |
| pe_default_on_shuffled | recombination high | 573 | 487 | 0.872 | 0.527 | 0.625 | 0.111 |  |
| pe_default_on_shuffled | all test rows | 2676 | 3740 | 0.947 | - | 0.721 | 0.040 | 0.812 |
| pe_weights | target rows | 1724 | 1494 | 0.860 | 0.495 | 0.728 | 0.168 |  |
| pe_weights | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.670 | 0.206 | 0.552 | 0.262 |  |
| pe_weights | twin genera | 860 | 735 | 0.841 | 0.478 | 0.656 | 0.150 |  |
| pe_weights | no allele genomes | 573 | 496 | 0.825 | 0.492 | 0.647 | 0.149 |  |
| pe_weights | recombination high | 575 | 492 | 0.897 | 0.506 | 0.810 | 0.189 |  |
| pe_weights | all test rows | 2684 | 3524 | 0.957 | - | 0.825 | 0.073 | 0.859 |
| pe_without weights | target rows | 1724 | 1494 | 0.920 | 0.751 | 0.814 | 0.096 |  |
| pe_without weights | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.797 | 0.477 | 0.642 | 0.140 |  |
| pe_without weights | twin genera | 860 | 735 | 0.933 | 0.748 | 0.826 | 0.098 |  |
| pe_without weights | no allele genomes | 573 | 496 | 0.816 | 0.499 | 0.637 | 0.141 |  |
| pe_without weights | recombination high | 575 | 492 | 0.937 | 0.767 | 0.840 | 0.114 |  |
| pe_without weights | all test rows | 2684 | 3524 | 0.974 | - | 0.881 | 0.044 | 0.909 |

The models' calls on the target rows by bin: 'called' the share called at the model's knob (Q: sensitivity; N: false-positive rate), 'at 5% FPR' the share above the score 5% of the model's N rows exceed.
| set | class | by | level | rows | base called | ancestry called | alleles called | default called | base at 5% FPR | ancestry at 5% FPR | alleles at 5% FPR | default at 5% FPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 278 | 0.709 | 0.727 | 0.813 | 0.802 | 0.511 | 0.525 | 0.755 | 0.734 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.391 | 0.391 | 0.266 | 0.297 | 0.172 | 0.172 | 0.203 | 0.188 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 573 | 0.634 | 0.637 | 0.644 | 0.614 | 0.403 | 0.435 | 0.517 | 0.492 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 809 | 0.800 | 0.801 | 0.989 | 0.985 | 0.595 | 0.597 | 0.962 | 0.962 |
| pe | Q | depth | 10 | 435 | 0.745 | 0.756 | 0.821 | 0.807 | 0.506 | 0.513 | 0.759 | 0.766 |
| pe | Q | depth | 100 | 426 | 0.737 | 0.744 | 0.843 | 0.831 | 0.577 | 0.601 | 0.779 | 0.768 |
| pe | Q | depth | 3 | 439 | 0.656 | 0.645 | 0.779 | 0.765 | 0.401 | 0.410 | 0.688 | 0.663 |
| pe | Q | depth | 30 | 424 | 0.722 | 0.733 | 0.835 | 0.825 | 0.526 | 0.542 | 0.785 | 0.767 |
| pe | Q | distance to rep | 0.005-0.015 | 790 | 0.720 | 0.722 | 0.881 | 0.844 | 0.481 | 0.491 | 0.789 | 0.762 |
| pe | Q | distance to rep | 0.015-0.03 | 433 | 0.649 | 0.667 | 0.764 | 0.788 | 0.372 | 0.395 | 0.684 | 0.688 |
| pe | Q | distance to rep | <0.005 | 298 | 0.983 | 0.987 | 0.990 | 0.990 | 0.940 | 0.963 | 0.980 | 0.973 |
| pe | Q | distance to rep | >=0.03 | 203 | 0.438 | 0.429 | 0.443 | 0.433 | 0.217 | 0.212 | 0.424 | 0.424 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 570 | 0.763 | 0.779 | 0.837 | 0.825 | 0.482 | 0.516 | 0.751 | 0.753 |
| pe | Q | genes with a congener's segment | <=0.1 | 454 | 0.659 | 0.661 | 0.709 | 0.700 | 0.515 | 0.520 | 0.656 | 0.659 |
| pe | Q | genes with a congener's segment | >0.3 | 251 | 0.737 | 0.737 | 0.861 | 0.861 | 0.486 | 0.510 | 0.805 | 0.789 |
| pe | Q | genes with a congener's segment | none | 449 | 0.697 | 0.693 | 0.884 | 0.862 | 0.521 | 0.514 | 0.822 | 0.780 |
| pe | Q | regime | high | 575 | 0.790 | 0.810 | 0.833 | 0.835 | 0.501 | 0.543 | 0.757 | 0.757 |
| pe | Q | regime | low | 574 | 0.657 | 0.653 | 0.726 | 0.707 | 0.493 | 0.497 | 0.666 | 0.664 |
| pe | Q | regime | none | 575 | 0.697 | 0.694 | 0.897 | 0.878 | 0.511 | 0.508 | 0.835 | 0.800 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 933 | 0.773 | 0.782 | 0.977 | 0.976 | 0.565 | 0.569 | 0.950 | 0.952 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 1.000 | 1.000 | 0.750 | 1.000 | 1.000 | 1.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 589 | 0.642 | 0.645 | 0.652 | 0.618 | 0.414 | 0.443 | 0.526 | 0.496 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 198 | 0.652 | 0.636 | 0.566 | 0.566 | 0.460 | 0.470 | 0.490 | 0.465 |
| pe | Q | role | Q_deep | 310 | 0.535 | 0.545 | 0.652 | 0.626 | 0.248 | 0.265 | 0.555 | 0.523 |
| pe | Q | role | Q_ils | 200 | 0.465 | 0.465 | 0.570 | 0.540 | 0.165 | 0.150 | 0.450 | 0.420 |
| pe | Q | role | Q_imp | 402 | 0.741 | 0.734 | 0.883 | 0.871 | 0.500 | 0.517 | 0.821 | 0.811 |
| pe | Q | role | Q_lone | 277 | 0.715 | 0.726 | 0.773 | 0.762 | 0.527 | 0.552 | 0.675 | 0.664 |
| pe | Q | role | Q_near | 334 | 0.835 | 0.850 | 0.991 | 0.994 | 0.668 | 0.668 | 0.973 | 0.973 |
| pe | Q | role | Q_rep | 201 | 0.985 | 0.985 | 0.975 | 0.975 | 0.920 | 0.960 | 0.960 | 0.970 |
| pe | Q | twin genus | 0 | 864 | 0.793 | 0.794 | 0.799 | 0.807 | 0.580 | 0.582 | 0.748 | 0.740 |
| pe | Q | twin genus | 1 | 860 | 0.636 | 0.644 | 0.840 | 0.807 | 0.423 | 0.449 | 0.757 | 0.741 |
| pe | N | depth | 10 | 382 | 0.170 | 0.178 | 0.128 | 0.118 | 0.042 | 0.055 | 0.068 | 0.073 |
| pe | N | depth | 100 | 377 | 0.162 | 0.162 | 0.093 | 0.072 | 0.056 | 0.042 | 0.037 | 0.029 |
| pe | N | depth | 3 | 354 | 0.167 | 0.169 | 0.116 | 0.102 | 0.056 | 0.048 | 0.051 | 0.065 |
| pe | N | depth | 30 | 381 | 0.178 | 0.168 | 0.100 | 0.089 | 0.047 | 0.055 | 0.045 | 0.034 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 861 | 0.195 | 0.189 | 0.093 | 0.091 | 0.057 | 0.064 | 0.037 | 0.042 |
| pe | N | genes with a segment of S or a congener | >0.3 | 266 | 0.060 | 0.064 | 0.071 | 0.053 | 0.011 | 0.011 | 0.026 | 0.038 |
| pe | N | genes with a segment of S or a congener | none | 367 | 0.188 | 0.199 | 0.174 | 0.136 | 0.063 | 0.046 | 0.098 | 0.079 |
| pe | N | regime | high | 492 | 0.199 | 0.193 | 0.118 | 0.112 | 0.047 | 0.059 | 0.055 | 0.057 |
| pe | N | regime | low | 496 | 0.149 | 0.155 | 0.054 | 0.054 | 0.054 | 0.054 | 0.016 | 0.024 |
| pe | N | regime | none | 506 | 0.160 | 0.160 | 0.154 | 0.119 | 0.049 | 0.038 | 0.079 | 0.069 |
| pe | N | role | N_0 | 196 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 203 | 0.020 | 0.010 | 0.005 | 0.005 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.8 | 215 | 0.177 | 0.200 | 0.102 | 0.098 | 0.037 | 0.023 | 0.042 | 0.047 |
| pe | N | role | N_0.95 | 227 | 0.357 | 0.352 | 0.220 | 0.185 | 0.110 | 0.106 | 0.097 | 0.097 |
| pe | N | role | N_1 | 228 | 0.469 | 0.456 | 0.289 | 0.254 | 0.167 | 0.180 | 0.149 | 0.136 |
| pe | N | role | N_imp0.2 | 213 | 0.042 | 0.042 | 0.042 | 0.042 | 0.009 | 0.014 | 0.028 | 0.023 |
| pe | N | role | N_imp0.4 | 212 | 0.066 | 0.071 | 0.071 | 0.052 | 0.009 | 0.009 | 0.019 | 0.033 |
| pe | N | twin genus | 0 | 759 | 0.174 | 0.183 | 0.090 | 0.099 | 0.058 | 0.063 | 0.045 | 0.053 |
| pe | N | twin genus | 1 | 735 | 0.165 | 0.155 | 0.129 | 0.091 | 0.042 | 0.037 | 0.056 | 0.048 |
| hifi | Q | MRCA with nearest allele genome | 0.002-0.01 | 97 | 0.588 |  |  | 0.825 | 0.433 |  |  | 0.753 |
| hifi | Q | MRCA with nearest allele genome | long ago (>=0.01) | 20 | 0.500 |  |  | 0.300 | 0.300 |  |  | 0.200 |
| hifi | Q | MRCA with nearest allele genome | no allele genome | 192 | 0.625 |  |  | 0.714 | 0.536 |  |  | 0.625 |
| hifi | Q | MRCA with nearest allele genome | recent (<0.002) | 267 | 0.809 |  |  | 0.993 | 0.689 |  |  | 0.993 |
| hifi | Q | depth | 0.5 | 143 | 0.636 |  |  | 0.839 | 0.524 |  |  | 0.797 |
| hifi | Q | depth | 1 | 159 | 0.679 |  |  | 0.830 | 0.560 |  |  | 0.780 |
| hifi | Q | depth | 2 | 136 | 0.750 |  |  | 0.868 | 0.625 |  |  | 0.846 |
| hifi | Q | depth | 4 | 138 | 0.739 |  |  | 0.855 | 0.623 |  |  | 0.790 |
| hifi | Q | distance to rep | 0.005-0.015 | 274 | 0.719 |  |  | 0.931 | 0.591 |  |  | 0.872 |
| hifi | Q | distance to rep | 0.015-0.03 | 148 | 0.595 |  |  | 0.736 | 0.439 |  |  | 0.669 |
| hifi | Q | distance to rep | <0.005 | 94 | 1.000 |  |  | 1.000 | 1.000 |  |  | 1.000 |
| hifi | Q | distance to rep | >=0.03 | 60 | 0.400 |  |  | 0.500 | 0.233 |  |  | 0.500 |
| hifi | Q | genes with a congener's segment | 0.1-0.3 | 190 | 0.811 |  |  | 0.900 | 0.695 |  |  | 0.858 |
| hifi | Q | genes with a congener's segment | <=0.1 | 152 | 0.559 |  |  | 0.691 | 0.461 |  |  | 0.632 |
| hifi | Q | genes with a congener's segment | >0.3 | 82 | 0.854 |  |  | 0.902 | 0.659 |  |  | 0.854 |
| hifi | Q | genes with a congener's segment | none | 152 | 0.618 |  |  | 0.908 | 0.520 |  |  | 0.875 |
| hifi | Q | regime | high | 192 | 0.891 |  |  | 0.906 | 0.771 |  |  | 0.865 |
| hifi | Q | regime | low | 192 | 0.568 |  |  | 0.708 | 0.443 |  |  | 0.656 |
| hifi | Q | regime | none | 192 | 0.641 |  |  | 0.927 | 0.531 |  |  | 0.885 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 314 | 0.755 |  |  | 0.984 | 0.611 |  |  | 0.981 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | no such site | 2 | 1.000 |  |  | 1.000 | 1.000 |  |  | 1.000 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 198 | 0.636 |  |  | 0.722 | 0.551 |  |  | 0.636 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | partly | 62 | 0.613 |  |  | 0.548 | 0.516 |  |  | 0.419 |
| hifi | Q | role | Q_deep | 101 | 0.485 |  |  | 0.723 | 0.386 |  |  | 0.604 |
| hifi | Q | role | Q_ils | 74 | 0.432 |  |  | 0.676 | 0.189 |  |  | 0.635 |
| hifi | Q | role | Q_imp | 136 | 0.779 |  |  | 0.904 | 0.610 |  |  | 0.868 |
| hifi | Q | role | Q_lone | 92 | 0.685 |  |  | 0.804 | 0.609 |  |  | 0.739 |
| hifi | Q | role | Q_near | 100 | 0.810 |  |  | 0.980 | 0.740 |  |  | 0.980 |
| hifi | Q | role | Q_rep | 73 | 0.986 |  |  | 0.959 | 0.945 |  |  | 0.959 |
| hifi | Q | twin genus | 0 | 288 | 0.771 |  |  | 0.812 | 0.649 |  |  | 0.771 |
| hifi | Q | twin genus | 1 | 288 | 0.628 |  |  | 0.882 | 0.514 |  |  | 0.833 |
| hifi | N | depth | 0.5 | 129 | 0.124 |  |  | 0.078 | 0.047 |  |  | 0.054 |
| hifi | N | depth | 1 | 135 | 0.126 |  |  | 0.052 | 0.067 |  |  | 0.044 |
| hifi | N | depth | 2 | 122 | 0.156 |  |  | 0.090 | 0.074 |  |  | 0.066 |
| hifi | N | depth | 4 | 127 | 0.134 |  |  | 0.102 | 0.016 |  |  | 0.039 |
| hifi | N | genes with a segment of S or a congener | <=0.3 | 289 | 0.163 |  |  | 0.073 | 0.059 |  |  | 0.038 |
| hifi | N | genes with a segment of S or a congener | >0.3 | 98 | 0.051 |  |  | 0.010 | 0.010 |  |  | 0.000 |
| hifi | N | genes with a segment of S or a congener | none | 126 | 0.135 |  |  | 0.151 | 0.063 |  |  | 0.119 |
| hifi | N | regime | high | 166 | 0.211 |  |  | 0.090 | 0.066 |  |  | 0.048 |
| hifi | N | regime | low | 168 | 0.089 |  |  | 0.036 | 0.042 |  |  | 0.018 |
| hifi | N | regime | none | 179 | 0.106 |  |  | 0.112 | 0.045 |  |  | 0.084 |
| hifi | N | role | N_0 | 75 | 0.000 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_0.5 | 74 | 0.000 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_0.8 | 65 | 0.138 |  |  | 0.154 | 0.000 |  |  | 0.092 |
| hifi | N | role | N_0.95 | 73 | 0.370 |  |  | 0.164 | 0.164 |  |  | 0.110 |
| hifi | N | role | N_1 | 72 | 0.375 |  |  | 0.250 | 0.167 |  |  | 0.167 |
| hifi | N | role | N_imp0.2 | 77 | 0.013 |  |  | 0.000 | 0.013 |  |  | 0.000 |
| hifi | N | role | N_imp0.4 | 77 | 0.065 |  |  | 0.013 | 0.013 |  |  | 0.000 |
| hifi | N | twin genus | 0 | 259 | 0.151 |  |  | 0.031 | 0.050 |  |  | 0.012 |
| hifi | N | twin genus | 1 | 254 | 0.118 |  |  | 0.130 | 0.051 |  |  | 0.091 |
