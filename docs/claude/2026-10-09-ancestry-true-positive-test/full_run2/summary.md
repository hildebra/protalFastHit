# Ancestry true-positive test: /home/falk/atp_full2

Worlds: test (seed 11), train (seed 23); 36 test genera. See docs/claude/2026-10-09-ancestry-true-positive-test/README.md for the design and the predictions.

## Verdicts

| hypothesis | prediction | outcome | verdict |
|---|---|---|---|
| H1a table | stored alleles are the allele genomes' differences (>= 0.95 exact) | min share exact 1.000 | PASS |
| H3 sites | consensus sites = stem + representative's lineage (no recombination, >= 3 congeners): >= 0.98 of them, recall >= 0.95 | share 0.930, recall 0.866 | FAIL |
| H3 fallback | with 2 congeners the nearest congener's own derived states enter (share in history clearly below 1) | mean share 0.502 | PASS |
| H7 novel reads | the allele scores draw some of a novel species' reads onto the target (it carries the alleles' ancestral base at the representative's own sites) | 6118 of 115947 reads (0.053) | INFO |
| H7 settling | the strains' reads the allele scores move go to their species (>= 0.9) | 1.000 of 3239 moved | PASS |
| H1b/H2b oracle | protal's ancestry, allele and polymorphic features equal the oracle's (error-free reads, every row within 1e-6) | worst share within 1e-6: 1.000 | PASS |
| H1 leak | an allele genome's own reads: allele_explained_share >= 0.95, polymorphic_novel_share <= 0.02 | 0.956, 0.000 | PASS |
| H2 calibration | planted A_a: ancestry_agreement = a (+-0.03) | worst |diff| 0.015 | PASS |
| H1 calibration | planted P_f: allele_explained_share rises with f (a step near 0.5: a read is explained once it shares more than half of the allele's edits in its span) | f 0: 0.000, f 0.25: 0.051, f 0.5: 0.286, f 0.75: 0.539, f 1: 0.705 | PASS |
| H4 ancestry, clean strata | the same with >= 3 congeners and no recombination (the consensus applies) | b 0: 0.999, b 0.5: 0.854 | INFO |
| H5 fixed sites | deep strains with their representative-lineage sites covered (>= 0.8) against late novel species (b 0.8, 0.95), genera with alleles: ancestry_fixed_gain adds to ancestry_agreement (CI above 0) | AUC 0.655 -> 0.884, gain 0.230 [0.153, 0.313]; ancestry_fixed_agreement alone 0.837, the true fixed-site agreement alone 0.837 (160 Q, 297 N rows) | PASS |
| H4 ancestry | ancestry_agreement separates deep strains from novel species that left the stem early (b <= 0.5): AUC >= 0.9 | b 0: 0.906, b 0.5: 0.816 | FAIL |
| H6 confound | allele_explained_share favours strains with a stored relative (Q_near) over late twins, but not lone strains | Q_near 0.908, Q_lone 0.320 | PASS |
| H9 decoy | the decoy (a tip of the deep strains' clade labelled a species) against Q_deep: no feature separates (AUC 0.5 +- 0.1) | largest |AUC - 0.5|: top_identity 0.472 (310 Q, 200 decoy rows) | PASS |
| H9b crown lineage | N_1, a novel lineage leaving at S's crown, against Q_deep: the ancestry features cannot tell them apart; the allele features can where Q_deep shares its side of the crown with allele genomes | ancestry_agreement 0.529, allele_explained_share 0.641 | INFO |
| H4/H5 gain | the three groups add separation over identity and depth (pe, CI above 0) | 0.116 [0.067, 0.162] | PASS |
| H5 ceiling | the true fixed-site agreement alone (the oracle's) adds separation; compare with '+all three' for the gap | 0.000 [-0.006, 0.006] | INFO |
| H5 ceiling, genera with alleles | the true fixed-site agreement (the oracle's) over identity and depth, where the fixed sites exist; compare with '+all three' | 0.074 [0.047, 0.096] | INFO |
| H8 calls pe target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.077 [0.045, 0.113] | PASS |
| H8 calls pe hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.156 [0.096, 0.221] | PASS |
| H8 calls hifi target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.070 [0.041, 0.100] | PASS |
| H8 calls hifi hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.172 [0.096, 0.239] | PASS |
| H9 shuffled | shuffled alleles cost the default model separation (the alleles carry information) | AUC 0.923 -> 0.830 | PASS |

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

The build's line: Strain alleles: 7640 alleles (106298 edits) of 2776 gene copies of 24 species, up to 4 each (of 38880 full-reference copies of the database's species and genes: 0 of genomes outside --allele_genome_share 1, 26017 identical to the representative's, 2245 repeated, 45 covering less than 0.5 of it or more than 0.1 apart, 361 as far as the nearest congener's copy or farther; 10212 alleles offered): /home/falk/atp_full2/test/db/strain_alleles.tsv

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
| Q | 0 | 66611 | 0.938 | 0.946 | 589 | 589 | 1.000 |
| Q | 1 | 65713 | 0.887 | 0.927 | 2650 | 2650 | 1.000 |
| N | 0 | 57255 | 0.833 | 0.862 | 1670 | 1670 | 1.000 |
| N | 1 | 58692 | 0.757 | 0.833 | 4448 | 4448 | 1.000 |
| class | twin | from | to | reads |
|---|---|---|---|---|
| D | 0 | congener | target | 11 |
| D | 0 | none | target | 1 |
| D | 0 | sister | target | 104 |
| D | 1 | congener | target | 10 |
| D | 1 | sister | target | 440 |
| N | 0 | congener | target | 149 |
| N | 0 | sister | target | 1521 |
| N | 1 | congener | target | 87 |
| N | 1 | sister | target | 4361 |
| Q | 0 | congener | target | 89 |
| Q | 0 | none | target | 2 |
| Q | 0 | sister | target | 498 |
| Q | 1 | congener | target | 107 |
| Q | 1 | sister | target | 2543 |
The protal log's 'strain alleles:' lines (every sample of the default run): 332593 unsure reads took the shifts, 26074 of them to another species.

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
| default | exact | allele_copy_share | 864 | 0.000 | 1.000 |
| default | exact | allele_sites_per_kb | 864 | 0.000 | 1.000 |
| default | exact | ancestry_fixed_share | 864 | 0.000 | 1.000 |
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
| default | pe | allele_copy_share | 3418 | 0.000 | 1.000 |
| default | pe | allele_sites_per_kb | 3418 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_share | 3418 | 0.000 | 1.000 |
| shuffled | exact | ancestry_sites_per_record | 864 | 0.000 | 1.000 |
| shuffled | exact | ancestry_agreement | 864 | 0.000 | 1.000 |
| shuffled | exact | ancestry_congener_share | 864 | 0.000 | 1.000 |
| shuffled | exact | ancestry_indel_sites_per_record | 864 | 0.000 | 1.000 |
| shuffled | exact | ancestry_indel_congener_share | 864 | 0.000 | 1.000 |
| shuffled | exact | allele_explained_share | 864 | 0.000 | 1.000 |
| shuffled | exact | allele_identity_gain | 864 | 0.000 | 1.000 |
| shuffled | exact | polymorphic_known_share | 864 | 0.000 | 1.000 |
| shuffled | exact | polymorphic_novel_share | 864 | 0.000 | 1.000 |
| shuffled | exact | ancestry_fixed_gain | 864 | 0.000 | 1.000 |
| shuffled | exact | ancestry_fixed_agreement | 864 | 0.000 | 1.000 |
| shuffled | exact | allele_copy_share | 864 | 0.000 | 1.000 |
(16 more rows in l2_oracle.tsv)

## L2b: the planted genomes and the leaky control (error-free reads)

Mean features of the planted genomes (P_f: a share f of the farthest allele genome's differences; A_a: the congeners' base at a share 1 - a of the consensus sites) and of Q_leak (an allele genome itself). truth_agreement: the genome's agreement at the consensus sites by position.
| role | rows | truth_agreement | ancestry_agreement | allele_explained_share | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | identity |
|---|---|---|---|---|---|---|---|---|
| A_a0 | 7 | 0.001 | 0.002 | 0.145 | 0.309 | 0.025 | 0.001 | 0.976 |
| A_a0.25 | 7 | 0.252 | 0.265 | 0.103 | 0.250 | 0.017 | 0.001 | 0.978 |
| A_a0.5 | 6 | 0.500 | 0.502 | 0.057 | 0.169 | 0.015 | -0.000 | 0.980 |
| A_a0.75 | 7 | 0.751 | 0.749 | 0.022 | 0.083 | 0.008 | -0.001 | 0.980 |
| A_a1 | 6 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0 | 7 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0.25 | 6 | 0.865 | 0.926 | 0.051 | 0.187 | 0.000 | 0.071 | 0.981 |
| P_f0.5 | 7 | 0.943 | 0.840 | 0.286 | 0.432 | 0.000 | 0.161 | 0.980 |
| P_f0.75 | 6 | 0.828 | 0.786 | 0.539 | 0.593 | 0.000 | 0.214 | 0.979 |
| P_f1 | 6 | 0.479 | 0.724 | 0.705 | 0.805 | 0.000 | 0.278 | 0.978 |
| Q_leak | 110 | 0.412 | 0.805 | 0.956 | 0.521 | 0.000 | 0.189 | 0.989 |

## L2c: separation, pe (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1724 | 1494 | 0.780 | 0.685 | 0.799 | 0.554 | 0.500 | 0.815 | 0.186 | 0.623 | 0.550 | 0.537 | 0.427 | 0.486 | 0.652 | 0.821 | 0.650 |
| Q_deep vs N_0 | 310 | 196 | 0.797 | 0.782 | 0.874 | 0.524 | 0.554 | 0.906 | 0.091 | 0.708 | 0.650 | 0.643 | 0.462 | 0.683 | 0.754 | 0.928 | 0.742 |
| Q_deep vs N_0.5 | 310 | 203 | 0.715 | 0.674 | 0.807 | 0.537 | 0.507 | 0.816 | 0.184 | 0.688 | 0.636 | 0.640 | 0.447 | 0.593 | 0.709 | 0.835 | 0.709 |
| Q_deep vs N_0.8 | 310 | 215 | 0.613 | 0.614 | 0.675 | 0.530 | 0.510 | 0.685 | 0.321 | 0.659 | 0.626 | 0.627 | 0.486 | 0.552 | 0.656 | 0.719 | 0.656 |
| Q_deep vs N_0.95 | 310 | 227 | 0.529 | 0.544 | 0.540 | 0.527 | 0.518 | 0.576 | 0.423 | 0.656 | 0.637 | 0.646 | 0.488 | 0.541 | 0.595 | 0.623 | 0.595 |
| Q_deep vs N_1 | 310 | 228 | 0.486 | 0.511 | 0.483 | 0.501 | 0.513 | 0.529 | 0.473 | 0.641 | 0.630 | 0.635 | 0.482 | 0.529 | 0.549 | 0.564 | 0.549 |
| Q_near vs N_0.95 | 334 | 227 | 0.745 | 0.645 | 0.727 | 0.560 | 0.494 | 0.762 | 0.239 | 0.908 | 0.783 | 0.752 | 0.464 | 0.598 | 0.887 | 0.925 | 0.887 |
| Q_lone vs N_0.95 | 277 | 227 | 0.762 | 0.664 | 0.773 | 0.565 | 0.490 | 0.774 | 0.227 | 0.320 | 0.295 | 0.299 | 0.363 | 0.282 | 0.345 | 0.589 | 0.345 |
| Q_rep vs N_0 | 201 | 196 | 0.993 | 0.891 | 0.987 | 0.605 | 0.515 | 1.000 | 0.000 | 0.460 | 0.351 | 0.316 | 0.401 | 0.385 | 0.707 | 1.000 | 0.691 |
| Q_imp vs N (b <= 0.95) | 402 | 841 | 0.818 | 0.738 | 0.842 | 0.561 | 0.502 | 0.858 | 0.143 | 0.667 | 0.587 | 0.566 | 0.402 | 0.498 | 0.675 | 0.859 | 0.672 |
| Q (no import) vs N_imp | 1322 | 425 | 0.801 | 0.660 | 0.822 | 0.559 | 0.500 | 0.843 | 0.157 | 0.613 | 0.533 | 0.525 | 0.439 | 0.492 | 0.673 | 0.868 | 0.673 |
| Q_ils vs N_0.95 | 200 | 227 | 0.528 | 0.537 | 0.510 | 0.517 | 0.501 | 0.488 | 0.512 | 0.586 | 0.571 | 0.592 | 0.478 | 0.486 | 0.465 | 0.479 | 0.465 |
| decoy: Q_deep vs decoy | 310 | 200 | 0.478 | 0.472 | 0.485 | 0.488 | 0.520 | 0.495 | 0.503 | 0.524 | 0.524 | 0.513 | 0.520 | 0.512 | 0.526 | 0.520 | 0.527 |
| allele genomes 0 | 573 | 496 | 0.767 | 0.684 | 0.798 | 0.534 | 0.510 | 0.807 | 0.193 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.807 | 0.500 |
| allele genomes 2 | 576 | 497 | 0.803 | 0.680 | 0.808 | 0.570 | 0.504 | 0.826 | 0.177 | 0.788 | 0.648 | 0.616 | 0.397 | 0.469 | 0.839 | 0.866 | 0.839 |
| allele genomes 8 | 575 | 501 | 0.772 | 0.692 | 0.794 | 0.559 | 0.488 | 0.814 | 0.185 | 0.770 | 0.585 | 0.575 | 0.338 | 0.482 | 0.834 | 0.893 | 0.834 |
| congeners 2 | 862 | 746 | 0.796 | 0.702 | 0.795 | 0.553 | 0.480 | 0.838 | 0.164 | 0.628 | 0.563 | 0.529 | 0.422 | 0.470 | 0.646 | 0.819 | 0.646 |
| congeners 6 | 862 | 748 | 0.770 | 0.669 | 0.802 | 0.557 | 0.525 | 0.850 | 0.153 | 0.619 | 0.541 | 0.545 | 0.431 | 0.502 | 0.665 | 0.844 | 0.661 |
| twin 0 | 864 | 759 | 0.874 | 0.772 | 0.869 | 0.667 | 0.496 | 0.827 | 0.173 | 0.634 | 0.545 | 0.510 | 0.377 | 0.483 | 0.679 | 0.840 | 0.677 |
| twin 1 | 860 | 735 | 0.789 | 0.607 | 0.713 | 0.510 | 0.517 | 0.803 | 0.197 | 0.615 | 0.558 | 0.562 | 0.469 | 0.487 | 0.632 | 0.806 | 0.630 |
| regime low | 574 | 496 | 0.772 | 0.694 | 0.798 | 0.540 | 0.499 | 0.802 | 0.199 | 0.638 | 0.559 | 0.563 | 0.410 | 0.496 | 0.655 | 0.784 | 0.654 |
| regime high | 575 | 492 | 0.796 | 0.694 | 0.833 | 0.574 | 0.501 | 0.854 | 0.144 | 0.647 | 0.541 | 0.493 | 0.439 | 0.460 | 0.653 | 0.839 | 0.652 |
| regime none | 575 | 506 | 0.774 | 0.670 | 0.769 | 0.551 | 0.499 | 0.797 | 0.205 | 0.588 | 0.551 | 0.555 | 0.431 | 0.502 | 0.651 | 0.850 | 0.648 |
| depth 3 | 439 | 354 | 0.758 | 0.754 | 0.753 | 0.542 | 0.475 | 0.783 | 0.219 | 0.616 | 0.550 | 0.528 | 0.465 | 0.481 | 0.630 | 0.781 | 0.627 |
| depth 10 | 435 | 382 | 0.779 | 0.717 | 0.805 | 0.566 | 0.521 | 0.837 | 0.163 | 0.613 | 0.552 | 0.540 | 0.446 | 0.491 | 0.655 | 0.833 | 0.652 |
| depth 30 | 424 | 381 | 0.794 | 0.680 | 0.829 | 0.552 | 0.500 | 0.819 | 0.181 | 0.626 | 0.546 | 0.534 | 0.390 | 0.484 | 0.661 | 0.834 | 0.660 |
| depth 100 | 426 | 377 | 0.791 | 0.630 | 0.816 | 0.554 | 0.503 | 0.821 | 0.178 | 0.641 | 0.555 | 0.548 | 0.358 | 0.492 | 0.677 | 0.844 | 0.677 |
| all, within identity bands of 0.005 | 1724 | 1494 | 0.561 | 0.545 | 0.643 | 0.642 | 0.688 | 0.731 | 0.271 | 0.677 | 0.674 | 0.626 | 0.502 | 0.552 | 0.620 | 0.762 | 0.620 |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 158 | 97 | 0.796 | 0.860 | 0.969 | 0.743 | 0.648 | 0.776 | 0.949 | 0.776 |
| Q_deep vs N_0 | congeners 6 | 152 | 99 | 0.801 | 0.890 | 0.990 | 0.677 | 0.636 | 0.766 | 0.995 | 0.741 |
| Q_deep vs N_0 | no twin | 162 | 100 | 0.988 | 0.979 | 0.907 | 0.769 | 0.667 | 0.806 | 0.953 | 0.796 |
| Q_deep vs N_0 | twin | 148 | 96 | 0.926 | 0.790 | 0.909 | 0.658 | 0.611 | 0.711 | 0.906 | 0.698 |
| Q_deep vs N_0 | no allele genomes | 82 | 65 | 0.821 | 0.855 | 0.866 | 0.500 | 0.500 | 0.500 | 0.866 | 0.500 |
| Q_deep vs N_0 | allele genomes | 228 | 131 | 0.787 | 0.886 | 0.922 | 0.839 | 0.713 | 0.922 | 0.965 | 0.922 |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 22 | 0.735 | 0.842 | 0.999 | 0.679 | 0.711 | 1.000 | 1.000 | 1.000 |
| Q_deep vs N_0 | recombination none | 106 | 69 | 0.776 | 0.866 | 0.868 | 0.670 | 0.659 | 0.770 | 0.940 | 0.754 |
| Q_deep vs N_0 | recombination high | 102 | 60 | 0.792 | 0.875 | 0.965 | 0.762 | 0.626 | 0.776 | 0.980 | 0.766 |
| Q_deep vs N_0.5 | congeners 2 | 158 | 102 | 0.705 | 0.802 | 0.864 | 0.728 | 0.656 | 0.735 | 0.839 | 0.735 |
| Q_deep vs N_0.5 | congeners 6 | 152 | 101 | 0.733 | 0.812 | 0.907 | 0.648 | 0.623 | 0.709 | 0.895 | 0.709 |
| Q_deep vs N_0.5 | no twin | 162 | 102 | 0.922 | 0.905 | 0.816 | 0.734 | 0.663 | 0.776 | 0.887 | 0.776 |
| Q_deep vs N_0.5 | twin | 148 | 101 | 0.816 | 0.721 | 0.811 | 0.653 | 0.612 | 0.656 | 0.781 | 0.655 |
| Q_deep vs N_0.5 | no allele genomes | 82 | 69 | 0.741 | 0.804 | 0.768 | 0.500 | 0.500 | 0.500 | 0.768 | 0.500 |
| Q_deep vs N_0.5 | allele genomes | 228 | 134 | 0.708 | 0.814 | 0.843 | 0.807 | 0.711 | 0.853 | 0.890 | 0.853 |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.580 | 0.709 | 0.854 | 0.640 | 0.731 | 0.868 | 0.899 | 0.871 |
| Q_deep vs N_0.5 | recombination none | 106 | 68 | 0.690 | 0.764 | 0.799 | 0.656 | 0.670 | 0.710 | 0.855 | 0.711 |
| Q_deep vs N_0.5 | recombination high | 102 | 73 | 0.702 | 0.825 | 0.854 | 0.743 | 0.619 | 0.732 | 0.893 | 0.732 |
| Q_deep vs N_0.8 | congeners 2 | 158 | 109 | 0.607 | 0.645 | 0.676 | 0.706 | 0.648 | 0.679 | 0.721 | 0.679 |
| Q_deep vs N_0.8 | congeners 6 | 152 | 106 | 0.628 | 0.703 | 0.753 | 0.610 | 0.602 | 0.643 | 0.738 | 0.644 |
| Q_deep vs N_0.8 | no twin | 162 | 112 | 0.711 | 0.718 | 0.665 | 0.695 | 0.652 | 0.703 | 0.749 | 0.703 |
| Q_deep vs N_0.8 | twin | 148 | 103 | 0.670 | 0.632 | 0.714 | 0.627 | 0.599 | 0.624 | 0.685 | 0.624 |
| Q_deep vs N_0.8 | no allele genomes | 82 | 67 | 0.637 | 0.664 | 0.679 | 0.500 | 0.500 | 0.500 | 0.679 | 0.500 |
| Q_deep vs N_0.8 | allele genomes | 228 | 148 | 0.606 | 0.686 | 0.687 | 0.769 | 0.707 | 0.762 | 0.789 | 0.762 |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 24 | 0.582 | 0.625 | 0.669 | 0.612 | 0.684 | 0.741 | 0.759 | 0.747 |
| Q_deep vs N_0.8 | recombination none | 106 | 74 | 0.603 | 0.590 | 0.657 | 0.637 | 0.654 | 0.657 | 0.733 | 0.657 |
| Q_deep vs N_0.8 | recombination high | 102 | 74 | 0.614 | 0.759 | 0.743 | 0.694 | 0.604 | 0.682 | 0.754 | 0.682 |
| Q_deep vs N_0.95 | congeners 2 | 158 | 112 | 0.523 | 0.536 | 0.557 | 0.695 | 0.656 | 0.604 | 0.643 | 0.605 |
| Q_deep vs N_0.95 | congeners 6 | 152 | 115 | 0.537 | 0.546 | 0.588 | 0.617 | 0.636 | 0.576 | 0.608 | 0.576 |
| Q_deep vs N_0.95 | no twin | 162 | 116 | 0.551 | 0.577 | 0.571 | 0.675 | 0.661 | 0.633 | 0.635 | 0.633 |
| Q_deep vs N_0.95 | twin | 148 | 111 | 0.554 | 0.500 | 0.584 | 0.643 | 0.632 | 0.571 | 0.608 | 0.571 |
| Q_deep vs N_0.95 | no allele genomes | 82 | 78 | 0.523 | 0.511 | 0.526 | 0.500 | 0.500 | 0.500 | 0.526 | 0.500 |
| Q_deep vs N_0.95 | allele genomes | 228 | 149 | 0.544 | 0.559 | 0.603 | 0.739 | 0.722 | 0.615 | 0.665 | 0.615 |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.499 | 0.616 | 0.552 | 0.601 | 0.719 | 0.572 | 0.577 | 0.569 |
| Q_deep vs N_0.95 | recombination none | 106 | 76 | 0.508 | 0.526 | 0.559 | 0.648 | 0.673 | 0.607 | 0.632 | 0.606 |
| Q_deep vs N_0.95 | recombination high | 102 | 73 | 0.533 | 0.544 | 0.598 | 0.685 | 0.642 | 0.619 | 0.660 | 0.619 |

## L2c: separation, se (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1714 | 1491 | 0.759 | 0.675 | 0.769 | 0.564 | 0.505 | 0.784 | 0.217 | 0.621 | 0.556 | 0.539 | 0.441 | 0.487 | 0.644 | 0.803 | - |
| Q_deep vs N_0 | 309 | 199 | 0.762 | 0.752 | 0.843 | 0.554 | 0.567 | 0.875 | 0.139 | 0.691 | 0.641 | 0.631 | 0.492 | 0.651 | 0.746 | 0.909 | - |
| Q_deep vs N_0.5 | 309 | 202 | 0.667 | 0.653 | 0.728 | 0.553 | 0.523 | 0.765 | 0.227 | 0.680 | 0.641 | 0.631 | 0.481 | 0.602 | 0.711 | 0.838 | - |
| Q_deep vs N_0.8 | 309 | 213 | 0.604 | 0.584 | 0.652 | 0.548 | 0.510 | 0.656 | 0.342 | 0.652 | 0.619 | 0.626 | 0.493 | 0.546 | 0.632 | 0.711 | - |
| Q_deep vs N_0.95 | 309 | 226 | 0.526 | 0.532 | 0.532 | 0.538 | 0.512 | 0.562 | 0.434 | 0.639 | 0.621 | 0.622 | 0.492 | 0.531 | 0.585 | 0.609 | - |
| Q_deep vs N_1 | 309 | 227 | 0.465 | 0.486 | 0.467 | 0.518 | 0.521 | 0.504 | 0.489 | 0.642 | 0.633 | 0.638 | 0.485 | 0.549 | 0.560 | 0.565 | - |
| Q_near vs N_0.95 | 333 | 226 | 0.735 | 0.644 | 0.702 | 0.578 | 0.484 | 0.734 | 0.267 | 0.889 | 0.770 | 0.731 | 0.456 | 0.579 | 0.865 | 0.886 | - |
| Q_lone vs N_0.95 | 276 | 226 | 0.759 | 0.663 | 0.748 | 0.567 | 0.486 | 0.737 | 0.260 | 0.321 | 0.295 | 0.301 | 0.390 | 0.284 | 0.328 | 0.555 | - |
| Q_rep vs N_0 | 201 | 199 | 0.984 | 0.883 | 0.983 | 0.615 | 0.539 | 0.999 | 0.016 | 0.466 | 0.369 | 0.322 | 0.448 | 0.374 | 0.700 | 0.991 | - |
| Q_imp vs N (b <= 0.95) | 398 | 840 | 0.794 | 0.721 | 0.815 | 0.568 | 0.508 | 0.837 | 0.167 | 0.659 | 0.588 | 0.565 | 0.410 | 0.493 | 0.666 | 0.850 | - |
| Q (no import) vs N_imp | 1316 | 424 | 0.779 | 0.658 | 0.784 | 0.567 | 0.503 | 0.798 | 0.201 | 0.615 | 0.544 | 0.528 | 0.449 | 0.494 | 0.662 | 0.833 | - |
| Q_ils vs N_0.95 | 197 | 226 | 0.529 | 0.550 | 0.514 | 0.507 | 0.489 | 0.490 | 0.508 | 0.559 | 0.553 | 0.565 | 0.492 | 0.473 | 0.463 | 0.463 | - |
| decoy: Q_deep vs decoy | 309 | 198 | 0.460 | 0.462 | 0.455 | 0.503 | 0.512 | 0.493 | 0.503 | 0.521 | 0.527 | 0.515 | 0.513 | 0.522 | 0.522 | 0.536 | - |
| allele genomes 0 | 570 | 497 | 0.751 | 0.680 | 0.775 | 0.557 | 0.518 | 0.774 | 0.226 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.774 | - |
| allele genomes 2 | 572 | 496 | 0.773 | 0.676 | 0.772 | 0.574 | 0.512 | 0.800 | 0.201 | 0.791 | 0.667 | 0.628 | 0.410 | 0.481 | 0.810 | 0.854 | - |
| allele genomes 8 | 572 | 498 | 0.757 | 0.670 | 0.765 | 0.563 | 0.488 | 0.780 | 0.220 | 0.747 | 0.583 | 0.565 | 0.392 | 0.474 | 0.827 | 0.865 | - |
| congeners 2 | 857 | 747 | 0.771 | 0.685 | 0.767 | 0.570 | 0.480 | 0.802 | 0.202 | 0.627 | 0.569 | 0.533 | 0.441 | 0.476 | 0.638 | 0.804 | - |
| congeners 6 | 857 | 744 | 0.750 | 0.667 | 0.772 | 0.557 | 0.538 | 0.818 | 0.186 | 0.616 | 0.547 | 0.544 | 0.441 | 0.499 | 0.657 | 0.824 | - |
| twin 0 | 858 | 756 | 0.854 | 0.753 | 0.841 | 0.647 | 0.509 | 0.809 | 0.191 | 0.634 | 0.552 | 0.517 | 0.392 | 0.486 | 0.675 | 0.832 | - |
| twin 1 | 856 | 735 | 0.743 | 0.602 | 0.682 | 0.550 | 0.514 | 0.760 | 0.242 | 0.611 | 0.563 | 0.559 | 0.484 | 0.488 | 0.621 | 0.777 | - |
| regime low | 573 | 498 | 0.743 | 0.673 | 0.760 | 0.547 | 0.508 | 0.772 | 0.231 | 0.638 | 0.569 | 0.567 | 0.422 | 0.496 | 0.641 | 0.767 | - |
| regime high | 571 | 489 | 0.778 | 0.676 | 0.802 | 0.586 | 0.503 | 0.822 | 0.176 | 0.636 | 0.543 | 0.497 | 0.449 | 0.462 | 0.650 | 0.824 | - |
| regime none | 570 | 504 | 0.758 | 0.678 | 0.750 | 0.560 | 0.506 | 0.768 | 0.234 | 0.593 | 0.555 | 0.552 | 0.451 | 0.501 | 0.644 | 0.826 | - |
| depth 3 | 429 | 351 | 0.719 | 0.725 | 0.712 | 0.544 | 0.487 | 0.723 | 0.281 | 0.602 | 0.567 | 0.536 | 0.484 | 0.494 | 0.620 | 0.741 | - |
| depth 10 | 435 | 382 | 0.764 | 0.727 | 0.774 | 0.569 | 0.523 | 0.817 | 0.184 | 0.611 | 0.550 | 0.534 | 0.458 | 0.486 | 0.650 | 0.826 | - |
| depth 30 | 424 | 381 | 0.778 | 0.686 | 0.802 | 0.568 | 0.505 | 0.800 | 0.200 | 0.639 | 0.554 | 0.541 | 0.405 | 0.482 | 0.649 | 0.822 | - |
| depth 100 | 426 | 377 | 0.779 | 0.626 | 0.801 | 0.564 | 0.508 | 0.804 | 0.195 | 0.646 | 0.558 | 0.549 | 0.386 | 0.490 | 0.677 | 0.843 | - |
| all, within identity bands of 0.005 | 1714 | 1491 | 0.553 | 0.550 | 0.632 | 0.652 | 0.660 | 0.700 | 0.304 | 0.678 | 0.676 | 0.633 | 0.507 | 0.561 | 0.618 | 0.742 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 158 | 99 | 0.750 | 0.824 | 0.916 | 0.729 | 0.649 | 0.759 | 0.921 | - |
| Q_deep vs N_0 | congeners 6 | 151 | 100 | 0.775 | 0.862 | 0.974 | 0.652 | 0.612 | 0.760 | 0.979 | - |
| Q_deep vs N_0 | no twin | 161 | 99 | 0.974 | 0.938 | 0.889 | 0.750 | 0.668 | 0.803 | 0.948 | - |
| Q_deep vs N_0 | twin | 148 | 100 | 0.847 | 0.781 | 0.864 | 0.639 | 0.594 | 0.701 | 0.876 | - |
| Q_deep vs N_0 | no allele genomes | 82 | 68 | 0.783 | 0.848 | 0.847 | 0.500 | 0.500 | 0.500 | 0.847 | - |
| Q_deep vs N_0 | allele genomes | 227 | 131 | 0.758 | 0.857 | 0.894 | 0.809 | 0.692 | 0.909 | 0.943 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 22 | 0.730 | 0.865 | 0.999 | 0.669 | 0.626 | 1.000 | 1.000 | - |
| Q_deep vs N_0 | recombination none | 105 | 70 | 0.765 | 0.839 | 0.847 | 0.671 | 0.648 | 0.776 | 0.939 | - |
| Q_deep vs N_0 | recombination high | 101 | 61 | 0.735 | 0.826 | 0.912 | 0.722 | 0.626 | 0.764 | 0.950 | - |
| Q_deep vs N_0.5 | congeners 2 | 158 | 100 | 0.648 | 0.706 | 0.757 | 0.726 | 0.652 | 0.726 | 0.839 | - |
| Q_deep vs N_0.5 | congeners 6 | 151 | 102 | 0.692 | 0.748 | 0.872 | 0.635 | 0.609 | 0.718 | 0.892 | - |
| Q_deep vs N_0.5 | no twin | 161 | 101 | 0.861 | 0.831 | 0.779 | 0.731 | 0.674 | 0.777 | 0.884 | - |
| Q_deep vs N_0.5 | twin | 148 | 101 | 0.704 | 0.646 | 0.749 | 0.636 | 0.587 | 0.660 | 0.792 | - |
| Q_deep vs N_0.5 | no allele genomes | 82 | 68 | 0.669 | 0.758 | 0.714 | 0.500 | 0.500 | 0.500 | 0.714 | - |
| Q_deep vs N_0.5 | allele genomes | 227 | 134 | 0.675 | 0.726 | 0.799 | 0.794 | 0.697 | 0.858 | 0.897 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.605 | 0.624 | 0.834 | 0.650 | 0.697 | 0.909 | 0.947 | - |
| Q_deep vs N_0.5 | recombination none | 105 | 66 | 0.656 | 0.691 | 0.759 | 0.652 | 0.650 | 0.729 | 0.860 | - |
| Q_deep vs N_0.5 | recombination high | 101 | 72 | 0.659 | 0.761 | 0.775 | 0.721 | 0.619 | 0.734 | 0.890 | - |
| Q_deep vs N_0.8 | congeners 2 | 158 | 111 | 0.589 | 0.633 | 0.609 | 0.707 | 0.662 | 0.648 | 0.717 | - |
| Q_deep vs N_0.8 | congeners 6 | 151 | 102 | 0.622 | 0.673 | 0.743 | 0.594 | 0.586 | 0.627 | 0.732 | - |
| Q_deep vs N_0.8 | no twin | 161 | 112 | 0.696 | 0.694 | 0.656 | 0.694 | 0.648 | 0.692 | 0.749 | - |
| Q_deep vs N_0.8 | twin | 148 | 101 | 0.615 | 0.609 | 0.658 | 0.618 | 0.601 | 0.587 | 0.672 | - |
| Q_deep vs N_0.8 | no allele genomes | 82 | 66 | 0.628 | 0.694 | 0.657 | 0.500 | 0.500 | 0.500 | 0.657 | - |
| Q_deep vs N_0.8 | allele genomes | 227 | 147 | 0.599 | 0.640 | 0.658 | 0.758 | 0.705 | 0.718 | 0.768 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.606 | 0.585 | 0.681 | 0.601 | 0.632 | 0.675 | 0.731 | - |
| Q_deep vs N_0.8 | recombination none | 105 | 72 | 0.585 | 0.587 | 0.638 | 0.641 | 0.648 | 0.624 | 0.731 | - |
| Q_deep vs N_0.8 | recombination high | 101 | 74 | 0.616 | 0.717 | 0.710 | 0.676 | 0.604 | 0.648 | 0.746 | - |
| Q_deep vs N_0.95 | congeners 2 | 158 | 111 | 0.530 | 0.515 | 0.531 | 0.681 | 0.639 | 0.597 | 0.612 | - |
| Q_deep vs N_0.95 | congeners 6 | 151 | 115 | 0.525 | 0.551 | 0.578 | 0.597 | 0.601 | 0.568 | 0.606 | - |
| Q_deep vs N_0.95 | no twin | 161 | 116 | 0.549 | 0.567 | 0.568 | 0.679 | 0.654 | 0.632 | 0.651 | - |
| Q_deep vs N_0.95 | twin | 148 | 110 | 0.534 | 0.501 | 0.556 | 0.605 | 0.591 | 0.552 | 0.567 | - |
| Q_deep vs N_0.95 | no allele genomes | 82 | 76 | 0.548 | 0.568 | 0.529 | 0.500 | 0.500 | 0.500 | 0.529 | - |
| Q_deep vs N_0.95 | allele genomes | 227 | 150 | 0.525 | 0.525 | 0.580 | 0.725 | 0.685 | 0.602 | 0.640 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.470 | 0.561 | 0.557 | 0.602 | 0.663 | 0.559 | 0.555 | - |
| Q_deep vs N_0.95 | recombination none | 105 | 77 | 0.503 | 0.516 | 0.546 | 0.630 | 0.650 | 0.603 | 0.623 | - |
| Q_deep vs N_0.95 | recombination high | 101 | 71 | 0.521 | 0.524 | 0.577 | 0.673 | 0.616 | 0.592 | 0.628 | - |

## L2c: separation, hifi (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 576 | 513 | 0.814 | 0.758 | 0.877 | 0.599 | 0.500 | 0.840 | 0.162 | 0.662 | 0.600 | 0.540 | 0.297 | 0.491 | 0.669 | 0.829 | - |
| Q_deep vs N_0 | 101 | 75 | 0.805 | 0.789 | 0.924 | 0.665 | 0.489 | 0.939 | 0.057 | 0.745 | 0.712 | 0.708 | 0.458 | 0.709 | 0.747 | 0.946 | - |
| Q_deep vs N_0.5 | 101 | 74 | 0.777 | 0.756 | 0.967 | 0.602 | 0.453 | 0.849 | 0.153 | 0.682 | 0.670 | 0.669 | 0.362 | 0.608 | 0.716 | 0.827 | - |
| Q_deep vs N_0.8 | 101 | 65 | 0.657 | 0.716 | 0.791 | 0.600 | 0.488 | 0.708 | 0.298 | 0.673 | 0.671 | 0.686 | 0.387 | 0.572 | 0.672 | 0.724 | - |
| Q_deep vs N_0.95 | 101 | 73 | 0.584 | 0.637 | 0.629 | 0.555 | 0.456 | 0.568 | 0.431 | 0.643 | 0.652 | 0.668 | 0.416 | 0.525 | 0.608 | 0.613 | - |
| Q_deep vs N_1 | 101 | 72 | 0.538 | 0.617 | 0.579 | 0.523 | 0.476 | 0.531 | 0.470 | 0.625 | 0.643 | 0.660 | 0.367 | 0.469 | 0.549 | 0.524 | - |
| Q_near vs N_0.95 | 100 | 73 | 0.812 | 0.833 | 0.805 | 0.619 | 0.484 | 0.785 | 0.217 | 1.000 | 0.883 | 0.717 | 0.285 | 0.592 | 0.976 | 0.981 | - |
| Q_lone vs N_0.95 | 92 | 73 | 0.748 | 0.770 | 0.769 | 0.531 | 0.456 | 0.758 | 0.244 | 0.338 | 0.305 | 0.305 | 0.217 | 0.274 | 0.374 | 0.601 | - |
| Q_rep vs N_0 | 73 | 75 | 0.991 | 0.960 | 0.999 | 0.718 | 0.514 | 1.000 | 0.000 | 0.652 | 0.500 | 0.380 | 0.380 | 0.477 | 0.783 | 1.000 | - |
| Q_imp vs N (b <= 0.95) | 136 | 287 | 0.825 | 0.822 | 0.908 | 0.612 | 0.548 | 0.877 | 0.124 | 0.660 | 0.618 | 0.534 | 0.288 | 0.478 | 0.625 | 0.833 | - |
| Q (no import) vs N_imp | 440 | 154 | 0.857 | 0.657 | 0.931 | 0.607 | 0.482 | 0.877 | 0.125 | 0.665 | 0.587 | 0.539 | 0.286 | 0.511 | 0.714 | 0.893 | - |
| Q_ils vs N_0.95 | 74 | 73 | 0.547 | 0.587 | 0.535 | 0.470 | 0.483 | 0.502 | 0.501 | 0.630 | 0.631 | 0.635 | 0.393 | 0.511 | 0.485 | 0.489 | - |
| decoy: Q_deep vs decoy | 101 | 63 | 0.476 | 0.468 | 0.487 | 0.527 | 0.490 | 0.447 | 0.551 | 0.529 | 0.534 | 0.534 | 0.534 | 0.514 | 0.521 | 0.485 | - |
| allele genomes 0 | 192 | 172 | 0.806 | 0.740 | 0.862 | 0.585 | 0.495 | 0.829 | 0.171 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.829 | - |
| allele genomes 2 | 192 | 171 | 0.824 | 0.753 | 0.886 | 0.585 | 0.505 | 0.848 | 0.153 | 0.869 | 0.775 | 0.660 | 0.139 | 0.496 | 0.869 | 0.873 | - |
| allele genomes 8 | 192 | 170 | 0.818 | 0.780 | 0.883 | 0.634 | 0.500 | 0.850 | 0.151 | 0.876 | 0.676 | 0.530 | 0.064 | 0.470 | 0.886 | 0.906 | - |
| congeners 2 | 288 | 257 | 0.834 | 0.792 | 0.900 | 0.605 | 0.501 | 0.898 | 0.109 | 0.672 | 0.622 | 0.528 | 0.317 | 0.480 | 0.673 | 0.828 | - |
| congeners 6 | 288 | 256 | 0.795 | 0.732 | 0.858 | 0.596 | 0.509 | 0.863 | 0.139 | 0.655 | 0.583 | 0.555 | 0.279 | 0.508 | 0.677 | 0.851 | - |
| twin 0 | 288 | 259 | 0.916 | 0.851 | 0.919 | 0.753 | 0.497 | 0.859 | 0.141 | 0.672 | 0.605 | 0.509 | 0.276 | 0.477 | 0.683 | 0.856 | - |
| twin 1 | 288 | 254 | 0.871 | 0.743 | 0.828 | 0.525 | 0.517 | 0.824 | 0.176 | 0.654 | 0.599 | 0.569 | 0.304 | 0.503 | 0.657 | 0.804 | - |
| regime low | 192 | 168 | 0.787 | 0.718 | 0.865 | 0.573 | 0.506 | 0.814 | 0.186 | 0.675 | 0.613 | 0.596 | 0.310 | 0.520 | 0.666 | 0.778 | - |
| regime high | 192 | 166 | 0.854 | 0.836 | 0.952 | 0.627 | 0.486 | 0.905 | 0.097 | 0.703 | 0.624 | 0.491 | 0.301 | 0.470 | 0.687 | 0.879 | - |
| regime none | 192 | 179 | 0.806 | 0.730 | 0.833 | 0.594 | 0.507 | 0.809 | 0.193 | 0.617 | 0.557 | 0.539 | 0.280 | 0.497 | 0.675 | 0.851 | - |
| depth 1 | 159 | 135 | 0.806 | 0.757 | 0.884 | 0.611 | 0.521 | 0.852 | 0.148 | 0.675 | 0.599 | 0.547 | 0.301 | 0.494 | 0.686 | 0.830 | - |
| depth 2 | 136 | 122 | 0.815 | 0.754 | 0.874 | 0.590 | 0.494 | 0.828 | 0.175 | 0.624 | 0.555 | 0.497 | 0.248 | 0.450 | 0.632 | 0.810 | - |
| depth 4 | 138 | 127 | 0.829 | 0.791 | 0.889 | 0.572 | 0.481 | 0.836 | 0.164 | 0.672 | 0.629 | 0.563 | 0.309 | 0.511 | 0.676 | 0.838 | - |
| depth 0.5 | 143 | 129 | 0.809 | 0.734 | 0.860 | 0.602 | 0.506 | 0.843 | 0.158 | 0.675 | 0.613 | 0.550 | 0.323 | 0.509 | 0.678 | 0.840 | - |
| all, within identity bands of 0.005 | 576 | 513 | 0.620 | 0.579 | 0.723 | 0.646 | 0.688 | 0.717 | 0.284 | 0.635 | 0.637 | 0.583 | 0.332 | 0.494 | 0.577 | 0.698 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | ancestry_fixed_agreement | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 44 | 34 | 0.779 | 0.900 | 1.000 | 0.800 | 0.753 | 0.803 | 0.989 | - |
| Q_deep vs N_0 | congeners 6 | 57 | 41 | 0.819 | 0.940 | 1.000 | 0.707 | 0.683 | 0.730 | 0.998 | - |
| Q_deep vs N_0 | no twin | 48 | 38 | 0.997 | 0.993 | 0.985 | 0.796 | 0.735 | 0.820 | 0.988 | - |
| Q_deep vs N_0 | twin | 53 | 37 | 0.916 | 0.829 | 0.904 | 0.701 | 0.673 | 0.680 | 0.903 | - |
| Q_deep vs N_0 | no allele genomes | 27 | 31 | 0.827 | 0.892 | 0.916 | 0.500 | 0.500 | 0.500 | 0.916 | - |
| Q_deep vs N_0 | allele genomes | 74 | 44 | 0.807 | 0.939 | 0.948 | 0.887 | 0.815 | 0.904 | 0.972 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 11 | 7 | 0.870 | 0.987 | 1.000 | 0.688 | 0.636 | 1.000 | 1.000 | - |
| Q_deep vs N_0 | recombination none | 35 | 25 | 0.835 | 0.915 | 0.930 | 0.706 | 0.680 | 0.727 | 0.953 | - |
| Q_deep vs N_0 | recombination high | 29 | 23 | 0.838 | 0.970 | 1.000 | 0.783 | 0.700 | 0.790 | 0.999 | - |
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
| Q_deep vs N_0.8 | twin | 53 | 35 | 0.805 | 0.767 | 0.740 | 0.646 | 0.684 | 0.643 | 0.682 | - |
| Q_deep vs N_0.8 | no allele genomes | 27 | 23 | 0.713 | 0.760 | 0.681 | 0.500 | 0.500 | 0.500 | 0.681 | - |
| Q_deep vs N_0.8 | allele genomes | 74 | 42 | 0.645 | 0.822 | 0.714 | 0.774 | 0.801 | 0.771 | 0.798 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 11 | 9 | 0.596 | 0.758 | 0.737 | 0.707 | 0.747 | 0.808 | 0.838 | - |
| Q_deep vs N_0.8 | recombination none | 35 | 23 | 0.629 | 0.684 | 0.719 | 0.635 | 0.648 | 0.658 | 0.735 | - |
| Q_deep vs N_0.8 | recombination high | 29 | 17 | 0.625 | 0.982 | 0.878 | 0.765 | 0.722 | 0.759 | 0.824 | - |
| Q_deep vs N_0.95 | congeners 2 | 44 | 38 | 0.566 | 0.678 | 0.623 | 0.706 | 0.671 | 0.689 | 0.697 | - |
| Q_deep vs N_0.95 | congeners 6 | 57 | 35 | 0.595 | 0.607 | 0.591 | 0.618 | 0.665 | 0.601 | 0.592 | - |
| Q_deep vs N_0.95 | no twin | 48 | 37 | 0.651 | 0.680 | 0.617 | 0.681 | 0.683 | 0.634 | 0.655 | - |
| Q_deep vs N_0.95 | twin | 53 | 36 | 0.593 | 0.593 | 0.534 | 0.630 | 0.660 | 0.592 | 0.585 | - |
| Q_deep vs N_0.95 | no allele genomes | 27 | 23 | 0.652 | 0.651 | 0.578 | 0.500 | 0.500 | 0.500 | 0.578 | - |
| Q_deep vs N_0.95 | allele genomes | 74 | 50 | 0.574 | 0.623 | 0.559 | 0.737 | 0.786 | 0.667 | 0.666 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 11 | 7 | 0.623 | 0.688 | 0.636 | 0.727 | 0.727 | 0.831 | 0.831 | - |
| Q_deep vs N_0.95 | recombination none | 35 | 25 | 0.609 | 0.618 | 0.581 | 0.640 | 0.648 | 0.654 | 0.677 | - |
| Q_deep vs N_0.95 | recombination high | 29 | 26 | 0.576 | 0.768 | 0.672 | 0.714 | 0.718 | 0.677 | 0.699 | - |

## L2d: the gain of each group over identity and depth (logistic models, genera held out)

AUC of grouped-CV logistic models on the target rows; gain over identity, top identity and fragments with a 95% interval over resampled genera.
| set | model | auc_base | auc | gain | gain_lo | gain_hi |
|---|---|---|---|---|---|---|
| pe | identity, depth +ancestry | 0.780 | 0.822 | 0.043 | 0.022 | 0.064 |
| pe | identity, depth +alleles | 0.780 | 0.855 | 0.076 | 0.036 | 0.122 |
| pe | identity, depth +polymorphic | 0.780 | 0.799 | 0.019 | -0.014 | 0.047 |
| pe | identity, depth +all three | 0.780 | 0.895 | 0.116 | 0.067 | 0.162 |
| pe | identity, depth +ancestry, polymorphic | 0.780 | 0.870 | 0.090 | 0.060 | 0.119 |
| pe | identity, depth +oracle fixed agreement | 0.780 | 0.780 | 0.000 | -0.006 | 0.006 |
| pe, genera with alleles | identity, depth +ancestry | 0.784 | 0.828 | 0.044 | 0.015 | 0.070 |
| pe, genera with alleles | identity, depth +alleles | 0.784 | 0.927 | 0.143 | 0.084 | 0.200 |
| pe, genera with alleles | identity, depth +polymorphic | 0.784 | 0.891 | 0.107 | 0.071 | 0.144 |
| pe, genera with alleles | identity, depth +all three | 0.784 | 0.950 | 0.167 | 0.107 | 0.218 |
| pe, genera with alleles | identity, depth +ancestry, polymorphic | 0.784 | 0.915 | 0.131 | 0.092 | 0.172 |
| pe, genera with alleles | identity, depth +oracle fixed agreement | 0.784 | 0.858 | 0.074 | 0.047 | 0.096 |
| se | identity, depth +ancestry | 0.760 | 0.798 | 0.038 | 0.018 | 0.056 |
| se | identity, depth +alleles | 0.760 | 0.846 | 0.086 | 0.044 | 0.134 |
| se | identity, depth +polymorphic | 0.760 | 0.781 | 0.021 | -0.001 | 0.044 |
| se | identity, depth +all three | 0.760 | 0.884 | 0.124 | 0.081 | 0.169 |
| se | identity, depth +ancestry, polymorphic | 0.760 | 0.848 | 0.088 | 0.060 | 0.115 |
| se | identity, depth +oracle fixed agreement | 0.760 | 0.760 | 0.000 | 0.000 | 0.000 |
| se, genera with alleles | identity, depth +ancestry | 0.761 | 0.801 | 0.040 | 0.015 | 0.061 |
| se, genera with alleles | identity, depth +alleles | 0.761 | 0.914 | 0.152 | 0.092 | 0.210 |
| se, genera with alleles | identity, depth +polymorphic | 0.761 | 0.865 | 0.103 | 0.072 | 0.139 |
| se, genera with alleles | identity, depth +all three | 0.761 | 0.937 | 0.176 | 0.115 | 0.232 |
| se, genera with alleles | identity, depth +ancestry, polymorphic | 0.761 | 0.892 | 0.131 | 0.094 | 0.173 |
| se, genera with alleles | identity, depth +oracle fixed agreement | 0.761 | 0.761 | 0.000 | 0.000 | 0.000 |
| hifi | identity, depth +ancestry | 0.805 | 0.853 | 0.047 | 0.016 | 0.074 |
| hifi | identity, depth +alleles | 0.805 | 0.884 | 0.079 | 0.038 | 0.126 |
| hifi | identity, depth +polymorphic | 0.805 | 0.885 | 0.080 | 0.030 | 0.121 |
| hifi | identity, depth +all three | 0.805 | 0.942 | 0.137 | 0.089 | 0.179 |
| hifi | identity, depth +ancestry, polymorphic | 0.805 | 0.916 | 0.110 | 0.070 | 0.149 |
| hifi | identity, depth +oracle fixed agreement | 0.805 | 0.805 | 0.000 | 0.000 | 0.000 |
| hifi, genera with alleles | identity, depth +ancestry | 0.805 | 0.850 | 0.045 | 0.010 | 0.077 |
| hifi, genera with alleles | identity, depth +alleles | 0.805 | 0.953 | 0.147 | 0.083 | 0.199 |
| hifi, genera with alleles | identity, depth +polymorphic | 0.805 | 0.936 | 0.131 | 0.083 | 0.174 |
| hifi, genera with alleles | identity, depth +all three | 0.805 | 0.974 | 0.169 | 0.111 | 0.214 |
| hifi, genera with alleles | identity, depth +ancestry, polymorphic | 0.805 | 0.951 | 0.145 | 0.098 | 0.187 |
| hifi, genera with alleles | identity, depth +oracle fixed agreement | 0.805 | 0.805 | 0.000 | 0.000 | 0.000 |

## Sensitivity: by divergence of the lineage and by recombination with congeners

Per feature, the threshold 5% of all N rows exceed (the feature oriented by its AUC); rate_above: the share of Q rows above it (sensitivity) or of N rows (false-positive rate) in each bin. 'MRCA with nearest allele genome' is how recently the strain shared an ancestor with a strain the database knows; 'genes with a congener's segment' how much of it came from a congener.
| set | class | by | level | rows | identity | ancestry_agreement | allele_explained_share | polymorphic_known_share | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 278 | 0.295 | 0.392 | 0.403 | 0.295 | 0.457 | 0.000 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.078 | 0.125 | 0.078 | 0.094 | 0.078 | 0.000 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 573 | 0.351 | 0.527 | 0.000 | 0.000 | 0.092 | 0.000 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 809 | 0.387 | 0.524 | 0.740 | 0.472 | 0.696 | 0.000 |
| pe | Q | depth | 10 | 435 | 0.349 | 0.499 | 0.405 | 0.255 | 0.421 | 0.000 |
| pe | Q | depth | 100 | 426 | 0.329 | 0.469 | 0.423 | 0.293 | 0.404 | 0.000 |
| pe | Q | depth | 3 | 439 | 0.405 | 0.528 | 0.424 | 0.264 | 0.492 | 0.000 |
| pe | Q | depth | 30 | 424 | 0.309 | 0.458 | 0.410 | 0.278 | 0.417 | 0.000 |
| pe | Q | distance to rep | 0.005-0.015 | 790 | 0.363 | 0.516 | 0.473 | 0.272 | 0.420 | 0.000 |
| pe | Q | distance to rep | 0.015-0.03 | 433 | 0.051 | 0.330 | 0.527 | 0.367 | 0.365 | 0.000 |
| pe | Q | distance to rep | <0.005 | 298 | 0.977 | 0.966 | 0.074 | 0.070 | 0.735 | 0.000 |
| pe | Q | distance to rep | >=0.03 | 203 | 0.005 | 0.020 | 0.453 | 0.369 | 0.192 | 0.000 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 570 | 0.318 | 0.488 | 0.428 | 0.267 | 0.361 | 0.000 |
| pe | Q | genes with a congener's segment | <=0.1 | 454 | 0.326 | 0.478 | 0.410 | 0.229 | 0.463 | 0.000 |
| pe | Q | genes with a congener's segment | >0.3 | 251 | 0.303 | 0.514 | 0.502 | 0.307 | 0.335 | 0.000 |
| pe | Q | genes with a congener's segment | none | 449 | 0.437 | 0.488 | 0.356 | 0.305 | 0.552 | 0.000 |
| pe | Q | regime | high | 575 | 0.353 | 0.520 | 0.445 | 0.233 | 0.348 | 0.000 |
| pe | Q | regime | low | 574 | 0.286 | 0.456 | 0.425 | 0.267 | 0.432 | 0.000 |
| pe | Q | regime | none | 575 | 0.407 | 0.490 | 0.376 | 0.318 | 0.522 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 933 | 0.339 | 0.477 | 0.740 | 0.464 | 0.684 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 589 | 0.367 | 0.518 | 0.000 | 0.002 | 0.090 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 198 | 0.328 | 0.449 | 0.131 | 0.182 | 0.268 | 0.000 |
| pe | Q | role | Q_deep | 310 | 0.126 | 0.165 | 0.432 | 0.335 | 0.306 | 0.000 |
| pe | Q | role | Q_ils | 200 | 0.090 | 0.105 | 0.370 | 0.280 | 0.075 | 0.000 |
| pe | Q | role | Q_imp | 402 | 0.308 | 0.515 | 0.507 | 0.331 | 0.386 | 0.000 |
| pe | Q | role | Q_lone | 277 | 0.397 | 0.632 | 0.101 | 0.079 | 0.220 | 0.000 |
| pe | Q | role | Q_near | 334 | 0.395 | 0.569 | 0.763 | 0.464 | 0.832 | 0.000 |
| pe | Q | role | Q_rep | 201 | 0.886 | 0.990 | 0.104 | 0.000 | 0.716 | 0.000 |
| pe | Q | twin genus | 0 | 864 | 0.241 | 0.542 | 0.407 | 0.223 | 0.444 | 0.000 |
| pe | Q | twin genus | 1 | 860 | 0.457 | 0.436 | 0.423 | 0.322 | 0.423 | 0.000 |
| pe | N | depth | 10 | 382 | 0.052 | 0.026 | 0.058 | 0.071 | 0.047 | 0.000 |
| pe | N | depth | 100 | 377 | 0.011 | 0.029 | 0.013 | 0.005 | 0.016 | 0.000 |
| pe | N | depth | 3 | 354 | 0.124 | 0.124 | 0.093 | 0.105 | 0.113 | 0.000 |
| pe | N | depth | 30 | 381 | 0.018 | 0.026 | 0.039 | 0.024 | 0.029 | 0.000 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 861 | 0.043 | 0.056 | 0.042 | 0.052 | 0.048 | 0.000 |
| pe | N | genes with a segment of S or a congener | >0.3 | 266 | 0.045 | 0.015 | 0.049 | 0.064 | 0.019 | 0.000 |
| pe | N | genes with a segment of S or a congener | none | 367 | 0.071 | 0.063 | 0.071 | 0.035 | 0.079 | 0.000 |
| pe | N | regime | high | 492 | 0.063 | 0.057 | 0.053 | 0.061 | 0.035 | 0.000 |
| pe | N | regime | low | 496 | 0.024 | 0.044 | 0.036 | 0.056 | 0.054 | 0.000 |
| pe | N | regime | none | 506 | 0.063 | 0.049 | 0.061 | 0.034 | 0.061 | 0.000 |
| pe | N | role | N_0 | 196 | 0.010 | 0.000 | 0.031 | 0.056 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 203 | 0.025 | 0.010 | 0.020 | 0.039 | 0.015 | 0.000 |
| pe | N | role | N_0.8 | 215 | 0.042 | 0.028 | 0.047 | 0.051 | 0.028 | 0.000 |
| pe | N | role | N_0.95 | 227 | 0.079 | 0.115 | 0.070 | 0.031 | 0.079 | 0.000 |
| pe | N | role | N_1 | 228 | 0.110 | 0.154 | 0.105 | 0.053 | 0.184 | 0.000 |
| pe | N | role | N_imp0.2 | 213 | 0.033 | 0.014 | 0.033 | 0.066 | 0.014 | 0.000 |
| pe | N | role | N_imp0.4 | 212 | 0.042 | 0.014 | 0.038 | 0.057 | 0.014 | 0.000 |
| pe | N | twin genus | 0 | 759 | 0.003 | 0.063 | 0.026 | 0.028 | 0.043 | 0.000 |
| pe | N | twin genus | 1 | 735 | 0.099 | 0.037 | 0.075 | 0.073 | 0.057 | 0.000 |
| se | Q | MRCA with nearest allele genome | 0.002-0.01 | 275 | 0.280 | 0.364 | 0.404 | 0.240 | 0.076 |  |
| se | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.078 | 0.109 | 0.141 | 0.062 | 0.000 |  |
| se | Q | MRCA with nearest allele genome | no allele genome | 570 | 0.332 | 0.453 | 0.000 | 0.000 | 0.000 |  |
| se | Q | MRCA with nearest allele genome | recent (<0.002) | 805 | 0.354 | 0.477 | 0.717 | 0.405 | 0.130 |  |
| se | Q | depth | 10 | 435 | 0.340 | 0.444 | 0.393 | 0.216 | 0.120 |  |
| se | Q | depth | 100 | 426 | 0.300 | 0.401 | 0.427 | 0.242 | 0.038 |  |
| se | Q | depth | 3 | 429 | 0.366 | 0.485 | 0.382 | 0.231 | 0.051 |  |
| se | Q | depth | 30 | 424 | 0.290 | 0.417 | 0.425 | 0.236 | 0.085 |  |
| se | Q | distance to rep | 0.005-0.015 | 785 | 0.324 | 0.460 | 0.461 | 0.232 | 0.073 |  |
| se | Q | distance to rep | 0.015-0.03 | 429 | 0.054 | 0.249 | 0.503 | 0.312 | 0.054 |  |
| se | Q | distance to rep | <0.005 | 298 | 0.933 | 0.933 | 0.091 | 0.060 | 0.141 |  |
| se | Q | distance to rep | >=0.03 | 202 | 0.005 | 0.015 | 0.455 | 0.307 | 0.020 |  |
| se | Q | genes with a congener's segment | 0.1-0.3 | 567 | 0.325 | 0.429 | 0.420 | 0.208 | 0.067 |  |
| se | Q | genes with a congener's segment | <=0.1 | 453 | 0.272 | 0.406 | 0.400 | 0.185 | 0.071 |  |
| se | Q | genes with a congener's segment | >0.3 | 248 | 0.282 | 0.476 | 0.476 | 0.286 | 0.048 |  |
| se | Q | genes with a congener's segment | none | 446 | 0.401 | 0.457 | 0.359 | 0.276 | 0.099 |  |
| se | Q | regime | high | 571 | 0.340 | 0.473 | 0.427 | 0.186 | 0.060 |  |
| se | Q | regime | low | 573 | 0.258 | 0.389 | 0.414 | 0.215 | 0.072 |  |
| se | Q | regime | none | 570 | 0.375 | 0.449 | 0.379 | 0.293 | 0.089 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 928 | 0.314 | 0.430 | 0.719 | 0.399 | 0.128 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 586 | 0.340 | 0.447 | 0.000 | 0.000 | 0.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | partly | 196 | 0.316 | 0.429 | 0.153 | 0.133 | 0.036 |  |
| se | Q | role | Q_deep | 309 | 0.091 | 0.142 | 0.427 | 0.262 | 0.052 |  |
| se | Q | role | Q_ils | 197 | 0.071 | 0.091 | 0.335 | 0.254 | 0.015 |  |
| se | Q | role | Q_imp | 398 | 0.307 | 0.462 | 0.482 | 0.294 | 0.068 |  |
| se | Q | role | Q_lone | 276 | 0.366 | 0.496 | 0.109 | 0.072 | 0.040 |  |
| se | Q | role | Q_near | 333 | 0.354 | 0.523 | 0.751 | 0.381 | 0.156 |  |
| se | Q | role | Q_rep | 201 | 0.861 | 0.955 | 0.134 | 0.005 | 0.085 |  |
| se | Q | twin genus | 0 | 858 | 0.228 | 0.456 | 0.402 | 0.182 | 0.056 |  |
| se | Q | twin genus | 1 | 856 | 0.421 | 0.418 | 0.411 | 0.280 | 0.091 |  |
| se | N | depth | 10 | 382 | 0.050 | 0.045 | 0.071 | 0.071 | 0.016 |  |
| se | N | depth | 100 | 377 | 0.005 | 0.008 | 0.005 | 0.008 | 0.000 |  |
| se | N | depth | 3 | 351 | 0.131 | 0.145 | 0.094 | 0.100 | 0.023 |  |
| se | N | depth | 30 | 381 | 0.021 | 0.010 | 0.034 | 0.013 | 0.008 |  |
| se | N | genes with a segment of S or a congener | <=0.3 | 863 | 0.051 | 0.053 | 0.042 | 0.049 | 0.009 |  |
| se | N | genes with a segment of S or a congener | >0.3 | 264 | 0.042 | 0.030 | 0.038 | 0.053 | 0.008 |  |
| se | N | genes with a segment of S or a congener | none | 364 | 0.055 | 0.058 | 0.080 | 0.038 | 0.019 |  |
| se | N | regime | high | 489 | 0.063 | 0.053 | 0.051 | 0.057 | 0.008 |  |
| se | N | regime | low | 498 | 0.034 | 0.046 | 0.036 | 0.048 | 0.012 |  |
| se | N | regime | none | 504 | 0.054 | 0.052 | 0.063 | 0.036 | 0.014 |  |
| se | N | role | N_0 | 199 | 0.025 | 0.000 | 0.030 | 0.065 | 0.005 |  |
| se | N | role | N_0.5 | 202 | 0.030 | 0.010 | 0.025 | 0.035 | 0.005 |  |
| se | N | role | N_0.8 | 213 | 0.033 | 0.028 | 0.047 | 0.038 | 0.009 |  |
| se | N | role | N_0.95 | 226 | 0.066 | 0.106 | 0.080 | 0.049 | 0.018 |  |
| se | N | role | N_1 | 227 | 0.106 | 0.132 | 0.101 | 0.040 | 0.031 |  |
| se | N | role | N_imp0.2 | 212 | 0.038 | 0.038 | 0.028 | 0.052 | 0.005 |  |
| se | N | role | N_imp0.4 | 212 | 0.047 | 0.024 | 0.033 | 0.052 | 0.005 |  |
| se | N | twin genus | 0 | 756 | 0.003 | 0.038 | 0.022 | 0.033 | 0.003 |  |
| se | N | twin genus | 1 | 735 | 0.099 | 0.063 | 0.079 | 0.061 | 0.020 |  |
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
| hifi_base | target rows | 576 | 513 | 0.870 | 0.521 | 0.679 | 0.129 |  |
| hifi_base | hard: Q_deep vs N_0.8/0.95 | 101 | 138 | 0.682 | 0.257 | 0.455 | 0.239 |  |
| hifi_base | twin genera | 288 | 254 | 0.842 | 0.472 | 0.642 | 0.110 |  |
| hifi_base | no allele genomes | 192 | 172 | 0.848 | 0.510 | 0.661 | 0.157 |  |
| hifi_base | recombination high | 192 | 166 | 0.935 | 0.635 | 0.906 | 0.217 |  |
| hifi_base | all test rows | 896 | 936 | 0.948 | - | 0.794 | 0.071 | 0.850 |
| hifi_default | target rows | 576 | 513 | 0.940 | 0.797 | 0.845 | 0.090 |  |
| hifi_default | hard: Q_deep vs N_0.8/0.95 | 101 | 138 | 0.854 | 0.535 | 0.713 | 0.167 |  |
| hifi_default | twin genera | 288 | 254 | 0.928 | 0.795 | 0.875 | 0.126 |  |
| hifi_default | no allele genomes | 192 | 172 | 0.854 | 0.604 | 0.703 | 0.180 |  |
| hifi_default | recombination high | 192 | 166 | 0.974 | 0.875 | 0.906 | 0.114 |  |
| hifi_default | all test rows | 896 | 936 | 0.975 | - | 0.901 | 0.049 | 0.923 |
| pe_alleles | target rows | 1724 | 1494 | 0.926 | 0.745 | 0.811 | 0.108 |  |
| pe_alleles | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.810 | 0.510 | 0.642 | 0.147 |  |
| pe_alleles | twin genera | 860 | 735 | 0.934 | 0.743 | 0.823 | 0.128 |  |
| pe_alleles | no allele genomes | 573 | 496 | 0.827 | 0.504 | 0.625 | 0.147 |  |
| pe_alleles | recombination high | 575 | 492 | 0.941 | 0.758 | 0.828 | 0.112 |  |
| pe_alleles | all test rows | 2684 | 3730 | 0.978 | - | 0.879 | 0.044 | 0.906 |
| pe_ancestry+alleles | target rows | 1724 | 1494 | 0.916 | 0.742 | 0.796 | 0.101 |  |
| pe_ancestry+alleles | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.788 | 0.474 | 0.626 | 0.152 |  |
| pe_ancestry+alleles | twin genera | 860 | 735 | 0.924 | 0.750 | 0.793 | 0.094 |  |
| pe_ancestry+alleles | no allele genomes | 573 | 496 | 0.816 | 0.503 | 0.578 | 0.099 |  |
| pe_ancestry+alleles | recombination high | 575 | 492 | 0.935 | 0.755 | 0.810 | 0.100 |  |
| pe_ancestry+alleles | all test rows | 2684 | 3730 | 0.974 | - | 0.869 | 0.041 | 0.902 |
| pe_ancestry | target rows | 1724 | 1494 | 0.845 | 0.478 | 0.665 | 0.163 |  |
| pe_ancestry | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.647 | 0.181 | 0.455 | 0.269 |  |
| pe_ancestry | twin genera | 860 | 735 | 0.838 | 0.457 | 0.636 | 0.147 |  |
| pe_ancestry | no allele genomes | 573 | 496 | 0.821 | 0.504 | 0.667 | 0.173 |  |
| pe_ancestry | recombination high | 575 | 492 | 0.884 | 0.503 | 0.809 | 0.220 |  |
| pe_ancestry | all test rows | 2684 | 3730 | 0.957 | - | 0.784 | 0.066 | 0.836 |
| pe_base | target rows | 1724 | 1494 | 0.846 | 0.494 | 0.665 | 0.165 |  |
| pe_base | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.647 | 0.200 | 0.458 | 0.276 |  |
| pe_base | twin genera | 860 | 735 | 0.829 | 0.469 | 0.630 | 0.161 |  |
| pe_base | no allele genomes | 573 | 496 | 0.824 | 0.501 | 0.656 | 0.169 |  |
| pe_base | recombination high | 575 | 492 | 0.882 | 0.480 | 0.802 | 0.209 |  |
| pe_base | all test rows | 2684 | 3730 | 0.956 | - | 0.785 | 0.068 | 0.835 |
| pe_default | target rows | 1724 | 1494 | 0.923 | 0.738 | 0.809 | 0.106 |  |
| pe_default | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.803 | 0.474 | 0.632 | 0.145 |  |
| pe_default | twin genera | 860 | 735 | 0.935 | 0.752 | 0.819 | 0.113 |  |
| pe_default | no allele genomes | 573 | 496 | 0.817 | 0.475 | 0.621 | 0.143 |  |
| pe_default | recombination high | 575 | 492 | 0.937 | 0.744 | 0.831 | 0.122 |  |
| pe_default | all test rows | 2684 | 3730 | 0.976 | - | 0.877 | 0.043 | 0.906 |
| pe_default_on_noallele | target rows | 1714 | 1466 | 0.917 | 0.716 | 0.790 | 0.106 |  |
| pe_default_on_noallele | hard: Q_deep vs N_0.8/0.95 | 307 | 436 | 0.801 | 0.404 | 0.612 | 0.147 |  |
| pe_default_on_noallele | twin genera | 850 | 709 | 0.923 | 0.708 | 0.788 | 0.107 |  |
| pe_default_on_noallele | no allele genomes | 573 | 496 | 0.818 | 0.476 | 0.618 | 0.133 |  |
| pe_default_on_noallele | recombination high | 572 | 486 | 0.933 | 0.745 | 0.811 | 0.113 |  |
| pe_default_on_noallele | all test rows | 2674 | 3883 | 0.976 | - | 0.865 | 0.041 | 0.899 |
| pe_default_on_shuffled | target rows | 1716 | 1468 | 0.830 | 0.480 | 0.579 | 0.107 |  |
| pe_default_on_shuffled | hard: Q_deep vs N_0.8/0.95 | 307 | 437 | 0.629 | 0.111 | 0.319 | 0.151 |  |
| pe_default_on_shuffled | twin genera | 852 | 711 | 0.816 | 0.481 | 0.560 | 0.114 |  |
| pe_default_on_shuffled | no allele genomes | 573 | 496 | 0.818 | 0.480 | 0.618 | 0.135 |  |
| pe_default_on_shuffled | recombination high | 573 | 487 | 0.873 | 0.553 | 0.635 | 0.115 |  |
| pe_default_on_shuffled | all test rows | 2676 | 3868 | 0.950 | - | 0.730 | 0.042 | 0.816 |

The models' calls on the target rows by bin: 'called' the share called at the model's knob (Q: sensitivity; N: false-positive rate), 'at 5% FPR' the share above the score 5% of the model's N rows exceed.
| set | class | by | level | rows | base called | ancestry called | alleles called | default called | base at 5% FPR | ancestry at 5% FPR | alleles at 5% FPR | default at 5% FPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 278 | 0.608 | 0.622 | 0.817 | 0.824 | 0.450 | 0.442 | 0.723 | 0.709 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.312 | 0.266 | 0.297 | 0.250 | 0.141 | 0.125 | 0.188 | 0.188 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 573 | 0.656 | 0.667 | 0.625 | 0.621 | 0.496 | 0.483 | 0.506 | 0.501 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 809 | 0.719 | 0.710 | 0.983 | 0.980 | 0.536 | 0.514 | 0.967 | 0.960 |
| pe | Q | depth | 10 | 435 | 0.694 | 0.692 | 0.821 | 0.821 | 0.501 | 0.480 | 0.752 | 0.761 |
| pe | Q | depth | 100 | 426 | 0.674 | 0.671 | 0.819 | 0.819 | 0.538 | 0.535 | 0.772 | 0.765 |
| pe | Q | depth | 3 | 439 | 0.624 | 0.622 | 0.768 | 0.765 | 0.437 | 0.392 | 0.688 | 0.658 |
| pe | Q | depth | 30 | 424 | 0.670 | 0.675 | 0.840 | 0.830 | 0.502 | 0.507 | 0.771 | 0.771 |
| pe | Q | distance to rep | 0.005-0.015 | 790 | 0.711 | 0.716 | 0.859 | 0.856 | 0.509 | 0.481 | 0.782 | 0.768 |
| pe | Q | distance to rep | 0.015-0.03 | 433 | 0.573 | 0.570 | 0.776 | 0.774 | 0.342 | 0.342 | 0.672 | 0.670 |
| pe | Q | distance to rep | <0.005 | 298 | 0.987 | 0.983 | 0.990 | 0.990 | 0.960 | 0.946 | 0.973 | 0.973 |
| pe | Q | distance to rep | >=0.03 | 203 | 0.212 | 0.197 | 0.438 | 0.433 | 0.079 | 0.069 | 0.424 | 0.424 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 570 | 0.749 | 0.751 | 0.825 | 0.832 | 0.507 | 0.516 | 0.747 | 0.749 |
| pe | Q | genes with a congener's segment | <=0.1 | 454 | 0.586 | 0.590 | 0.698 | 0.703 | 0.465 | 0.436 | 0.648 | 0.659 |
| pe | Q | genes with a congener's segment | >0.3 | 251 | 0.701 | 0.693 | 0.853 | 0.845 | 0.510 | 0.486 | 0.801 | 0.785 |
| pe | Q | genes with a congener's segment | none | 449 | 0.619 | 0.615 | 0.886 | 0.866 | 0.499 | 0.468 | 0.811 | 0.780 |
| pe | Q | regime | high | 575 | 0.802 | 0.809 | 0.828 | 0.831 | 0.558 | 0.558 | 0.744 | 0.750 |
| pe | Q | regime | low | 574 | 0.575 | 0.573 | 0.707 | 0.711 | 0.439 | 0.416 | 0.660 | 0.662 |
| pe | Q | regime | none | 575 | 0.619 | 0.612 | 0.899 | 0.883 | 0.485 | 0.459 | 0.831 | 0.803 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 933 | 0.675 | 0.674 | 0.972 | 0.975 | 0.494 | 0.480 | 0.950 | 0.949 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 589 | 0.662 | 0.672 | 0.632 | 0.625 | 0.504 | 0.487 | 0.513 | 0.497 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 198 | 0.621 | 0.591 | 0.586 | 0.566 | 0.455 | 0.429 | 0.470 | 0.460 |
| pe | Q | role | Q_deep | 310 | 0.458 | 0.455 | 0.642 | 0.632 | 0.239 | 0.206 | 0.539 | 0.516 |
| pe | Q | role | Q_ils | 200 | 0.415 | 0.430 | 0.580 | 0.545 | 0.145 | 0.130 | 0.430 | 0.425 |
| pe | Q | role | Q_imp | 402 | 0.687 | 0.682 | 0.866 | 0.866 | 0.490 | 0.485 | 0.823 | 0.816 |
| pe | Q | role | Q_lone | 277 | 0.740 | 0.736 | 0.755 | 0.762 | 0.603 | 0.581 | 0.661 | 0.653 |
| pe | Q | role | Q_near | 334 | 0.728 | 0.731 | 0.991 | 0.991 | 0.584 | 0.557 | 0.979 | 0.973 |
| pe | Q | role | Q_rep | 201 | 0.985 | 0.980 | 0.975 | 0.990 | 0.945 | 0.955 | 0.950 | 0.965 |
| pe | Q | twin genus | 0 | 864 | 0.700 | 0.693 | 0.800 | 0.799 | 0.537 | 0.527 | 0.740 | 0.731 |
| pe | Q | twin genus | 1 | 860 | 0.630 | 0.636 | 0.823 | 0.819 | 0.451 | 0.429 | 0.751 | 0.745 |
| pe | N | depth | 10 | 382 | 0.175 | 0.165 | 0.126 | 0.118 | 0.068 | 0.050 | 0.073 | 0.068 |
| pe | N | depth | 100 | 377 | 0.154 | 0.151 | 0.085 | 0.093 | 0.040 | 0.050 | 0.029 | 0.029 |
| pe | N | depth | 3 | 354 | 0.169 | 0.169 | 0.124 | 0.116 | 0.048 | 0.051 | 0.065 | 0.059 |
| pe | N | depth | 30 | 381 | 0.160 | 0.168 | 0.100 | 0.100 | 0.045 | 0.050 | 0.034 | 0.045 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 861 | 0.180 | 0.187 | 0.096 | 0.102 | 0.063 | 0.064 | 0.039 | 0.045 |
| pe | N | genes with a segment of S or a congener | >0.3 | 266 | 0.075 | 0.068 | 0.071 | 0.071 | 0.011 | 0.008 | 0.034 | 0.030 |
| pe | N | genes with a segment of S or a congener | none | 367 | 0.193 | 0.177 | 0.163 | 0.142 | 0.049 | 0.049 | 0.087 | 0.076 |
| pe | N | regime | high | 492 | 0.209 | 0.220 | 0.112 | 0.122 | 0.079 | 0.083 | 0.041 | 0.055 |
| pe | N | regime | low | 496 | 0.121 | 0.119 | 0.069 | 0.073 | 0.030 | 0.032 | 0.034 | 0.030 |
| pe | N | regime | none | 506 | 0.164 | 0.152 | 0.144 | 0.125 | 0.042 | 0.036 | 0.075 | 0.065 |
| pe | N | role | N_0 | 196 | 0.000 | 0.000 | 0.005 | 0.005 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 203 | 0.015 | 0.025 | 0.015 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.8 | 215 | 0.195 | 0.186 | 0.088 | 0.088 | 0.028 | 0.033 | 0.033 | 0.033 |
| pe | N | role | N_0.95 | 227 | 0.352 | 0.348 | 0.203 | 0.198 | 0.101 | 0.106 | 0.097 | 0.110 |
| pe | N | role | N_1 | 228 | 0.417 | 0.412 | 0.289 | 0.298 | 0.184 | 0.184 | 0.149 | 0.145 |
| pe | N | role | N_imp0.2 | 213 | 0.038 | 0.042 | 0.052 | 0.056 | 0.009 | 0.005 | 0.028 | 0.023 |
| pe | N | role | N_imp0.4 | 212 | 0.085 | 0.080 | 0.075 | 0.066 | 0.009 | 0.005 | 0.028 | 0.024 |
| pe | N | twin genus | 0 | 759 | 0.169 | 0.179 | 0.090 | 0.100 | 0.054 | 0.059 | 0.047 | 0.051 |
| pe | N | twin genus | 1 | 735 | 0.161 | 0.147 | 0.128 | 0.113 | 0.046 | 0.041 | 0.053 | 0.049 |
| hifi | Q | MRCA with nearest allele genome | 0.002-0.01 | 97 | 0.536 |  |  | 0.804 | 0.340 |  |  | 0.753 |
| hifi | Q | MRCA with nearest allele genome | long ago (>=0.01) | 20 | 0.600 |  |  | 0.400 | 0.350 |  |  | 0.200 |
| hifi | Q | MRCA with nearest allele genome | no allele genome | 192 | 0.661 |  |  | 0.703 | 0.526 |  |  | 0.609 |
| hifi | Q | MRCA with nearest allele genome | recent (<0.002) | 267 | 0.749 |  |  | 0.996 | 0.596 |  |  | 0.993 |
| hifi | Q | depth | 0.5 | 143 | 0.650 |  |  | 0.832 | 0.497 |  |  | 0.783 |
| hifi | Q | depth | 1 | 159 | 0.667 |  |  | 0.824 | 0.509 |  |  | 0.780 |
| hifi | Q | depth | 2 | 136 | 0.706 |  |  | 0.875 | 0.581 |  |  | 0.824 |
| hifi | Q | depth | 4 | 138 | 0.696 |  |  | 0.855 | 0.500 |  |  | 0.804 |
| hifi | Q | distance to rep | 0.005-0.015 | 274 | 0.734 |  |  | 0.923 | 0.529 |  |  | 0.858 |
| hifi | Q | distance to rep | 0.015-0.03 | 148 | 0.547 |  |  | 0.743 | 0.351 |  |  | 0.676 |
| hifi | Q | distance to rep | <0.005 | 94 | 1.000 |  |  | 1.000 | 1.000 |  |  | 1.000 |
| hifi | Q | distance to rep | >=0.03 | 60 | 0.250 |  |  | 0.500 | 0.150 |  |  | 0.500 |
| hifi | Q | genes with a congener's segment | 0.1-0.3 | 190 | 0.811 |  |  | 0.895 | 0.605 |  |  | 0.853 |
| hifi | Q | genes with a congener's segment | <=0.1 | 152 | 0.507 |  |  | 0.684 | 0.395 |  |  | 0.625 |
| hifi | Q | genes with a congener's segment | >0.3 | 82 | 0.829 |  |  | 0.902 | 0.622 |  |  | 0.878 |
| hifi | Q | genes with a congener's segment | none | 152 | 0.605 |  |  | 0.914 | 0.487 |  |  | 0.855 |
| hifi | Q | regime | high | 192 | 0.906 |  |  | 0.906 | 0.693 |  |  | 0.870 |
| hifi | Q | regime | low | 192 | 0.505 |  |  | 0.703 | 0.385 |  |  | 0.646 |
| hifi | Q | regime | none | 192 | 0.625 |  |  | 0.927 | 0.484 |  |  | 0.875 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 314 | 0.688 |  |  | 0.987 | 0.519 |  |  | 0.984 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | no such site | 2 | 1.000 |  |  | 1.000 | 1.000 |  |  | 1.000 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 198 | 0.672 |  |  | 0.712 | 0.535 |  |  | 0.621 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | partly | 62 | 0.645 |  |  | 0.548 | 0.468 |  |  | 0.403 |
| hifi | Q | role | Q_deep | 101 | 0.455 |  |  | 0.713 | 0.267 |  |  | 0.594 |
| hifi | Q | role | Q_ils | 74 | 0.419 |  |  | 0.676 | 0.149 |  |  | 0.608 |
| hifi | Q | role | Q_imp | 136 | 0.743 |  |  | 0.897 | 0.574 |  |  | 0.875 |
| hifi | Q | role | Q_lone | 92 | 0.728 |  |  | 0.804 | 0.576 |  |  | 0.728 |
| hifi | Q | role | Q_near | 100 | 0.740 |  |  | 0.990 | 0.620 |  |  | 0.980 |
| hifi | Q | role | Q_rep | 73 | 0.986 |  |  | 0.959 | 0.945 |  |  | 0.959 |
| hifi | Q | twin genus | 0 | 288 | 0.715 |  |  | 0.816 | 0.594 |  |  | 0.771 |
| hifi | Q | twin genus | 1 | 288 | 0.642 |  |  | 0.875 | 0.448 |  |  | 0.823 |
| hifi | N | depth | 0.5 | 129 | 0.116 |  |  | 0.070 | 0.039 |  |  | 0.039 |
| hifi | N | depth | 1 | 135 | 0.111 |  |  | 0.074 | 0.044 |  |  | 0.037 |
| hifi | N | depth | 2 | 122 | 0.164 |  |  | 0.123 | 0.074 |  |  | 0.074 |
| hifi | N | depth | 4 | 127 | 0.126 |  |  | 0.094 | 0.047 |  |  | 0.055 |
| hifi | N | genes with a segment of S or a congener | <=0.3 | 289 | 0.156 |  |  | 0.076 | 0.045 |  |  | 0.024 |
| hifi | N | genes with a segment of S or a congener | >0.3 | 98 | 0.061 |  |  | 0.051 | 0.041 |  |  | 0.020 |
| hifi | N | genes with a segment of S or a congener | none | 126 | 0.119 |  |  | 0.151 | 0.071 |  |  | 0.135 |
| hifi | N | regime | high | 166 | 0.217 |  |  | 0.114 | 0.078 |  |  | 0.036 |
| hifi | N | regime | low | 168 | 0.071 |  |  | 0.030 | 0.018 |  |  | 0.006 |
| hifi | N | regime | none | 179 | 0.101 |  |  | 0.123 | 0.056 |  |  | 0.106 |
| hifi | N | role | N_0 | 75 | 0.000 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_0.5 | 74 | 0.014 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_0.8 | 65 | 0.108 |  |  | 0.185 | 0.000 |  |  | 0.062 |
| hifi | N | role | N_0.95 | 73 | 0.356 |  |  | 0.151 | 0.123 |  |  | 0.123 |
| hifi | N | role | N_1 | 72 | 0.333 |  |  | 0.236 | 0.167 |  |  | 0.153 |
| hifi | N | role | N_imp0.2 | 77 | 0.026 |  |  | 0.013 | 0.013 |  |  | 0.000 |
| hifi | N | role | N_imp0.4 | 77 | 0.078 |  |  | 0.065 | 0.052 |  |  | 0.026 |
| hifi | N | twin genus | 0 | 259 | 0.147 |  |  | 0.054 | 0.062 |  |  | 0.019 |
| hifi | N | twin genus | 1 | 254 | 0.110 |  |  | 0.126 | 0.039 |  |  | 0.083 |
