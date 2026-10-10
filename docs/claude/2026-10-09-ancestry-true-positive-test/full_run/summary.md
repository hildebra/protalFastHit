# Ancestry true-positive test: /home/falk/atp_full

Worlds: test (seed 11), train (seed 23); 36 test genera. See docs/claude/2026-10-09-ancestry-true-positive-test/README.md for the design and the predictions.

## Verdicts

| hypothesis | prediction | outcome | verdict |
|---|---|---|---|
| H1a table | stored alleles are the allele genomes' differences (>= 0.95 exact) | min share exact 1.000 | PASS |
| H3 sites | consensus sites = stem + representative's lineage (no recombination, >= 3 congeners): >= 0.98 of them, recall >= 0.95 | share 0.965, recall 0.798 | FAIL |
| H3 fallback | with 2 congeners the nearest congener's own derived states enter (share in history clearly below 1) | mean share 0.502 | PASS |
| H7 novel reads | the allele scores draw some of a novel species' reads onto the target (it carries the alleles' ancestral base at the representative's own sites) | 6072 of 115881 reads (0.052) | INFO |
| H7 settling | the strains' reads the allele scores move go to their species (>= 0.9) | 1.000 of 3187 moved | PASS |
| H1b/H2b oracle | protal's ancestry, allele and polymorphic features equal the oracle's (error-free reads, every row within 1e-6) | worst share within 1e-6: 1.000 | PASS |
| H1 leak | an allele genome's own reads: allele_explained_share >= 0.95, polymorphic_novel_share <= 0.02 | 0.965, 0.000 | PASS |
| H2 calibration | planted A_a: ancestry_agreement = a (+-0.03) | worst |diff| 0.009 | PASS |
| H1 calibration | planted P_f: allele_explained_share rises with f (a step near 0.5: a read is explained once it shares more than half of the allele's edits in its span) | f 0: 0.000, f 0.25: 0.052, f 0.5: 0.285, f 0.75: 0.551, f 1: 0.691 | PASS |
| H4 ancestry, clean strata | the same with >= 3 congeners and no recombination (the consensus applies) | b 0: 1.000, b 0.5: 0.901 | INFO |
| H5 fixed sites | deep strains with their representative-lineage sites covered (>= 0.8) against late novel species (b 0.8, 0.95), genera with alleles: ancestry_fixed_gain adds to ancestry_agreement (CI above 0) | AUC 0.647 -> 0.879, gain 0.232 [0.149, 0.312]; the true fixed-site agreement alone 0.834 (160 Q, 297 N rows) | PASS |
| H4 ancestry | ancestry_agreement separates deep strains from novel species that left the stem early (b <= 0.5): AUC >= 0.9 | b 0: 0.904, b 0.5: 0.807 | FAIL |
| H6 confound | allele_explained_share favours strains with a stored relative (Q_near) over late twins, but not lone strains | Q_near 0.904, Q_lone 0.320 | PASS |
| H9 decoy | the decoy (a tip of the deep strains' clade labelled a species) against Q_deep: no feature separates (AUC 0.5 +- 0.1) | largest |AUC - 0.5|: top_identity 0.463 (310 Q, 200 decoy rows) | PASS |
| H9b crown lineage | N_1, a novel lineage leaving at S's crown, against Q_deep: the ancestry features cannot tell them apart; the allele features can where Q_deep shares its side of the crown with allele genomes | ancestry_agreement 0.523, allele_explained_share 0.640 | INFO |
| H4/H5 gain | the three groups add separation over identity and depth (pe, CI above 0) | 0.106 [0.058, 0.156] | PASS |
| H5 ceiling | the true fixed-site agreement alone (the oracle's) adds separation; compare with '+all three' for the gap | 0.000 [-0.006, 0.006] | INFO |
| H5 ceiling, genera with alleles | the true fixed-site agreement (the oracle's) over identity and depth, where the fixed sites exist; compare with '+all three' | 0.070 [0.043, 0.092] | INFO |
| H8 calls pe target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.073 [0.040, 0.112] | PASS |
| H8 calls pe hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.150 [0.077, 0.224] | PASS |
| H8 calls hifi target rows | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.067 [0.040, 0.099] | PASS |
| H8 calls hifi hard | the default set's model separates better than base (AUC gain, 95% CI by genus above 0) | 0.162 [0.082, 0.234] | PASS |
| H9 shuffled | shuffled alleles cost the default model separation (the alleles carry information) | AUC 0.918 -> 0.829 | PASS |

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

Per strain with an allele genome in its species: the share of genes in which its nearest allele genome is among the stored alleles (K = 4 farthest first; the build rejects alleles as far as the nearest congener's copy).
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
| GCA_998007004.1 | 7 | Q_deep | 8 | 0 | 0.011217 | 0.642 |
| GCA_998007005.1 | 7 | Q_deep | 8 | 0 | 0.011217 | 0.642 |
| GCA_998007017.1 | 7 | Q_ils | 8 | 0 | 0.011217 | 0.642 |
| GCA_998007013.1 | 7 | Q_imp0.2 | 8 | 0 | 0.000103 | 0.642 |
| GCA_998007015.1 | 7 | Q_imp0.4 | 8 | 0 | 0.000103 | 0.642 |
| GCA_998007001.1 | 7 | Q_near | 8 | 0 | 0.000103 | 0.642 |
| GCA_998007002.1 | 7 | Q_near | 8 | 0 | 2.3e-05 | 0.500 |
| GCA_998007003.1 | 7 | Q_near | 8 | 0 | 6.1e-05 | 0.450 |
| GCA_998007006.1 | 7 | Q_rep | 8 | 0 | 0.000784 | 0.450 |
| GCA_998008004.1 | 8 | Q_deep | 8 | 0 | 0.008927 | 0.567 |
| GCA_998008005.1 | 8 | Q_deep | 8 | 0 | 0.008927 | 0.567 |
| GCA_998008017.1 | 8 | Q_ils | 8 | 0 | 0.008927 | 0.567 |
| GCA_998008013.1 | 8 | Q_imp0.2 | 8 | 0 | 8e-06 | 0.567 |
(171 more rows in l0_nearest_allele_stored.tsv)

The build's line: Strain alleles: 7640 alleles (116863 edits) of 2776 gene copies of 24 species, up to 4 each (of 38880 full-reference copies of the database's species and genes: 0 of genomes outside --allele_genome_share 1, 26017 identical to the representative's, 2245 repeated, 45 covering less than 0.5 of it or more than 0.1 apart, 361 as far as the nearest congener's copy or farther; 10212 alleles offered): /home/falk/atp_full/test/db/strain_alleles.tsv

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
| 19 | 6 | 0 | none | 3861 | 0.978 | 0.827 | 0.819 | 10 | 691 | 85 | 65 | 20 |
| 20 | 6 | 0 | low | 4376 | 0.942 | 0.812 | 0.807 | 21 | 803 | 110 | 224 | 29 |
| 21 | 6 | 0 | high | 3414 | 0.974 | 0.707 | 0.677 | 26 | 1210 | 98 | 77 | 13 |
| 22 | 6 | 0 | none | 2891 | 0.973 | 0.855 | 0.849 |  | 422 | 67 | 59 | 20 |
| 23 | 6 | 0 | low | 4030 | 0.949 | 0.793 | 0.771 | 7 | 903 | 115 | 183 | 23 |
| 24 | 6 | 0 | high | 3402 | 0.958 | 0.736 | 0.662 | 18 | 1144 | 79 | 128 | 15 |
| 25 | 6 | 0 | none | 4488 | 0.980 | 0.861 | 0.846 | 22 | 591 | 121 | 74 | 16 |
| 26 | 6 | 0 | low | 4266 | 0.970 | 0.806 | 0.790 | 11 | 857 | 109 | 113 | 17 |
| 27 | 6 | 0 | high | 3723 | 0.930 | 0.716 | 0.649 | 23 | 1298 | 107 | 252 | 8 |
| 28 | 6 | 1 | none | 1440 | 0.967 | 0.826 | 0.813 |  | 277 | 26 | 32 | 16 |
| 29 | 6 | 1 | low | 1458 | 0.945 | 0.816 | 0.832 | 1 | 260 | 23 | 69 | 11 |
| 30 | 6 | 1 | high | 1427 | 0.869 | 0.759 | 0.503 | 1 | 640 | 23 | 177 | 10 |
| 31 | 6 | 1 | none | 1288 | 0.982 | 0.799 | 0.798 | 1 | 290 | 28 | 21 | 2 |
| 32 | 6 | 1 | low | 1335 | 0.961 | 0.839 | 0.816 | 1 | 248 | 23 | 41 | 11 |
| 33 | 6 | 1 | high | 1666 | 0.837 | 0.769 | 0.682 | 8 | 452 | 37 | 252 | 20 |
| 34 | 6 | 1 | none | 1366 | 0.965 | 0.799 | 0.805 | 1 | 284 | 40 | 37 | 11 |
| 35 | 6 | 1 | low | 1565 | 0.910 | 0.796 | 0.825 |  | 308 | 35 | 122 | 19 |
| 36 | 6 | 1 | high | 1202 | 0.927 | 0.719 | 0.450 | 4 | 731 | 22 | 76 | 12 |

## L1: the focus genomes' reads, with and without the allele scores

The paired-end focus reads (best record of each mate): the share on the target without and with the allele scores, and the reads the scores moved (Q: the target is right; N: no species is).
| class | twin | reads | on_target_noallele | on_target_default | moved | moved_to_target | share_moves_to_target |
|---|---|---|---|---|---|---|---|
| Q | 0 | 66490 | 0.939 | 0.947 | 549 | 549 | 1.000 |
| Q | 1 | 65713 | 0.887 | 0.927 | 2638 | 2638 | 1.000 |
| N | 0 | 57189 | 0.833 | 0.862 | 1650 | 1650 | 1.000 |
| N | 1 | 58692 | 0.757 | 0.833 | 4422 | 4422 | 1.000 |
| class | twin | from | to | reads |
|---|---|---|---|---|
| D | 0 | congener | target | 9 |
| D | 0 | none | target | 1 |
| D | 0 | sister | target | 102 |
| D | 1 | congener | target | 10 |
| D | 1 | sister | target | 440 |
| N | 0 | congener | target | 157 |
| N | 0 | sister | target | 1493 |
| N | 1 | congener | target | 87 |
| N | 1 | sister | target | 4335 |
| Q | 0 | congener | target | 97 |
| Q | 0 | none | target | 2 |
| Q | 0 | sister | target | 450 |
| Q | 1 | congener | target | 107 |
| Q | 1 | sister | target | 2531 |
The protal log's 'strain alleles:' lines (every sample of the default run): 332282 unsure reads took the shifts, 25787 of them to another species.

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
| default | exact | allele_copy_share | 864 | 0.000 | 1.000 |
| default | exact | allele_sites_per_kb | 864 | 0.000 | 1.000 |
| default | exact | ancestry_fixed_share | 864 | 0.000 | 1.000 |
| default | pe | ancestry_sites_per_record | 3417 | 0.000 | 1.000 |
| default | pe | ancestry_agreement | 3417 | 0.000 | 1.000 |
| default | pe | ancestry_congener_share | 3417 | 0.000 | 1.000 |
| default | pe | ancestry_indel_sites_per_record | 3417 | 0.000 | 1.000 |
| default | pe | ancestry_indel_congener_share | 3417 | 0.000 | 1.000 |
| default | pe | allele_explained_share | 3417 | 0.000 | 1.000 |
| default | pe | allele_identity_gain | 3417 | 0.000 | 1.000 |
| default | pe | polymorphic_known_share | 3417 | 0.000 | 1.000 |
| default | pe | polymorphic_novel_share | 3417 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_gain | 3417 | 0.000 | 1.000 |
| default | pe | allele_copy_share | 3417 | 0.000 | 1.000 |
| default | pe | allele_sites_per_kb | 3417 | 0.000 | 1.000 |
| default | pe | ancestry_fixed_share | 3417 | 0.000 | 1.000 |
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
| shuffled | exact | allele_copy_share | 864 | 0.000 | 1.000 |
| shuffled | exact | allele_sites_per_kb | 864 | 0.000 | 1.000 |
| shuffled | exact | ancestry_fixed_share | 864 | 0.000 | 1.000 |
| noallele | exact | ancestry_sites_per_record | 864 | 0.000 | 1.000 |
(12 more rows in l2_oracle.tsv)

## L2b: the planted genomes and the leaky control (error-free reads)

Mean features of the planted genomes (P_f: a share f of the farthest allele genome's differences; A_a: the congeners' base at a share 1 - a of the consensus sites) and of Q_leak (an allele genome itself). truth_agreement: the genome's agreement at the consensus sites by position.
| role | rows | truth_agreement | ancestry_agreement | allele_explained_share | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | identity |
|---|---|---|---|---|---|---|---|---|
| A_a0 | 7 | 0.001 | 0.001 | 0.173 | 0.324 | 0.011 | -0.000 | 0.978 |
| A_a0.25 | 7 | 0.250 | 0.259 | 0.090 | 0.229 | 0.011 | 0.004 | 0.979 |
| A_a0.5 | 6 | 0.500 | 0.500 | 0.049 | 0.170 | 0.005 | 0.002 | 0.980 |
| A_a0.75 | 7 | 0.750 | 0.759 | 0.018 | 0.079 | 0.003 | -0.002 | 0.980 |
| A_a1 | 6 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0 | 7 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.980 |
| P_f0.25 | 6 | 0.867 | 0.924 | 0.052 | 0.202 | 0.000 | 0.075 | 0.980 |
| P_f0.5 | 7 | 0.942 | 0.824 | 0.285 | 0.402 | 0.000 | 0.171 | 0.980 |
| P_f0.75 | 6 | 0.828 | 0.778 | 0.551 | 0.589 | 0.000 | 0.221 | 0.979 |
| P_f1 | 6 | 0.465 | 0.718 | 0.691 | 0.800 | 0.000 | 0.281 | 0.977 |
| Q_leak | 110 | 0.397 | 0.803 | 0.965 | 0.510 | 0.000 | 0.191 | 0.989 |

## L2c: separation, pe (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1724 | 1493 | 0.780 | 0.679 | 0.798 | 0.551 | 0.502 | 0.810 | 0.193 | 0.620 | 0.548 | 0.533 | 0.422 | 0.489 | 0.817 | 0.649 |
| Q_deep vs N_0 | 310 | 195 | 0.798 | 0.780 | 0.881 | 0.518 | 0.558 | 0.904 | 0.100 | 0.701 | 0.645 | 0.626 | 0.450 | 0.692 | 0.923 | 0.742 |
| Q_deep vs N_0.5 | 310 | 203 | 0.717 | 0.666 | 0.808 | 0.531 | 0.510 | 0.807 | 0.197 | 0.681 | 0.630 | 0.624 | 0.445 | 0.592 | 0.821 | 0.706 |
| Q_deep vs N_0.8 | 310 | 215 | 0.614 | 0.612 | 0.678 | 0.524 | 0.508 | 0.679 | 0.327 | 0.652 | 0.617 | 0.606 | 0.467 | 0.549 | 0.713 | 0.654 |
| Q_deep vs N_0.95 | 310 | 227 | 0.532 | 0.525 | 0.535 | 0.522 | 0.519 | 0.573 | 0.427 | 0.649 | 0.630 | 0.630 | 0.482 | 0.538 | 0.619 | 0.592 |
| Q_deep vs N_1 | 310 | 228 | 0.484 | 0.504 | 0.482 | 0.507 | 0.510 | 0.523 | 0.478 | 0.640 | 0.630 | 0.623 | 0.476 | 0.532 | 0.567 | 0.551 |
| Q_near vs N_0.95 | 334 | 227 | 0.747 | 0.636 | 0.722 | 0.556 | 0.498 | 0.759 | 0.243 | 0.904 | 0.778 | 0.753 | 0.458 | 0.591 | 0.917 | 0.878 |
| Q_lone vs N_0.95 | 277 | 227 | 0.762 | 0.651 | 0.769 | 0.561 | 0.493 | 0.767 | 0.233 | 0.320 | 0.298 | 0.300 | 0.360 | 0.285 | 0.579 | 0.342 |
| Q_rep vs N_0 | 201 | 195 | 0.993 | 0.890 | 0.987 | 0.601 | 0.527 | 1.000 | 0.000 | 0.451 | 0.344 | 0.308 | 0.397 | 0.405 | 1.000 | 0.692 |
| Q_imp vs N (b <= 0.95) | 402 | 840 | 0.820 | 0.734 | 0.840 | 0.560 | 0.506 | 0.852 | 0.152 | 0.661 | 0.581 | 0.560 | 0.393 | 0.497 | 0.852 | 0.668 |
| Q (no import) vs N_imp | 1322 | 425 | 0.801 | 0.649 | 0.823 | 0.556 | 0.500 | 0.840 | 0.164 | 0.612 | 0.532 | 0.521 | 0.441 | 0.496 | 0.864 | 0.673 |
| Q_ils vs N_0.95 | 200 | 227 | 0.526 | 0.520 | 0.498 | 0.502 | 0.503 | 0.491 | 0.511 | 0.592 | 0.575 | 0.590 | 0.472 | 0.491 | 0.479 | 0.462 |
| decoy: Q_deep vs decoy | 310 | 200 | 0.475 | 0.463 | 0.474 | 0.487 | 0.520 | 0.498 | 0.499 | 0.518 | 0.520 | 0.508 | 0.519 | 0.514 | 0.522 | 0.529 |
| allele genomes 0 | 573 | 496 | 0.767 | 0.684 | 0.798 | 0.534 | 0.512 | 0.806 | 0.196 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.806 | 0.500 |
| allele genomes 2 | 576 | 496 | 0.803 | 0.660 | 0.807 | 0.561 | 0.506 | 0.821 | 0.184 | 0.783 | 0.643 | 0.615 | 0.390 | 0.476 | 0.858 | 0.833 |
| allele genomes 8 | 575 | 501 | 0.773 | 0.692 | 0.793 | 0.558 | 0.488 | 0.805 | 0.196 | 0.764 | 0.581 | 0.568 | 0.327 | 0.486 | 0.885 | 0.832 |
| congeners 2 | 862 | 746 | 0.797 | 0.702 | 0.795 | 0.553 | 0.480 | 0.837 | 0.166 | 0.626 | 0.560 | 0.524 | 0.417 | 0.470 | 0.818 | 0.646 |
| congeners 6 | 862 | 747 | 0.770 | 0.656 | 0.800 | 0.552 | 0.529 | 0.847 | 0.155 | 0.617 | 0.540 | 0.542 | 0.428 | 0.508 | 0.839 | 0.660 |
| twin 0 | 864 | 758 | 0.873 | 0.759 | 0.869 | 0.663 | 0.495 | 0.822 | 0.182 | 0.628 | 0.541 | 0.505 | 0.368 | 0.485 | 0.837 | 0.676 |
| twin 1 | 860 | 735 | 0.789 | 0.607 | 0.712 | 0.509 | 0.520 | 0.798 | 0.203 | 0.615 | 0.559 | 0.559 | 0.468 | 0.490 | 0.799 | 0.628 |
| regime low | 574 | 495 | 0.770 | 0.683 | 0.794 | 0.546 | 0.503 | 0.796 | 0.206 | 0.634 | 0.553 | 0.559 | 0.409 | 0.492 | 0.780 | 0.653 |
| regime high | 575 | 492 | 0.798 | 0.686 | 0.832 | 0.565 | 0.500 | 0.843 | 0.161 | 0.642 | 0.539 | 0.486 | 0.437 | 0.465 | 0.831 | 0.652 |
| regime none | 575 | 506 | 0.775 | 0.668 | 0.771 | 0.546 | 0.503 | 0.797 | 0.206 | 0.590 | 0.553 | 0.554 | 0.422 | 0.506 | 0.847 | 0.645 |
| depth 3 | 439 | 353 | 0.755 | 0.746 | 0.750 | 0.536 | 0.479 | 0.771 | 0.233 | 0.614 | 0.552 | 0.527 | 0.457 | 0.484 | 0.763 | 0.622 |
| depth 10 | 435 | 382 | 0.783 | 0.709 | 0.808 | 0.564 | 0.524 | 0.837 | 0.166 | 0.607 | 0.542 | 0.529 | 0.440 | 0.493 | 0.836 | 0.654 |
| depth 30 | 424 | 381 | 0.794 | 0.672 | 0.826 | 0.549 | 0.503 | 0.815 | 0.187 | 0.626 | 0.547 | 0.531 | 0.390 | 0.487 | 0.832 | 0.657 |
| depth 100 | 426 | 377 | 0.791 | 0.630 | 0.816 | 0.556 | 0.501 | 0.818 | 0.186 | 0.641 | 0.555 | 0.547 | 0.352 | 0.495 | 0.843 | 0.677 |
| all, within identity bands of 0.005 | 1724 | 1493 | 0.561 | 0.541 | 0.639 | 0.638 | 0.683 | 0.723 | 0.278 | 0.676 | 0.673 | 0.620 | 0.500 | 0.554 | 0.754 | 0.616 |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 158 | 97 | 0.799 | 0.863 | 0.971 | 0.748 | 0.635 | 0.948 | 0.777 |
| Q_deep vs N_0 | congeners 6 | 152 | 98 | 0.800 | 0.897 | 0.993 | 0.659 | 0.614 | 0.994 | 0.743 |
| Q_deep vs N_0 | no twin | 162 | 99 | 0.986 | 0.988 | 0.905 | 0.752 | 0.646 | 0.948 | 0.796 |
| Q_deep vs N_0 | twin | 148 | 96 | 0.930 | 0.794 | 0.905 | 0.661 | 0.601 | 0.903 | 0.699 |
| Q_deep vs N_0 | no allele genomes | 82 | 65 | 0.821 | 0.855 | 0.859 | 0.500 | 0.500 | 0.859 | 0.500 |
| Q_deep vs N_0 | allele genomes | 228 | 130 | 0.788 | 0.896 | 0.922 | 0.830 | 0.683 | 0.964 | 0.920 |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 22 | 0.716 | 0.845 | 1.000 | 0.678 | 0.709 | 1.000 | 1.000 |
| Q_deep vs N_0 | recombination none | 106 | 69 | 0.778 | 0.874 | 0.867 | 0.668 | 0.654 | 0.936 | 0.753 |
| Q_deep vs N_0 | recombination high | 102 | 60 | 0.800 | 0.887 | 0.956 | 0.750 | 0.597 | 0.974 | 0.767 |
| Q_deep vs N_0.5 | congeners 2 | 158 | 102 | 0.708 | 0.801 | 0.855 | 0.732 | 0.643 | 0.829 | 0.736 |
| Q_deep vs N_0.5 | congeners 6 | 152 | 101 | 0.733 | 0.816 | 0.917 | 0.629 | 0.604 | 0.884 | 0.706 |
| Q_deep vs N_0.5 | no twin | 162 | 102 | 0.920 | 0.912 | 0.817 | 0.715 | 0.644 | 0.883 | 0.774 |
| Q_deep vs N_0.5 | twin | 148 | 101 | 0.820 | 0.720 | 0.792 | 0.655 | 0.602 | 0.758 | 0.651 |
| Q_deep vs N_0.5 | no allele genomes | 82 | 69 | 0.741 | 0.804 | 0.751 | 0.500 | 0.500 | 0.751 | 0.500 |
| Q_deep vs N_0.5 | allele genomes | 228 | 134 | 0.710 | 0.816 | 0.833 | 0.796 | 0.679 | 0.868 | 0.847 |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.590 | 0.716 | 0.901 | 0.608 | 0.676 | 0.878 | 0.867 |
| Q_deep vs N_0.5 | recombination none | 106 | 68 | 0.697 | 0.769 | 0.807 | 0.648 | 0.655 | 0.847 | 0.707 |
| Q_deep vs N_0.5 | recombination high | 102 | 73 | 0.711 | 0.829 | 0.819 | 0.732 | 0.602 | 0.869 | 0.738 |
| Q_deep vs N_0.8 | congeners 2 | 158 | 109 | 0.607 | 0.647 | 0.676 | 0.709 | 0.631 | 0.720 | 0.681 |
| Q_deep vs N_0.8 | congeners 6 | 152 | 106 | 0.632 | 0.708 | 0.761 | 0.591 | 0.580 | 0.729 | 0.639 |
| Q_deep vs N_0.8 | no twin | 162 | 112 | 0.712 | 0.725 | 0.661 | 0.676 | 0.628 | 0.746 | 0.702 |
| Q_deep vs N_0.8 | twin | 148 | 103 | 0.670 | 0.633 | 0.708 | 0.630 | 0.580 | 0.676 | 0.620 |
| Q_deep vs N_0.8 | no allele genomes | 82 | 67 | 0.637 | 0.664 | 0.677 | 0.500 | 0.500 | 0.677 | 0.500 |
| Q_deep vs N_0.8 | allele genomes | 228 | 148 | 0.608 | 0.689 | 0.681 | 0.758 | 0.667 | 0.781 | 0.758 |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 24 | 0.596 | 0.650 | 0.715 | 0.593 | 0.663 | 0.762 | 0.747 |
| Q_deep vs N_0.8 | recombination none | 106 | 74 | 0.607 | 0.596 | 0.664 | 0.630 | 0.640 | 0.735 | 0.657 |
| Q_deep vs N_0.8 | recombination high | 102 | 74 | 0.623 | 0.761 | 0.723 | 0.681 | 0.576 | 0.747 | 0.676 |
| Q_deep vs N_0.95 | congeners 2 | 158 | 112 | 0.526 | 0.541 | 0.557 | 0.697 | 0.640 | 0.641 | 0.605 |
| Q_deep vs N_0.95 | congeners 6 | 152 | 115 | 0.539 | 0.534 | 0.593 | 0.601 | 0.620 | 0.607 | 0.570 |
| Q_deep vs N_0.95 | no twin | 162 | 116 | 0.550 | 0.566 | 0.568 | 0.657 | 0.643 | 0.633 | 0.630 |
| Q_deep vs N_0.95 | twin | 148 | 111 | 0.559 | 0.504 | 0.580 | 0.645 | 0.617 | 0.602 | 0.568 |
| Q_deep vs N_0.95 | no allele genomes | 82 | 78 | 0.523 | 0.511 | 0.526 | 0.500 | 0.500 | 0.526 | 0.500 |
| Q_deep vs N_0.95 | allele genomes | 228 | 149 | 0.546 | 0.553 | 0.598 | 0.729 | 0.690 | 0.657 | 0.609 |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 25 | 0.524 | 0.642 | 0.614 | 0.578 | 0.678 | 0.574 | 0.545 |
| Q_deep vs N_0.95 | recombination none | 106 | 76 | 0.512 | 0.529 | 0.568 | 0.641 | 0.659 | 0.627 | 0.599 |
| Q_deep vs N_0.95 | recombination high | 102 | 73 | 0.545 | 0.538 | 0.590 | 0.680 | 0.622 | 0.664 | 0.627 |

## L2c: separation, se (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 1716 | 1491 | 0.761 | 0.669 | 0.770 | 0.563 | 0.506 | 0.782 | 0.222 | 0.618 | 0.553 | 0.536 | 0.437 | 0.489 | 0.798 | - |
| Q_deep vs N_0 | 310 | 201 | 0.767 | 0.753 | 0.850 | 0.563 | 0.569 | 0.874 | 0.141 | 0.685 | 0.636 | 0.618 | 0.487 | 0.664 | 0.905 | - |
| Q_deep vs N_0.5 | 310 | 202 | 0.675 | 0.644 | 0.728 | 0.545 | 0.524 | 0.764 | 0.232 | 0.667 | 0.625 | 0.610 | 0.481 | 0.588 | 0.825 | - |
| Q_deep vs N_0.8 | 310 | 213 | 0.612 | 0.575 | 0.658 | 0.540 | 0.506 | 0.654 | 0.348 | 0.641 | 0.605 | 0.606 | 0.476 | 0.530 | 0.696 | - |
| Q_deep vs N_0.95 | 310 | 225 | 0.525 | 0.507 | 0.525 | 0.523 | 0.511 | 0.563 | 0.436 | 0.630 | 0.615 | 0.612 | 0.487 | 0.530 | 0.609 | - |
| Q_deep vs N_1 | 310 | 226 | 0.470 | 0.472 | 0.478 | 0.519 | 0.510 | 0.507 | 0.492 | 0.635 | 0.622 | 0.619 | 0.493 | 0.541 | 0.564 | - |
| Q_near vs N_0.95 | 333 | 225 | 0.735 | 0.626 | 0.693 | 0.567 | 0.487 | 0.730 | 0.274 | 0.888 | 0.774 | 0.742 | 0.449 | 0.586 | 0.883 | - |
| Q_lone vs N_0.95 | 277 | 225 | 0.758 | 0.644 | 0.739 | 0.559 | 0.485 | 0.728 | 0.262 | 0.330 | 0.305 | 0.306 | 0.382 | 0.293 | 0.553 | - |
| Q_rep vs N_0 | 201 | 201 | 0.984 | 0.881 | 0.983 | 0.637 | 0.552 | 0.997 | 0.017 | 0.461 | 0.368 | 0.320 | 0.445 | 0.403 | 0.991 | - |
| Q_imp vs N (b <= 0.95) | 398 | 841 | 0.797 | 0.714 | 0.811 | 0.571 | 0.510 | 0.835 | 0.172 | 0.656 | 0.583 | 0.560 | 0.407 | 0.492 | 0.843 | - |
| Q (no import) vs N_imp | 1318 | 424 | 0.781 | 0.653 | 0.786 | 0.562 | 0.506 | 0.796 | 0.208 | 0.610 | 0.539 | 0.525 | 0.447 | 0.497 | 0.830 | - |
| Q_ils vs N_0.95 | 197 | 225 | 0.520 | 0.525 | 0.497 | 0.489 | 0.492 | 0.482 | 0.515 | 0.572 | 0.565 | 0.573 | 0.478 | 0.489 | 0.469 | - |
| decoy: Q_deep vs decoy | 310 | 198 | 0.466 | 0.445 | 0.455 | 0.494 | 0.508 | 0.502 | 0.494 | 0.511 | 0.516 | 0.500 | 0.514 | 0.507 | 0.530 | - |
| allele genomes 0 | 570 | 497 | 0.750 | 0.680 | 0.775 | 0.558 | 0.520 | 0.774 | 0.228 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.774 | - |
| allele genomes 2 | 574 | 496 | 0.779 | 0.655 | 0.775 | 0.570 | 0.509 | 0.793 | 0.210 | 0.782 | 0.656 | 0.623 | 0.406 | 0.487 | 0.849 | - |
| allele genomes 8 | 572 | 498 | 0.757 | 0.671 | 0.763 | 0.563 | 0.491 | 0.780 | 0.225 | 0.743 | 0.581 | 0.559 | 0.382 | 0.477 | 0.855 | - |
| congeners 2 | 857 | 747 | 0.771 | 0.684 | 0.766 | 0.570 | 0.480 | 0.803 | 0.201 | 0.626 | 0.567 | 0.530 | 0.438 | 0.476 | 0.805 | - |
| congeners 6 | 859 | 744 | 0.754 | 0.654 | 0.774 | 0.557 | 0.544 | 0.815 | 0.188 | 0.612 | 0.542 | 0.540 | 0.436 | 0.504 | 0.816 | - |
| twin 0 | 860 | 756 | 0.857 | 0.740 | 0.843 | 0.648 | 0.506 | 0.805 | 0.199 | 0.627 | 0.546 | 0.510 | 0.386 | 0.486 | 0.826 | - |
| twin 1 | 856 | 735 | 0.742 | 0.602 | 0.680 | 0.551 | 0.519 | 0.759 | 0.245 | 0.611 | 0.564 | 0.559 | 0.483 | 0.492 | 0.772 | - |
| regime low | 573 | 499 | 0.744 | 0.662 | 0.756 | 0.554 | 0.510 | 0.768 | 0.236 | 0.632 | 0.562 | 0.557 | 0.419 | 0.493 | 0.763 | - |
| regime high | 573 | 489 | 0.780 | 0.670 | 0.802 | 0.586 | 0.503 | 0.819 | 0.187 | 0.639 | 0.546 | 0.500 | 0.452 | 0.467 | 0.813 | - |
| regime none | 570 | 503 | 0.761 | 0.674 | 0.756 | 0.552 | 0.506 | 0.766 | 0.236 | 0.587 | 0.550 | 0.550 | 0.440 | 0.504 | 0.822 | - |
| depth 3 | 431 | 351 | 0.724 | 0.718 | 0.710 | 0.543 | 0.485 | 0.719 | 0.288 | 0.600 | 0.565 | 0.534 | 0.481 | 0.493 | 0.728 | - |
| depth 10 | 435 | 382 | 0.769 | 0.717 | 0.780 | 0.571 | 0.529 | 0.818 | 0.185 | 0.605 | 0.542 | 0.526 | 0.447 | 0.491 | 0.824 | - |
| depth 30 | 424 | 381 | 0.777 | 0.678 | 0.799 | 0.563 | 0.505 | 0.797 | 0.206 | 0.638 | 0.554 | 0.540 | 0.414 | 0.486 | 0.819 | - |
| depth 100 | 426 | 377 | 0.779 | 0.626 | 0.802 | 0.566 | 0.508 | 0.802 | 0.201 | 0.644 | 0.558 | 0.548 | 0.375 | 0.493 | 0.841 | - |
| all, within identity bands of 0.005 | 1716 | 1491 | 0.557 | 0.545 | 0.631 | 0.649 | 0.656 | 0.696 | 0.308 | 0.677 | 0.674 | 0.629 | 0.505 | 0.561 | 0.735 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 158 | 99 | 0.750 | 0.824 | 0.919 | 0.734 | 0.639 | 0.922 | - |
| Q_deep vs N_0 | congeners 6 | 152 | 102 | 0.783 | 0.872 | 0.978 | 0.638 | 0.597 | 0.981 | - |
| Q_deep vs N_0 | no twin | 162 | 101 | 0.973 | 0.948 | 0.890 | 0.737 | 0.649 | 0.945 | - |
| Q_deep vs N_0 | twin | 148 | 100 | 0.846 | 0.777 | 0.860 | 0.640 | 0.587 | 0.871 | - |
| Q_deep vs N_0 | no allele genomes | 82 | 68 | 0.778 | 0.842 | 0.846 | 0.500 | 0.500 | 0.846 | - |
| Q_deep vs N_0 | allele genomes | 228 | 133 | 0.767 | 0.868 | 0.893 | 0.797 | 0.671 | 0.937 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 37 | 22 | 0.721 | 0.894 | 1.000 | 0.655 | 0.619 | 1.000 | - |
| Q_deep vs N_0 | recombination none | 105 | 70 | 0.764 | 0.856 | 0.849 | 0.665 | 0.640 | 0.933 | - |
| Q_deep vs N_0 | recombination high | 102 | 62 | 0.753 | 0.828 | 0.910 | 0.718 | 0.616 | 0.947 | - |
| Q_deep vs N_0.5 | congeners 2 | 158 | 100 | 0.648 | 0.699 | 0.757 | 0.732 | 0.644 | 0.836 | - |
| Q_deep vs N_0.5 | congeners 6 | 152 | 102 | 0.705 | 0.756 | 0.874 | 0.602 | 0.576 | 0.877 | - |
| Q_deep vs N_0.5 | no twin | 162 | 101 | 0.874 | 0.840 | 0.788 | 0.702 | 0.633 | 0.879 | - |
| Q_deep vs N_0.5 | twin | 148 | 101 | 0.701 | 0.638 | 0.737 | 0.638 | 0.586 | 0.772 | - |
| Q_deep vs N_0.5 | no allele genomes | 82 | 68 | 0.669 | 0.758 | 0.708 | 0.500 | 0.500 | 0.708 | - |
| Q_deep vs N_0.5 | allele genomes | 228 | 134 | 0.686 | 0.726 | 0.799 | 0.774 | 0.658 | 0.878 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.616 | 0.649 | 0.874 | 0.584 | 0.644 | 0.939 | - |
| Q_deep vs N_0.5 | recombination none | 105 | 66 | 0.676 | 0.711 | 0.768 | 0.631 | 0.631 | 0.855 | - |
| Q_deep vs N_0.5 | recombination high | 102 | 72 | 0.669 | 0.751 | 0.768 | 0.708 | 0.595 | 0.860 | - |
| Q_deep vs N_0.8 | congeners 2 | 158 | 111 | 0.589 | 0.631 | 0.609 | 0.710 | 0.645 | 0.713 | - |
| Q_deep vs N_0.8 | congeners 6 | 152 | 102 | 0.637 | 0.686 | 0.755 | 0.568 | 0.562 | 0.709 | - |
| Q_deep vs N_0.8 | no twin | 162 | 112 | 0.718 | 0.707 | 0.659 | 0.667 | 0.613 | 0.727 | - |
| Q_deep vs N_0.8 | twin | 148 | 101 | 0.614 | 0.606 | 0.650 | 0.619 | 0.594 | 0.666 | - |
| Q_deep vs N_0.8 | no allele genomes | 82 | 66 | 0.628 | 0.694 | 0.659 | 0.500 | 0.500 | 0.659 | - |
| Q_deep vs N_0.8 | allele genomes | 228 | 147 | 0.611 | 0.648 | 0.654 | 0.738 | 0.668 | 0.743 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 37 | 23 | 0.627 | 0.645 | 0.731 | 0.564 | 0.574 | 0.695 | - |
| Q_deep vs N_0.8 | recombination none | 105 | 72 | 0.596 | 0.604 | 0.649 | 0.632 | 0.631 | 0.718 | - |
| Q_deep vs N_0.8 | recombination high | 102 | 74 | 0.633 | 0.720 | 0.691 | 0.658 | 0.579 | 0.714 | - |
| Q_deep vs N_0.95 | congeners 2 | 158 | 111 | 0.529 | 0.512 | 0.542 | 0.685 | 0.627 | 0.631 | - |
| Q_deep vs N_0.95 | congeners 6 | 152 | 114 | 0.525 | 0.543 | 0.581 | 0.578 | 0.590 | 0.593 | - |
| Q_deep vs N_0.95 | no twin | 162 | 115 | 0.551 | 0.557 | 0.562 | 0.655 | 0.633 | 0.644 | - |
| Q_deep vs N_0.95 | twin | 148 | 110 | 0.533 | 0.497 | 0.563 | 0.608 | 0.590 | 0.576 | - |
| Q_deep vs N_0.95 | no allele genomes | 82 | 76 | 0.548 | 0.568 | 0.537 | 0.500 | 0.500 | 0.537 | - |
| Q_deep vs N_0.95 | allele genomes | 228 | 149 | 0.524 | 0.513 | 0.579 | 0.706 | 0.664 | 0.640 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 37 | 24 | 0.516 | 0.626 | 0.569 | 0.545 | 0.611 | 0.546 | - |
| Q_deep vs N_0.95 | recombination none | 105 | 76 | 0.507 | 0.524 | 0.547 | 0.620 | 0.637 | 0.625 | - |
| Q_deep vs N_0.95 | recombination high | 102 | 71 | 0.529 | 0.516 | 0.592 | 0.665 | 0.605 | 0.636 | - |

## L2c: separation, hifi (target rows: Q = its strain present, N = a novel species' reads)

AUC of each feature, Q against N (above 0.5: higher in Q). fixed_agreement_approx = ancestry_agreement + ancestry_fixed_gain; oracle_fixed_agreement = the agreement at the sites no stored allele varies at (the oracle's, unweighted; -1 without alleles): the ceiling of the polymorphic group's idea.
| comparison | Q | N | identity | top_identity | lu_per_kb | uniqueness | ancestry_sites_per_record | ancestry_agreement | ancestry_congener_share | allele_explained_share | allele_identity_gain | polymorphic_known_share | polymorphic_novel_share | ancestry_fixed_gain | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all Q vs all N | 576 | 513 | 0.814 | 0.761 | 0.877 | 0.598 | 0.501 | 0.836 | 0.169 | 0.660 | 0.596 | 0.538 | 0.296 | 0.495 | 0.827 | - |
| Q_deep vs N_0 | 101 | 75 | 0.804 | 0.793 | 0.924 | 0.677 | 0.488 | 0.935 | 0.076 | 0.744 | 0.711 | 0.707 | 0.460 | 0.725 | 0.940 | - |
| Q_deep vs N_0.5 | 101 | 74 | 0.777 | 0.767 | 0.967 | 0.604 | 0.456 | 0.845 | 0.168 | 0.681 | 0.668 | 0.662 | 0.358 | 0.609 | 0.823 | - |
| Q_deep vs N_0.8 | 101 | 65 | 0.657 | 0.717 | 0.791 | 0.600 | 0.489 | 0.698 | 0.314 | 0.672 | 0.673 | 0.684 | 0.385 | 0.573 | 0.717 | - |
| Q_deep vs N_0.95 | 101 | 73 | 0.584 | 0.638 | 0.628 | 0.557 | 0.455 | 0.561 | 0.442 | 0.640 | 0.651 | 0.664 | 0.404 | 0.525 | 0.614 | - |
| Q_deep vs N_1 | 101 | 72 | 0.538 | 0.617 | 0.582 | 0.518 | 0.477 | 0.522 | 0.480 | 0.622 | 0.644 | 0.652 | 0.358 | 0.469 | 0.520 | - |
| Q_near vs N_0.95 | 100 | 73 | 0.813 | 0.836 | 0.806 | 0.615 | 0.479 | 0.784 | 0.219 | 0.999 | 0.872 | 0.718 | 0.277 | 0.592 | 0.985 | - |
| Q_lone vs N_0.95 | 92 | 73 | 0.748 | 0.768 | 0.770 | 0.540 | 0.456 | 0.757 | 0.243 | 0.338 | 0.307 | 0.308 | 0.217 | 0.274 | 0.599 | - |
| Q_rep vs N_0 | 73 | 75 | 0.991 | 0.960 | 0.999 | 0.727 | 0.513 | 1.000 | 0.000 | 0.640 | 0.495 | 0.380 | 0.387 | 0.520 | 1.000 | - |
| Q_imp vs N (b <= 0.95) | 136 | 287 | 0.825 | 0.825 | 0.907 | 0.617 | 0.548 | 0.875 | 0.130 | 0.659 | 0.612 | 0.529 | 0.292 | 0.482 | 0.832 | - |
| Q (no import) vs N_imp | 440 | 154 | 0.857 | 0.662 | 0.932 | 0.602 | 0.482 | 0.872 | 0.137 | 0.664 | 0.584 | 0.538 | 0.281 | 0.511 | 0.892 | - |
| Q_ils vs N_0.95 | 74 | 73 | 0.546 | 0.578 | 0.537 | 0.474 | 0.489 | 0.503 | 0.504 | 0.627 | 0.628 | 0.628 | 0.387 | 0.512 | 0.490 | - |
| decoy: Q_deep vs decoy | 101 | 63 | 0.476 | 0.471 | 0.487 | 0.529 | 0.489 | 0.450 | 0.551 | 0.529 | 0.532 | 0.537 | 0.529 | 0.513 | 0.490 | - |
| allele genomes 0 | 192 | 172 | 0.806 | 0.740 | 0.862 | 0.585 | 0.496 | 0.827 | 0.178 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.827 | - |
| allele genomes 2 | 192 | 171 | 0.825 | 0.763 | 0.889 | 0.585 | 0.507 | 0.844 | 0.160 | 0.864 | 0.768 | 0.655 | 0.140 | 0.507 | 0.873 | - |
| allele genomes 8 | 192 | 170 | 0.818 | 0.780 | 0.883 | 0.634 | 0.501 | 0.845 | 0.162 | 0.872 | 0.665 | 0.524 | 0.054 | 0.473 | 0.902 | - |
| congeners 2 | 288 | 257 | 0.834 | 0.792 | 0.900 | 0.605 | 0.501 | 0.898 | 0.109 | 0.670 | 0.614 | 0.526 | 0.317 | 0.479 | 0.827 | - |
| congeners 6 | 288 | 256 | 0.796 | 0.737 | 0.859 | 0.595 | 0.509 | 0.863 | 0.140 | 0.652 | 0.582 | 0.552 | 0.275 | 0.515 | 0.850 | - |
| twin 0 | 288 | 259 | 0.917 | 0.855 | 0.921 | 0.751 | 0.496 | 0.855 | 0.150 | 0.668 | 0.598 | 0.508 | 0.273 | 0.480 | 0.855 | - |
| twin 1 | 288 | 254 | 0.871 | 0.743 | 0.828 | 0.525 | 0.518 | 0.822 | 0.181 | 0.653 | 0.597 | 0.567 | 0.303 | 0.506 | 0.802 | - |
| regime low | 192 | 168 | 0.787 | 0.718 | 0.866 | 0.575 | 0.505 | 0.811 | 0.194 | 0.673 | 0.617 | 0.592 | 0.305 | 0.518 | 0.775 | - |
| regime high | 192 | 166 | 0.854 | 0.843 | 0.951 | 0.630 | 0.484 | 0.899 | 0.111 | 0.702 | 0.609 | 0.492 | 0.300 | 0.468 | 0.874 | - |
| regime none | 192 | 179 | 0.806 | 0.729 | 0.835 | 0.588 | 0.511 | 0.808 | 0.197 | 0.614 | 0.550 | 0.538 | 0.280 | 0.511 | 0.852 | - |
| depth 1 | 159 | 135 | 0.807 | 0.754 | 0.885 | 0.610 | 0.523 | 0.848 | 0.158 | 0.672 | 0.596 | 0.545 | 0.297 | 0.497 | 0.828 | - |
| depth 2 | 136 | 122 | 0.815 | 0.751 | 0.876 | 0.593 | 0.493 | 0.824 | 0.183 | 0.622 | 0.551 | 0.496 | 0.250 | 0.457 | 0.807 | - |
| depth 4 | 138 | 127 | 0.829 | 0.795 | 0.889 | 0.571 | 0.481 | 0.833 | 0.170 | 0.669 | 0.626 | 0.561 | 0.311 | 0.515 | 0.837 | - |
| depth 0.5 | 143 | 129 | 0.809 | 0.747 | 0.862 | 0.601 | 0.506 | 0.842 | 0.164 | 0.672 | 0.605 | 0.547 | 0.321 | 0.510 | 0.839 | - |
| all, within identity bands of 0.005 | 576 | 513 | 0.617 | 0.578 | 0.724 | 0.645 | 0.688 | 0.714 | 0.288 | 0.634 | 0.638 | 0.580 | 0.331 | 0.495 | 0.696 | - |

The deep strains against the novel species by where they left S's stem, within the strata that decide the comparison (the consensus needs 3 congeners; the fixed sites need allele genomes).
| pair | stratum | Q | N | identity | lu_per_kb | ancestry_agreement | allele_explained_share | polymorphic_known_share | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|
| Q_deep vs N_0 | congeners 2 | 44 | 34 | 0.779 | 0.900 | 1.000 | 0.800 | 0.750 | 0.988 | - |
| Q_deep vs N_0 | congeners 6 | 57 | 41 | 0.817 | 0.940 | 1.000 | 0.707 | 0.683 | 0.998 | - |
| Q_deep vs N_0 | no twin | 48 | 38 | 0.997 | 0.993 | 0.981 | 0.797 | 0.734 | 0.987 | - |
| Q_deep vs N_0 | twin | 53 | 37 | 0.916 | 0.829 | 0.897 | 0.697 | 0.668 | 0.896 | - |
| Q_deep vs N_0 | no allele genomes | 27 | 31 | 0.827 | 0.892 | 0.912 | 0.500 | 0.500 | 0.912 | - |
| Q_deep vs N_0 | allele genomes | 74 | 44 | 0.806 | 0.939 | 0.943 | 0.886 | 0.811 | 0.970 | - |
| Q_deep vs N_0 | 6 congeners, allele genomes, no recombination | 11 | 7 | 0.870 | 0.987 | 1.000 | 0.688 | 0.636 | 1.000 | - |
| Q_deep vs N_0 | recombination none | 35 | 25 | 0.835 | 0.917 | 0.914 | 0.703 | 0.680 | 0.952 | - |
| Q_deep vs N_0 | recombination high | 29 | 23 | 0.835 | 0.970 | 1.000 | 0.785 | 0.699 | 0.991 | - |
| Q_deep vs N_0.5 | congeners 2 | 44 | 41 | 0.761 | 0.975 | 0.978 | 0.779 | 0.707 | 0.894 | - |
| Q_deep vs N_0.5 | congeners 6 | 57 | 33 | 0.791 | 0.960 | 0.953 | 0.596 | 0.618 | 0.895 | - |
| Q_deep vs N_0.5 | no twin | 48 | 38 | 0.986 | 0.997 | 0.859 | 0.729 | 0.668 | 0.897 | - |
| Q_deep vs N_0.5 | twin | 53 | 36 | 0.975 | 0.941 | 0.822 | 0.645 | 0.665 | 0.744 | - |
| Q_deep vs N_0.5 | no allele genomes | 27 | 24 | 0.843 | 0.946 | 0.733 | 0.500 | 0.500 | 0.733 | - |
| Q_deep vs N_0.5 | allele genomes | 74 | 50 | 0.759 | 0.974 | 0.884 | 0.808 | 0.769 | 0.894 | - |
| Q_deep vs N_0.5 | 6 congeners, allele genomes, no recombination | 11 | 6 | 0.773 | 0.773 | 0.924 | 0.727 | 0.742 | 0.985 | - |
| Q_deep vs N_0.5 | recombination none | 35 | 22 | 0.816 | 0.945 | 0.818 | 0.714 | 0.717 | 0.903 | - |
| Q_deep vs N_0.5 | recombination high | 29 | 25 | 0.761 | 1.000 | 0.960 | 0.770 | 0.670 | 0.928 | - |
| Q_deep vs N_0.8 | congeners 2 | 44 | 33 | 0.663 | 0.815 | 0.796 | 0.742 | 0.679 | 0.781 | - |
| Q_deep vs N_0.8 | congeners 6 | 57 | 32 | 0.640 | 0.770 | 0.758 | 0.660 | 0.691 | 0.731 | - |
| Q_deep vs N_0.8 | no twin | 48 | 30 | 0.797 | 0.821 | 0.674 | 0.713 | 0.695 | 0.771 | - |
| Q_deep vs N_0.8 | twin | 53 | 35 | 0.805 | 0.767 | 0.736 | 0.647 | 0.677 | 0.678 | - |
| Q_deep vs N_0.8 | no allele genomes | 27 | 23 | 0.713 | 0.760 | 0.671 | 0.500 | 0.500 | 0.671 | - |
| Q_deep vs N_0.8 | allele genomes | 74 | 42 | 0.645 | 0.823 | 0.703 | 0.772 | 0.796 | 0.786 | - |
| Q_deep vs N_0.8 | 6 congeners, allele genomes, no recombination | 11 | 9 | 0.596 | 0.758 | 0.737 | 0.697 | 0.737 | 0.859 | - |
| Q_deep vs N_0.8 | recombination none | 35 | 23 | 0.629 | 0.684 | 0.717 | 0.634 | 0.647 | 0.742 | - |
| Q_deep vs N_0.8 | recombination high | 29 | 17 | 0.623 | 0.978 | 0.848 | 0.757 | 0.728 | 0.805 | - |
| Q_deep vs N_0.95 | congeners 2 | 44 | 38 | 0.566 | 0.678 | 0.623 | 0.702 | 0.664 | 0.696 | - |
| Q_deep vs N_0.95 | congeners 6 | 57 | 35 | 0.596 | 0.604 | 0.591 | 0.615 | 0.665 | 0.595 | - |
| Q_deep vs N_0.95 | no twin | 48 | 37 | 0.652 | 0.676 | 0.610 | 0.675 | 0.684 | 0.660 | - |
| Q_deep vs N_0.95 | twin | 53 | 36 | 0.593 | 0.593 | 0.528 | 0.632 | 0.650 | 0.584 | - |
| Q_deep vs N_0.95 | no allele genomes | 27 | 23 | 0.652 | 0.651 | 0.583 | 0.500 | 0.500 | 0.583 | - |
| Q_deep vs N_0.95 | allele genomes | 74 | 50 | 0.575 | 0.621 | 0.547 | 0.731 | 0.780 | 0.669 | - |
| Q_deep vs N_0.95 | 6 congeners, allele genomes, no recombination | 11 | 7 | 0.623 | 0.701 | 0.636 | 0.727 | 0.727 | 0.831 | - |
| Q_deep vs N_0.95 | recombination none | 35 | 25 | 0.609 | 0.619 | 0.579 | 0.638 | 0.648 | 0.675 | - |
| Q_deep vs N_0.95 | recombination high | 29 | 26 | 0.576 | 0.767 | 0.655 | 0.717 | 0.725 | 0.702 | - |

## L2d: the gain of each group over identity and depth (logistic models, genera held out)

AUC of grouped-CV logistic models on the target rows; gain over identity, top identity and fragments with a 95% interval over resampled genera.
| set | model | auc_base | auc | gain | gain_lo | gain_hi |
|---|---|---|---|---|---|---|
| pe | identity, depth +ancestry | 0.780 | 0.820 | 0.040 | 0.020 | 0.060 |
| pe | identity, depth +alleles | 0.780 | 0.855 | 0.075 | 0.035 | 0.123 |
| pe | identity, depth +polymorphic | 0.780 | 0.798 | 0.018 | -0.014 | 0.046 |
| pe | identity, depth +all three | 0.780 | 0.886 | 0.106 | 0.058 | 0.156 |
| pe | identity, depth +ancestry, polymorphic | 0.780 | 0.852 | 0.072 | 0.046 | 0.101 |
| pe | identity, depth +oracle fixed agreement | 0.780 | 0.780 | 0.000 | -0.006 | 0.006 |
| pe, genera with alleles | identity, depth +ancestry | 0.784 | 0.824 | 0.040 | 0.013 | 0.066 |
| pe, genera with alleles | identity, depth +alleles | 0.784 | 0.927 | 0.142 | 0.082 | 0.202 |
| pe, genera with alleles | identity, depth +polymorphic | 0.784 | 0.815 | 0.031 | -0.008 | 0.075 |
| pe, genera with alleles | identity, depth +all three | 0.784 | 0.946 | 0.162 | 0.104 | 0.214 |
| pe, genera with alleles | identity, depth +ancestry, polymorphic | 0.784 | 0.908 | 0.124 | 0.085 | 0.168 |
| pe, genera with alleles | identity, depth +oracle fixed agreement | 0.784 | 0.854 | 0.070 | 0.043 | 0.092 |
| se | identity, depth +ancestry | 0.763 | 0.796 | 0.034 | 0.016 | 0.053 |
| se | identity, depth +alleles | 0.763 | 0.846 | 0.084 | 0.040 | 0.131 |
| se | identity, depth +polymorphic | 0.763 | 0.784 | 0.021 | 0.000 | 0.045 |
| se | identity, depth +all three | 0.763 | 0.871 | 0.109 | 0.062 | 0.161 |
| se | identity, depth +ancestry, polymorphic | 0.763 | 0.832 | 0.070 | 0.045 | 0.096 |
| se | identity, depth +oracle fixed agreement | 0.763 | 0.763 | 0.000 | 0.000 | 0.000 |
| se, genera with alleles | identity, depth +ancestry | 0.765 | 0.800 | 0.034 | 0.008 | 0.056 |
| se, genera with alleles | identity, depth +alleles | 0.765 | 0.914 | 0.149 | 0.084 | 0.209 |
| se, genera with alleles | identity, depth +polymorphic | 0.765 | 0.802 | 0.037 | 0.008 | 0.075 |
| se, genera with alleles | identity, depth +all three | 0.765 | 0.934 | 0.169 | 0.104 | 0.227 |
| se, genera with alleles | identity, depth +ancestry, polymorphic | 0.765 | 0.884 | 0.119 | 0.080 | 0.156 |
| se, genera with alleles | identity, depth +oracle fixed agreement | 0.765 | 0.765 | 0.000 | 0.000 | 0.000 |
| hifi | identity, depth +ancestry | 0.805 | 0.849 | 0.044 | 0.013 | 0.071 |
| hifi | identity, depth +alleles | 0.805 | 0.886 | 0.081 | 0.039 | 0.129 |
| hifi | identity, depth +polymorphic | 0.805 | 0.885 | 0.080 | 0.032 | 0.120 |
| hifi | identity, depth +all three | 0.805 | 0.937 | 0.131 | 0.081 | 0.174 |
| hifi | identity, depth +ancestry, polymorphic | 0.805 | 0.911 | 0.106 | 0.064 | 0.139 |
| hifi | identity, depth +oracle fixed agreement | 0.805 | 0.805 | 0.000 | 0.000 | 0.000 |
| hifi, genera with alleles | identity, depth +ancestry | 0.806 | 0.847 | 0.042 | 0.010 | 0.072 |
| hifi, genera with alleles | identity, depth +alleles | 0.806 | 0.952 | 0.147 | 0.082 | 0.200 |
| hifi, genera with alleles | identity, depth +polymorphic | 0.806 | 0.924 | 0.118 | 0.057 | 0.166 |
| hifi, genera with alleles | identity, depth +all three | 0.806 | 0.971 | 0.166 | 0.106 | 0.215 |
| hifi, genera with alleles | identity, depth +ancestry, polymorphic | 0.806 | 0.942 | 0.136 | 0.081 | 0.178 |
| hifi, genera with alleles | identity, depth +oracle fixed agreement | 0.806 | 0.806 | 0.000 | 0.000 | 0.000 |

## Sensitivity: by divergence of the lineage and by recombination with congeners

Per feature, the threshold 5% of all N rows exceed (the feature oriented by its AUC); rate_above: the share of Q rows above it (sensitivity) or of N rows (false-positive rate) in each bin. 'MRCA with nearest allele genome' is how recently the strain shared an ancestor with a strain the database knows; 'genes with a congener's segment' how much of it came from a congener.
| set | class | by | level | rows | identity | ancestry_agreement | allele_explained_share | polymorphic_known_share | fixed_agreement_approx | oracle_fixed_agreement |
|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 278 | 0.299 | 0.392 | 0.403 | 0.281 | 0.410 | 0.000 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.062 | 0.109 | 0.109 | 0.031 | 0.078 | 0.000 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 573 | 0.351 | 0.525 | 0.000 | 0.000 | 0.079 | 0.000 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 809 | 0.396 | 0.528 | 0.727 | 0.447 | 0.630 | 0.000 |
| pe | Q | depth | 10 | 435 | 0.359 | 0.506 | 0.405 | 0.241 | 0.377 | 0.000 |
| pe | Q | depth | 100 | 426 | 0.329 | 0.467 | 0.425 | 0.265 | 0.354 | 0.000 |
| pe | Q | depth | 3 | 439 | 0.405 | 0.531 | 0.408 | 0.251 | 0.483 | 0.000 |
| pe | Q | depth | 30 | 424 | 0.316 | 0.453 | 0.403 | 0.269 | 0.347 | 0.000 |
| pe | Q | distance to rep | 0.005-0.015 | 790 | 0.371 | 0.515 | 0.470 | 0.261 | 0.370 | 0.000 |
| pe | Q | distance to rep | 0.015-0.03 | 433 | 0.053 | 0.337 | 0.508 | 0.333 | 0.333 | 0.000 |
| pe | Q | distance to rep | <0.005 | 298 | 0.977 | 0.966 | 0.070 | 0.067 | 0.685 | 0.000 |
| pe | Q | distance to rep | >=0.03 | 203 | 0.005 | 0.015 | 0.468 | 0.355 | 0.167 | 0.000 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 570 | 0.319 | 0.486 | 0.428 | 0.258 | 0.305 | 0.000 |
| pe | Q | genes with a congener's segment | <=0.1 | 454 | 0.333 | 0.467 | 0.401 | 0.200 | 0.436 | 0.000 |
| pe | Q | genes with a congener's segment | >0.3 | 251 | 0.315 | 0.538 | 0.486 | 0.291 | 0.295 | 0.000 |
| pe | Q | genes with a congener's segment | none | 449 | 0.437 | 0.490 | 0.354 | 0.292 | 0.508 | 0.000 |
| pe | Q | regime | high | 575 | 0.358 | 0.520 | 0.442 | 0.233 | 0.290 | 0.000 |
| pe | Q | regime | low | 574 | 0.293 | 0.455 | 0.416 | 0.232 | 0.408 | 0.000 |
| pe | Q | regime | none | 575 | 0.407 | 0.494 | 0.372 | 0.304 | 0.475 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 933 | 0.345 | 0.475 | 0.728 | 0.442 | 0.613 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 589 | 0.367 | 0.516 | 0.000 | 0.002 | 0.076 | 0.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 198 | 0.333 | 0.470 | 0.141 | 0.146 | 0.268 | 0.000 |
| pe | Q | role | Q_deep | 310 | 0.129 | 0.165 | 0.432 | 0.300 | 0.281 | 0.000 |
| pe | Q | role | Q_ils | 200 | 0.090 | 0.100 | 0.375 | 0.280 | 0.050 | 0.000 |
| pe | Q | role | Q_imp | 402 | 0.316 | 0.535 | 0.493 | 0.308 | 0.341 | 0.000 |
| pe | Q | role | Q_lone | 277 | 0.397 | 0.617 | 0.105 | 0.087 | 0.191 | 0.000 |
| pe | Q | role | Q_near | 334 | 0.407 | 0.563 | 0.751 | 0.434 | 0.754 | 0.000 |
| pe | Q | role | Q_rep | 201 | 0.881 | 0.990 | 0.100 | 0.000 | 0.672 | 0.000 |
| pe | Q | twin genus | 0 | 864 | 0.248 | 0.543 | 0.397 | 0.211 | 0.407 | 0.000 |
| pe | Q | twin genus | 1 | 860 | 0.458 | 0.436 | 0.423 | 0.302 | 0.374 | 0.000 |
| pe | N | depth | 10 | 382 | 0.052 | 0.029 | 0.052 | 0.073 | 0.039 | 0.000 |
| pe | N | depth | 100 | 377 | 0.011 | 0.027 | 0.011 | 0.008 | 0.011 | 0.000 |
| pe | N | depth | 3 | 353 | 0.125 | 0.125 | 0.102 | 0.079 | 0.133 | 0.000 |
| pe | N | depth | 30 | 381 | 0.018 | 0.026 | 0.039 | 0.026 | 0.024 | 0.000 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 860 | 0.043 | 0.056 | 0.043 | 0.048 | 0.048 | 0.000 |
| pe | N | genes with a segment of S or a congener | >0.3 | 266 | 0.045 | 0.015 | 0.049 | 0.060 | 0.026 | 0.000 |
| pe | N | genes with a segment of S or a congener | none | 367 | 0.071 | 0.063 | 0.068 | 0.033 | 0.074 | 0.000 |
| pe | N | regime | high | 492 | 0.063 | 0.057 | 0.053 | 0.057 | 0.043 | 0.000 |
| pe | N | regime | low | 495 | 0.024 | 0.044 | 0.040 | 0.053 | 0.051 | 0.000 |
| pe | N | regime | none | 506 | 0.063 | 0.049 | 0.057 | 0.030 | 0.057 | 0.000 |
| pe | N | role | N_0 | 195 | 0.010 | 0.000 | 0.031 | 0.056 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 203 | 0.025 | 0.015 | 0.025 | 0.039 | 0.030 | 0.000 |
| pe | N | role | N_0.8 | 215 | 0.042 | 0.028 | 0.051 | 0.042 | 0.037 | 0.000 |
| pe | N | role | N_0.95 | 227 | 0.079 | 0.110 | 0.079 | 0.035 | 0.079 | 0.000 |
| pe | N | role | N_1 | 228 | 0.110 | 0.154 | 0.092 | 0.044 | 0.154 | 0.000 |
| pe | N | role | N_imp0.2 | 213 | 0.033 | 0.014 | 0.028 | 0.052 | 0.014 | 0.000 |
| pe | N | role | N_imp0.4 | 212 | 0.042 | 0.014 | 0.038 | 0.057 | 0.024 | 0.000 |
| pe | N | twin genus | 0 | 758 | 0.003 | 0.063 | 0.026 | 0.021 | 0.037 | 0.000 |
| pe | N | twin genus | 1 | 735 | 0.099 | 0.037 | 0.075 | 0.072 | 0.064 | 0.000 |
| se | Q | MRCA with nearest allele genome | 0.002-0.01 | 276 | 0.283 | 0.362 | 0.399 | 0.228 | 0.080 |  |
| se | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.062 | 0.094 | 0.109 | 0.031 | 0.000 |  |
| se | Q | MRCA with nearest allele genome | no allele genome | 570 | 0.326 | 0.456 | 0.000 | 0.000 | 0.000 |  |
| se | Q | MRCA with nearest allele genome | recent (<0.002) | 806 | 0.362 | 0.483 | 0.710 | 0.402 | 0.127 |  |
| se | Q | depth | 10 | 435 | 0.356 | 0.451 | 0.386 | 0.218 | 0.103 |  |
| se | Q | depth | 100 | 426 | 0.291 | 0.408 | 0.427 | 0.237 | 0.038 |  |
| se | Q | depth | 3 | 431 | 0.371 | 0.485 | 0.367 | 0.211 | 0.053 |  |
| se | Q | depth | 30 | 424 | 0.285 | 0.415 | 0.427 | 0.241 | 0.094 |  |
| se | Q | distance to rep | 0.005-0.015 | 786 | 0.322 | 0.463 | 0.457 | 0.229 | 0.078 |  |
| se | Q | distance to rep | 0.015-0.03 | 430 | 0.060 | 0.258 | 0.495 | 0.307 | 0.044 |  |
| se | Q | distance to rep | <0.005 | 298 | 0.940 | 0.930 | 0.084 | 0.064 | 0.134 |  |
| se | Q | distance to rep | >=0.03 | 202 | 0.005 | 0.015 | 0.455 | 0.287 | 0.020 |  |
| se | Q | genes with a congener's segment | 0.1-0.3 | 569 | 0.322 | 0.432 | 0.417 | 0.213 | 0.065 |  |
| se | Q | genes with a congener's segment | <=0.1 | 453 | 0.278 | 0.406 | 0.395 | 0.177 | 0.071 |  |
| se | Q | genes with a congener's segment | >0.3 | 248 | 0.282 | 0.488 | 0.472 | 0.278 | 0.048 |  |
| se | Q | genes with a congener's segment | none | 446 | 0.406 | 0.457 | 0.350 | 0.267 | 0.096 |  |
| se | Q | regime | high | 573 | 0.333 | 0.473 | 0.428 | 0.188 | 0.061 |  |
| se | Q | regime | low | 573 | 0.267 | 0.396 | 0.405 | 0.204 | 0.072 |  |
| se | Q | regime | none | 570 | 0.379 | 0.451 | 0.372 | 0.288 | 0.084 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 929 | 0.319 | 0.431 | 0.715 | 0.392 | 0.126 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 586 | 0.334 | 0.451 | 0.000 | 0.000 | 0.000 |  |
| se | Q | rep-lineage sites it lacks that the alleles cover | partly | 197 | 0.325 | 0.442 | 0.127 | 0.127 | 0.036 |  |
| se | Q | role | Q_deep | 310 | 0.087 | 0.145 | 0.419 | 0.245 | 0.058 |  |
| se | Q | role | Q_ils | 197 | 0.071 | 0.076 | 0.340 | 0.254 | 0.015 |  |
| se | Q | role | Q_imp | 398 | 0.307 | 0.477 | 0.477 | 0.291 | 0.063 |  |
| se | Q | role | Q_lone | 277 | 0.365 | 0.498 | 0.112 | 0.079 | 0.043 |  |
| se | Q | role | Q_near | 333 | 0.369 | 0.526 | 0.754 | 0.372 | 0.153 |  |
| se | Q | role | Q_rep | 201 | 0.861 | 0.955 | 0.100 | 0.005 | 0.075 |  |
| se | Q | twin genus | 0 | 860 | 0.237 | 0.465 | 0.393 | 0.176 | 0.052 |  |
| se | Q | twin genus | 1 | 856 | 0.416 | 0.415 | 0.410 | 0.278 | 0.092 |  |
| se | N | depth | 10 | 382 | 0.047 | 0.029 | 0.071 | 0.081 | 0.016 |  |
| se | N | depth | 100 | 377 | 0.005 | 0.013 | 0.005 | 0.008 | 0.000 |  |
| se | N | depth | 3 | 351 | 0.134 | 0.145 | 0.094 | 0.100 | 0.028 |  |
| se | N | depth | 30 | 381 | 0.021 | 0.010 | 0.034 | 0.013 | 0.005 |  |
| se | N | genes with a segment of S or a congener | <=0.3 | 865 | 0.052 | 0.053 | 0.039 | 0.049 | 0.009 |  |
| se | N | genes with a segment of S or a congener | >0.3 | 264 | 0.038 | 0.027 | 0.049 | 0.053 | 0.004 |  |
| se | N | genes with a segment of S or a congener | none | 362 | 0.055 | 0.050 | 0.077 | 0.050 | 0.025 |  |
| se | N | regime | high | 489 | 0.065 | 0.051 | 0.049 | 0.053 | 0.008 |  |
| se | N | regime | low | 499 | 0.034 | 0.048 | 0.038 | 0.052 | 0.010 |  |
| se | N | regime | none | 503 | 0.052 | 0.044 | 0.064 | 0.044 | 0.018 |  |
| se | N | role | N_0 | 201 | 0.025 | 0.000 | 0.030 | 0.060 | 0.005 |  |
| se | N | role | N_0.5 | 202 | 0.035 | 0.010 | 0.020 | 0.045 | 0.010 |  |
| se | N | role | N_0.8 | 213 | 0.033 | 0.033 | 0.052 | 0.047 | 0.014 |  |
| se | N | role | N_0.95 | 225 | 0.067 | 0.107 | 0.076 | 0.049 | 0.018 |  |
| se | N | role | N_1 | 226 | 0.106 | 0.115 | 0.088 | 0.035 | 0.031 |  |
| se | N | role | N_imp0.2 | 211 | 0.038 | 0.038 | 0.038 | 0.062 | 0.005 |  |
| se | N | role | N_imp0.4 | 213 | 0.042 | 0.019 | 0.042 | 0.052 | 0.000 |  |
| se | N | twin genus | 0 | 756 | 0.003 | 0.040 | 0.021 | 0.029 | 0.003 |  |
| se | N | twin genus | 1 | 735 | 0.099 | 0.056 | 0.080 | 0.071 | 0.022 |  |
| hifi | Q | MRCA with nearest allele genome | 0.002-0.01 | 97 | 0.361 | 0.371 | 0.485 | 0.392 | 0.454 |  |
| hifi | Q | MRCA with nearest allele genome | long ago (>=0.01) | 20 | 0.200 | 0.300 | 0.000 | 0.150 | 0.200 |  |
| hifi | Q | MRCA with nearest allele genome | no allele genome | 192 | 0.417 | 0.516 | 0.000 | 0.000 | 0.073 |  |
| hifi | Q | MRCA with nearest allele genome | recent (<0.002) | 267 | 0.543 | 0.588 | 0.933 | 0.513 | 0.573 |  |
| hifi | Q | depth | 0.5 | 143 | 0.469 | 0.517 | 0.517 | 0.329 | 0.371 |  |
| hifi | Q | depth | 1 | 159 | 0.453 | 0.484 | 0.522 | 0.296 | 0.384 |  |
| hifi | Q | depth | 2 | 136 | 0.507 | 0.544 | 0.493 | 0.287 | 0.338 |  |
| hifi | Q | depth | 4 | 138 | 0.406 | 0.529 | 0.522 | 0.326 | 0.399 |  |
| hifi | Q | distance to rep | 0.005-0.015 | 274 | 0.599 | 0.515 | 0.562 | 0.299 | 0.350 |  |
| hifi | Q | distance to rep | 0.015-0.03 | 148 | 0.041 | 0.426 | 0.500 | 0.392 | 0.230 |  |
| hifi | Q | distance to rep | <0.005 | 94 | 1.000 | 1.000 | 0.404 | 0.064 | 0.830 |  |
| hifi | Q | distance to rep | >=0.03 | 60 | 0.000 | 0.000 | 0.500 | 0.533 | 0.117 |  |
| hifi | Q | genes with a congener's segment | 0.1-0.3 | 190 | 0.463 | 0.505 | 0.542 | 0.358 | 0.263 |  |
| hifi | Q | genes with a congener's segment | <=0.1 | 152 | 0.441 | 0.454 | 0.500 | 0.263 | 0.467 |  |
| hifi | Q | genes with a congener's segment | >0.3 | 82 | 0.390 | 0.622 | 0.512 | 0.232 | 0.024 |  |
| hifi | Q | genes with a congener's segment | none | 152 | 0.507 | 0.539 | 0.493 | 0.336 | 0.605 |  |
| hifi | Q | regime | high | 192 | 0.505 | 0.583 | 0.557 | 0.292 | 0.266 |  |
| hifi | Q | regime | low | 192 | 0.375 | 0.438 | 0.516 | 0.302 | 0.380 |  |
| hifi | Q | regime | none | 192 | 0.495 | 0.531 | 0.469 | 0.333 | 0.474 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 314 | 0.497 | 0.525 | 0.901 | 0.516 | 0.611 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | no such site | 2 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 198 | 0.434 | 0.500 | 0.000 | 0.000 | 0.071 |  |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | partly | 62 | 0.323 | 0.516 | 0.210 | 0.258 | 0.113 |  |
| hifi | Q | role | Q_deep | 101 | 0.218 | 0.188 | 0.455 | 0.426 | 0.307 |  |
| hifi | Q | role | Q_ils | 74 | 0.149 | 0.122 | 0.473 | 0.405 | 0.027 |  |
| hifi | Q | role | Q_imp | 136 | 0.397 | 0.559 | 0.544 | 0.309 | 0.059 |  |
| hifi | Q | role | Q_lone | 92 | 0.489 | 0.620 | 0.174 | 0.141 | 0.283 |  |
| hifi | Q | role | Q_near | 100 | 0.620 | 0.640 | 1.000 | 0.500 | 0.900 |  |
| hifi | Q | role | Q_rep | 73 | 0.959 | 1.000 | 0.342 | 0.000 | 0.795 |  |
| hifi | Q | twin genus | 0 | 288 | 0.326 | 0.611 | 0.507 | 0.271 | 0.375 |  |
| hifi | Q | twin genus | 1 | 288 | 0.590 | 0.424 | 0.521 | 0.347 | 0.372 |  |
| hifi | N | depth | 0.5 | 129 | 0.047 | 0.054 | 0.078 | 0.085 | 0.047 |  |
| hifi | N | depth | 1 | 135 | 0.044 | 0.030 | 0.052 | 0.052 | 0.030 |  |
| hifi | N | depth | 2 | 122 | 0.074 | 0.066 | 0.041 | 0.033 | 0.057 |  |
| hifi | N | depth | 4 | 127 | 0.039 | 0.055 | 0.031 | 0.031 | 0.071 |  |
| hifi | N | genes with a segment of S or a congener | <=0.3 | 289 | 0.045 | 0.038 | 0.038 | 0.035 | 0.038 |  |
| hifi | N | genes with a segment of S or a congener | >0.3 | 98 | 0.000 | 0.000 | 0.020 | 0.082 | 0.000 |  |
| hifi | N | genes with a segment of S or a congener | none | 126 | 0.103 | 0.119 | 0.103 | 0.063 | 0.119 |  |
| hifi | N | regime | high | 166 | 0.030 | 0.012 | 0.036 | 0.072 | 0.012 |  |
| hifi | N | regime | low | 168 | 0.048 | 0.054 | 0.030 | 0.012 | 0.054 |  |
| hifi | N | regime | none | 179 | 0.073 | 0.084 | 0.084 | 0.067 | 0.084 |  |
| hifi | N | role | N_0 | 75 | 0.027 | 0.000 | 0.013 | 0.013 | 0.000 |  |
| hifi | N | role | N_0.5 | 74 | 0.000 | 0.000 | 0.000 | 0.081 | 0.000 |  |
| hifi | N | role | N_0.8 | 65 | 0.046 | 0.046 | 0.000 | 0.031 | 0.000 |  |
| hifi | N | role | N_0.95 | 73 | 0.110 | 0.137 | 0.096 | 0.055 | 0.041 |  |
| hifi | N | role | N_1 | 72 | 0.181 | 0.181 | 0.208 | 0.042 | 0.319 |  |
| hifi | N | role | N_imp0.2 | 77 | 0.000 | 0.000 | 0.013 | 0.078 | 0.000 |  |
| hifi | N | role | N_imp0.4 | 77 | 0.000 | 0.000 | 0.026 | 0.052 | 0.000 |  |
| hifi | N | twin genus | 0 | 259 | 0.004 | 0.073 | 0.012 | 0.081 | 0.039 |  |
| hifi | N | twin genus | 1 | 254 | 0.098 | 0.028 | 0.091 | 0.020 | 0.063 |  |

## L3: the presence models (trained on the train world, called on the test world)

Models by read type and feature set (base: the default set without ancestry, alleles and polymorphic; _on_noallele / _on_shuffled: the default model on the test world's runs without allele scores and with the shuffled table). Q and N of 'all test rows': present and absent taxa of every row.
| model | stratum | Q | N | auc | tpr_at_5pct_fpr | sensitivity_at_knob | fpr_at_knob | f1_at_knob |
|---|---|---|---|---|---|---|---|---|
| hifi_base | target rows | 576 | 513 | 0.872 | 0.564 | 0.682 | 0.131 |  |
| hifi_base | hard: Q_deep vs N_0.8/0.95 | 101 | 138 | 0.683 | 0.277 | 0.475 | 0.232 |  |
| hifi_base | twin genera | 288 | 254 | 0.843 | 0.476 | 0.649 | 0.134 |  |
| hifi_base | no allele genomes | 192 | 172 | 0.848 | 0.562 | 0.667 | 0.163 |  |
| hifi_base | recombination high | 192 | 166 | 0.935 | 0.682 | 0.891 | 0.205 |  |
| hifi_base | all test rows | 896 | 940 | 0.950 | - | 0.796 | 0.071 | 0.851 |
| hifi_default | target rows | 576 | 513 | 0.939 | 0.804 | 0.849 | 0.099 |  |
| hifi_default | hard: Q_deep vs N_0.8/0.95 | 101 | 138 | 0.845 | 0.535 | 0.723 | 0.196 |  |
| hifi_default | twin genera | 288 | 254 | 0.928 | 0.774 | 0.872 | 0.142 |  |
| hifi_default | no allele genomes | 192 | 172 | 0.856 | 0.562 | 0.714 | 0.163 |  |
| hifi_default | recombination high | 192 | 166 | 0.974 | 0.885 | 0.906 | 0.096 |  |
| hifi_default | all test rows | 896 | 940 | 0.976 | - | 0.903 | 0.054 | 0.921 |
| pe_alleles | target rows | 1724 | 1493 | 0.922 | 0.733 | 0.794 | 0.101 |  |
| pe_alleles | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.801 | 0.503 | 0.606 | 0.129 |  |
| pe_alleles | twin genera | 860 | 735 | 0.929 | 0.728 | 0.802 | 0.112 |  |
| pe_alleles | no allele genomes | 573 | 496 | 0.823 | 0.487 | 0.599 | 0.115 |  |
| pe_alleles | recombination high | 575 | 492 | 0.938 | 0.763 | 0.810 | 0.085 |  |
| pe_alleles | all test rows | 2684 | 3724 | 0.977 | - | 0.867 | 0.042 | 0.901 |
| pe_ancestry+alleles | target rows | 1724 | 1493 | 0.917 | 0.738 | 0.791 | 0.094 |  |
| pe_ancestry+alleles | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.791 | 0.481 | 0.594 | 0.122 |  |
| pe_ancestry+alleles | twin genera | 860 | 735 | 0.925 | 0.744 | 0.785 | 0.086 |  |
| pe_ancestry+alleles | no allele genomes | 573 | 496 | 0.816 | 0.508 | 0.579 | 0.083 |  |
| pe_ancestry+alleles | recombination high | 575 | 492 | 0.934 | 0.763 | 0.809 | 0.089 |  |
| pe_ancestry+alleles | all test rows | 2684 | 3724 | 0.975 | - | 0.865 | 0.039 | 0.901 |
| pe_ancestry | target rows | 1724 | 1493 | 0.844 | 0.475 | 0.668 | 0.157 |  |
| pe_ancestry | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.642 | 0.165 | 0.471 | 0.258 |  |
| pe_ancestry | twin genera | 860 | 735 | 0.837 | 0.459 | 0.641 | 0.137 |  |
| pe_ancestry | no allele genomes | 573 | 496 | 0.818 | 0.469 | 0.668 | 0.177 |  |
| pe_ancestry | recombination high | 575 | 492 | 0.883 | 0.496 | 0.810 | 0.209 |  |
| pe_ancestry | all test rows | 2684 | 3724 | 0.955 | - | 0.787 | 0.064 | 0.839 |
| pe_base | target rows | 1724 | 1493 | 0.845 | 0.486 | 0.668 | 0.151 |  |
| pe_base | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.641 | 0.194 | 0.452 | 0.249 |  |
| pe_base | twin genera | 860 | 735 | 0.832 | 0.455 | 0.642 | 0.148 |  |
| pe_base | no allele genomes | 573 | 496 | 0.823 | 0.476 | 0.661 | 0.149 |  |
| pe_base | recombination high | 575 | 492 | 0.885 | 0.506 | 0.800 | 0.177 |  |
| pe_base | all test rows | 2684 | 3724 | 0.956 | - | 0.785 | 0.062 | 0.839 |
| pe_default | target rows | 1724 | 1493 | 0.918 | 0.734 | 0.795 | 0.094 |  |
| pe_default | hard: Q_deep vs N_0.8/0.95 | 310 | 442 | 0.791 | 0.474 | 0.594 | 0.118 |  |
| pe_default | twin genera | 860 | 735 | 0.928 | 0.748 | 0.795 | 0.095 |  |
| pe_default | no allele genomes | 573 | 496 | 0.818 | 0.501 | 0.593 | 0.105 |  |
| pe_default | recombination high | 575 | 492 | 0.935 | 0.757 | 0.807 | 0.085 |  |
| pe_default | all test rows | 2684 | 3724 | 0.975 | - | 0.868 | 0.039 | 0.903 |
| pe_default_on_noallele | target rows | 1714 | 1465 | 0.913 | 0.723 | 0.780 | 0.092 |  |
| pe_default_on_noallele | hard: Q_deep vs N_0.8/0.95 | 307 | 436 | 0.792 | 0.433 | 0.593 | 0.115 |  |
| pe_default_on_noallele | twin genera | 850 | 709 | 0.919 | 0.713 | 0.765 | 0.090 |  |
| pe_default_on_noallele | no allele genomes | 573 | 496 | 0.818 | 0.494 | 0.593 | 0.099 |  |
| pe_default_on_noallele | recombination high | 572 | 486 | 0.932 | 0.731 | 0.792 | 0.088 |  |
| pe_default_on_noallele | all test rows | 2674 | 3876 | 0.974 | - | 0.859 | 0.037 | 0.898 |
| pe_default_on_shuffled | target rows | 1716 | 1467 | 0.829 | 0.484 | 0.595 | 0.106 |  |
| pe_default_on_shuffled | hard: Q_deep vs N_0.8/0.95 | 307 | 437 | 0.630 | 0.130 | 0.342 | 0.140 |  |
| pe_default_on_shuffled | twin genera | 852 | 711 | 0.824 | 0.481 | 0.576 | 0.110 |  |
| pe_default_on_shuffled | no allele genomes | 573 | 496 | 0.818 | 0.499 | 0.595 | 0.101 |  |
| pe_default_on_shuffled | recombination high | 573 | 487 | 0.876 | 0.538 | 0.649 | 0.094 |  |
| pe_default_on_shuffled | all test rows | 2676 | 3864 | 0.951 | - | 0.740 | 0.042 | 0.822 |

The models' calls on the target rows by bin: 'called' the share called at the model's knob (Q: sensitivity; N: false-positive rate), 'at 5% FPR' the share above the score 5% of the model's N rows exceed.
| set | class | by | level | rows | base called | ancestry called | alleles called | default called | base at 5% FPR | ancestry at 5% FPR | alleles at 5% FPR | default at 5% FPR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pe | Q | MRCA with nearest allele genome | 0.002-0.01 | 278 | 0.622 | 0.626 | 0.759 | 0.781 | 0.435 | 0.428 | 0.698 | 0.712 |
| pe | Q | MRCA with nearest allele genome | long ago (>=0.01) | 64 | 0.312 | 0.281 | 0.281 | 0.281 | 0.094 | 0.094 | 0.219 | 0.203 |
| pe | Q | MRCA with nearest allele genome | no allele genome | 573 | 0.661 | 0.668 | 0.599 | 0.593 | 0.483 | 0.475 | 0.496 | 0.487 |
| pe | Q | MRCA with nearest allele genome | recent (<0.002) | 809 | 0.716 | 0.713 | 0.984 | 0.984 | 0.536 | 0.522 | 0.954 | 0.959 |
| pe | Q | depth | 10 | 435 | 0.699 | 0.706 | 0.795 | 0.805 | 0.520 | 0.492 | 0.743 | 0.752 |
| pe | Q | depth | 100 | 426 | 0.669 | 0.676 | 0.812 | 0.812 | 0.523 | 0.507 | 0.765 | 0.758 |
| pe | Q | depth | 3 | 439 | 0.624 | 0.613 | 0.749 | 0.752 | 0.412 | 0.403 | 0.663 | 0.658 |
| pe | Q | depth | 30 | 424 | 0.679 | 0.679 | 0.818 | 0.814 | 0.491 | 0.500 | 0.764 | 0.771 |
| pe | Q | distance to rep | 0.005-0.015 | 790 | 0.720 | 0.724 | 0.843 | 0.832 | 0.497 | 0.476 | 0.761 | 0.751 |
| pe | Q | distance to rep | 0.015-0.03 | 433 | 0.573 | 0.575 | 0.732 | 0.755 | 0.339 | 0.323 | 0.663 | 0.679 |
| pe | Q | distance to rep | <0.005 | 298 | 0.990 | 0.987 | 0.990 | 0.993 | 0.946 | 0.960 | 0.966 | 0.977 |
| pe | Q | distance to rep | >=0.03 | 203 | 0.192 | 0.182 | 0.443 | 0.448 | 0.079 | 0.084 | 0.433 | 0.433 |
| pe | Q | genes with a congener's segment | 0.1-0.3 | 570 | 0.746 | 0.761 | 0.811 | 0.805 | 0.512 | 0.519 | 0.737 | 0.739 |
| pe | Q | genes with a congener's segment | <=0.1 | 454 | 0.593 | 0.597 | 0.670 | 0.685 | 0.447 | 0.436 | 0.630 | 0.643 |
| pe | Q | genes with a congener's segment | >0.3 | 251 | 0.709 | 0.689 | 0.845 | 0.857 | 0.506 | 0.470 | 0.777 | 0.797 |
| pe | Q | genes with a congener's segment | none | 449 | 0.621 | 0.610 | 0.869 | 0.860 | 0.481 | 0.461 | 0.808 | 0.786 |
| pe | Q | regime | high | 575 | 0.800 | 0.810 | 0.810 | 0.807 | 0.572 | 0.567 | 0.739 | 0.744 |
| pe | Q | regime | low | 574 | 0.585 | 0.584 | 0.688 | 0.702 | 0.420 | 0.409 | 0.641 | 0.653 |
| pe | Q | regime | none | 575 | 0.617 | 0.610 | 0.882 | 0.877 | 0.466 | 0.449 | 0.819 | 0.805 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 933 | 0.673 | 0.673 | 0.965 | 0.969 | 0.492 | 0.477 | 0.936 | 0.945 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | no such site | 4 | 1.000 | 1.000 | 1.000 | 1.000 | 0.750 | 1.000 | 1.000 | 1.000 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 589 | 0.669 | 0.674 | 0.606 | 0.601 | 0.491 | 0.479 | 0.503 | 0.489 |
| pe | Q | rep-lineage sites it lacks that the alleles cover | partly | 198 | 0.631 | 0.621 | 0.540 | 0.551 | 0.439 | 0.444 | 0.460 | 0.465 |
| pe | Q | role | Q_deep | 310 | 0.452 | 0.471 | 0.606 | 0.594 | 0.226 | 0.194 | 0.513 | 0.510 |
| pe | Q | role | Q_ils | 200 | 0.430 | 0.425 | 0.510 | 0.520 | 0.125 | 0.120 | 0.430 | 0.420 |
| pe | Q | role | Q_imp | 402 | 0.694 | 0.687 | 0.863 | 0.873 | 0.490 | 0.473 | 0.803 | 0.813 |
| pe | Q | role | Q_lone | 277 | 0.740 | 0.736 | 0.740 | 0.740 | 0.596 | 0.585 | 0.661 | 0.639 |
| pe | Q | role | Q_near | 334 | 0.722 | 0.728 | 0.991 | 0.991 | 0.578 | 0.572 | 0.967 | 0.973 |
| pe | Q | role | Q_rep | 201 | 0.995 | 0.985 | 0.970 | 0.975 | 0.935 | 0.955 | 0.945 | 0.970 |
| pe | Q | twin genus | 0 | 864 | 0.693 | 0.696 | 0.785 | 0.795 | 0.537 | 0.523 | 0.738 | 0.745 |
| pe | Q | twin genus | 1 | 860 | 0.642 | 0.641 | 0.802 | 0.795 | 0.435 | 0.427 | 0.728 | 0.723 |
| pe | N | depth | 10 | 382 | 0.160 | 0.168 | 0.115 | 0.110 | 0.065 | 0.058 | 0.073 | 0.068 |
| pe | N | depth | 100 | 377 | 0.133 | 0.138 | 0.064 | 0.072 | 0.040 | 0.040 | 0.029 | 0.024 |
| pe | N | depth | 3 | 353 | 0.164 | 0.167 | 0.125 | 0.113 | 0.057 | 0.051 | 0.057 | 0.074 |
| pe | N | depth | 30 | 381 | 0.150 | 0.157 | 0.102 | 0.084 | 0.039 | 0.052 | 0.042 | 0.037 |
| pe | N | genes with a segment of S or a congener | <=0.3 | 860 | 0.167 | 0.181 | 0.081 | 0.085 | 0.060 | 0.063 | 0.034 | 0.043 |
| pe | N | genes with a segment of S or a congener | >0.3 | 266 | 0.068 | 0.060 | 0.075 | 0.060 | 0.008 | 0.015 | 0.038 | 0.034 |
| pe | N | genes with a segment of S or a congener | none | 367 | 0.174 | 0.172 | 0.166 | 0.142 | 0.057 | 0.046 | 0.098 | 0.079 |
| pe | N | regime | high | 492 | 0.177 | 0.209 | 0.085 | 0.085 | 0.079 | 0.081 | 0.043 | 0.047 |
| pe | N | regime | low | 495 | 0.119 | 0.111 | 0.069 | 0.069 | 0.028 | 0.032 | 0.026 | 0.036 |
| pe | N | regime | none | 506 | 0.158 | 0.152 | 0.148 | 0.128 | 0.043 | 0.038 | 0.081 | 0.067 |
| pe | N | role | N_0 | 195 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.5 | 203 | 0.010 | 0.020 | 0.015 | 0.020 | 0.000 | 0.000 | 0.000 | 0.000 |
| pe | N | role | N_0.8 | 215 | 0.158 | 0.177 | 0.088 | 0.074 | 0.028 | 0.019 | 0.037 | 0.033 |
| pe | N | role | N_0.95 | 227 | 0.335 | 0.335 | 0.167 | 0.159 | 0.119 | 0.115 | 0.101 | 0.110 |
| pe | N | role | N_1 | 228 | 0.395 | 0.399 | 0.285 | 0.285 | 0.175 | 0.180 | 0.145 | 0.140 |
| pe | N | role | N_imp0.2 | 213 | 0.038 | 0.052 | 0.052 | 0.038 | 0.005 | 0.005 | 0.019 | 0.023 |
| pe | N | role | N_imp0.4 | 212 | 0.075 | 0.071 | 0.071 | 0.057 | 0.005 | 0.014 | 0.033 | 0.028 |
| pe | N | twin genus | 0 | 758 | 0.154 | 0.177 | 0.091 | 0.094 | 0.059 | 0.066 | 0.050 | 0.059 |
| pe | N | twin genus | 1 | 735 | 0.148 | 0.137 | 0.112 | 0.095 | 0.041 | 0.034 | 0.050 | 0.041 |
| hifi | Q | MRCA with nearest allele genome | 0.002-0.01 | 97 | 0.557 |  |  | 0.825 | 0.412 |  |  | 0.804 |
| hifi | Q | MRCA with nearest allele genome | long ago (>=0.01) | 20 | 0.550 |  |  | 0.300 | 0.450 |  |  | 0.200 |
| hifi | Q | MRCA with nearest allele genome | no allele genome | 192 | 0.667 |  |  | 0.714 | 0.557 |  |  | 0.599 |
| hifi | Q | MRCA with nearest allele genome | recent (<0.002) | 267 | 0.749 |  |  | 0.996 | 0.633 |  |  | 0.996 |
| hifi | Q | depth | 0.5 | 143 | 0.650 |  |  | 0.839 | 0.538 |  |  | 0.776 |
| hifi | Q | depth | 1 | 159 | 0.679 |  |  | 0.830 | 0.585 |  |  | 0.799 |
| hifi | Q | depth | 2 | 136 | 0.713 |  |  | 0.875 | 0.596 |  |  | 0.846 |
| hifi | Q | depth | 4 | 138 | 0.688 |  |  | 0.855 | 0.536 |  |  | 0.797 |
| hifi | Q | distance to rep | 0.005-0.015 | 274 | 0.741 |  |  | 0.931 | 0.588 |  |  | 0.861 |
| hifi | Q | distance to rep | 0.015-0.03 | 148 | 0.541 |  |  | 0.743 | 0.405 |  |  | 0.696 |
| hifi | Q | distance to rep | <0.005 | 94 | 1.000 |  |  | 1.000 | 1.000 |  |  | 1.000 |
| hifi | Q | distance to rep | >=0.03 | 60 | 0.267 |  |  | 0.500 | 0.167 |  |  | 0.500 |
| hifi | Q | genes with a congener's segment | 0.1-0.3 | 190 | 0.805 |  |  | 0.895 | 0.668 |  |  | 0.847 |
| hifi | Q | genes with a congener's segment | <=0.1 | 152 | 0.520 |  |  | 0.691 | 0.428 |  |  | 0.645 |
| hifi | Q | genes with a congener's segment | >0.3 | 82 | 0.817 |  |  | 0.890 | 0.671 |  |  | 0.841 |
| hifi | Q | genes with a congener's segment | none | 152 | 0.618 |  |  | 0.928 | 0.513 |  |  | 0.888 |
| hifi | Q | regime | high | 192 | 0.891 |  |  | 0.906 | 0.760 |  |  | 0.870 |
| hifi | Q | regime | low | 192 | 0.516 |  |  | 0.703 | 0.417 |  |  | 0.651 |
| hifi | Q | regime | none | 192 | 0.641 |  |  | 0.938 | 0.516 |  |  | 0.891 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | all (>=0.9) | 314 | 0.691 |  |  | 0.987 | 0.564 |  |  | 0.987 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | no such site | 2 | 1.000 |  |  | 1.000 | 1.000 |  |  | 1.000 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | none (<0.1) | 198 | 0.677 |  |  | 0.722 | 0.566 |  |  | 0.611 |
| hifi | Q | rep-lineage sites it lacks that the alleles cover | partly | 62 | 0.645 |  |  | 0.548 | 0.548 |  |  | 0.484 |
| hifi | Q | role | Q_deep | 101 | 0.475 |  |  | 0.723 | 0.307 |  |  | 0.624 |
| hifi | Q | role | Q_ils | 74 | 0.405 |  |  | 0.703 | 0.243 |  |  | 0.635 |
| hifi | Q | role | Q_imp | 136 | 0.750 |  |  | 0.890 | 0.625 |  |  | 0.846 |
| hifi | Q | role | Q_lone | 92 | 0.728 |  |  | 0.804 | 0.609 |  |  | 0.750 |
| hifi | Q | role | Q_near | 100 | 0.740 |  |  | 0.990 | 0.650 |  |  | 0.990 |
| hifi | Q | role | Q_rep | 73 | 0.986 |  |  | 0.959 | 0.959 |  |  | 0.959 |
| hifi | Q | twin genus | 0 | 288 | 0.715 |  |  | 0.826 | 0.625 |  |  | 0.795 |
| hifi | Q | twin genus | 1 | 288 | 0.649 |  |  | 0.872 | 0.503 |  |  | 0.812 |
| hifi | N | depth | 0.5 | 129 | 0.116 |  |  | 0.085 | 0.047 |  |  | 0.054 |
| hifi | N | depth | 1 | 135 | 0.126 |  |  | 0.096 | 0.044 |  |  | 0.052 |
| hifi | N | depth | 2 | 122 | 0.148 |  |  | 0.115 | 0.082 |  |  | 0.066 |
| hifi | N | depth | 4 | 127 | 0.134 |  |  | 0.102 | 0.031 |  |  | 0.031 |
| hifi | N | genes with a segment of S or a congener | <=0.3 | 289 | 0.145 |  |  | 0.069 | 0.055 |  |  | 0.017 |
| hifi | N | genes with a segment of S or a congener | >0.3 | 98 | 0.061 |  |  | 0.061 | 0.051 |  |  | 0.010 |
| hifi | N | genes with a segment of S or a congener | none | 126 | 0.151 |  |  | 0.198 | 0.040 |  |  | 0.159 |
| hifi | N | regime | high | 166 | 0.205 |  |  | 0.096 | 0.078 |  |  | 0.018 |
| hifi | N | regime | low | 168 | 0.071 |  |  | 0.036 | 0.036 |  |  | 0.012 |
| hifi | N | regime | none | 179 | 0.117 |  |  | 0.162 | 0.039 |  |  | 0.117 |
| hifi | N | role | N_0 | 75 | 0.000 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_0.5 | 74 | 0.000 |  |  | 0.000 | 0.000 |  |  | 0.000 |
| hifi | N | role | N_0.8 | 65 | 0.092 |  |  | 0.154 | 0.000 |  |  | 0.062 |
| hifi | N | role | N_0.95 | 73 | 0.356 |  |  | 0.233 | 0.151 |  |  | 0.123 |
| hifi | N | role | N_1 | 72 | 0.361 |  |  | 0.236 | 0.125 |  |  | 0.167 |
| hifi | N | role | N_imp0.2 | 77 | 0.039 |  |  | 0.013 | 0.013 |  |  | 0.000 |
| hifi | N | role | N_imp0.4 | 77 | 0.078 |  |  | 0.078 | 0.065 |  |  | 0.013 |
| hifi | N | twin genus | 0 | 259 | 0.127 |  |  | 0.058 | 0.046 |  |  | 0.012 |
| hifi | N | twin genus | 1 | 254 | 0.134 |  |  | 0.142 | 0.055 |  |  | 0.091 |
