# Where protal's false positives come from at GTDB scale: the denoising ideas revisited, the sample's depth as a feature, and what GTDB's species definition can add

- **Date**: 2026-10-03.
- **Question** (the user's): the amplicon-denoising work of the last two days
  ([ideas](../2026-10-02-amplicon-denoising/README.md), [implementation](../2026-10-03-denoising-implementation/README.md))
  was less successful than hoped. Is there a better way to reduce false positives? And can GTDB's species definition
  and its cutoffs be exploited against false positives from congeners the database lacks?
- **Code**: protal `27423c6` (`audit-fixes`): nothing changed in protal. The trainer is `scripts/` of `27423c6` with
  two more feature sets ([`model_features_depth.patch`](model_features_depth.patch),
  [`model_features_genus.patch`](model_features_genus.patch)).
- **Data**: the r226 v5 run (`local/v5/`, git-ignored: `27423c6`, `--features normalized+adjacency+distance
  --call-mode fdr`, congener groups 0.25:2-5; 198 paired-end training samples of 1,000 to 10M read pairs, 78 test
  samples of 500 to 5M; 20% of species and 32 clades held out), with its models' predictions (species held out, and
  on the test set), and the v6 run (`local/v6/`, the defaults) for comparison. Another report of today compares the two
  runs' models ([r226 v5 and v6](../2026-10-03-r226-v5-v6-training/README.md)); this one looks at the errors.
- **Run**: in WSL with the protal-db-build env (Python 3.14, scikit-learn 1.9.1). [`fp_anatomy.py`](fp_anatomy.py) and
  [`fp_sources.py`](fp_sources.py) read the v5 tables and predictions (outputs [`fp_anatomy_output.txt`](fp_anatomy_output.txt),
  [`fp_sources_output.txt`](fp_sources_output.txt)); [`run_depth.sh`](run_depth.sh) and [`run_genus.sh`](run_genus.sh) add
  the columns ([`add_sample_depth.py`](add_sample_depth.py), [`add_genus_size.py`](add_genus_size.py)) and train the
  models of the experiment, [`summarize_depth.py`](summarize_depth.py) tabulates them
  ([`depth_feature_output.txt`](depth_feature_output.txt)). Models in WSL `~/fpdepth`, about a minute each.

## Summary

- **Half of the false positives are single reads that match the wrong species perfectly.** Of 764 paired-end false
  positives (v5 model, species held out plus the test set), 399 have one fragment, 348 of them at identity 0.97 or
  more. Their features equal those of true single-read species (identity 0.990 against 0.990, MAPQ 87 against 92).
  The only thing that separates them is the sample's depth: true singletons sit in samples of 50,000 pairs or fewer,
  false ones in samples of 100,000 or more. No per-taxon rule can touch them; a prior can.
- **The sample's depth is not a feature.** The forest sees no sample-level quantity; the knob curve supplies the depth
  from outside, one threshold per sample. Given as a feature (`sample_log_fragments`, exactly what `DepthKnobAt`
  reads), the forest halves the false positives at knob 0.5 (pe test 218 → 122, se 231 → 118), cuts the log loss
  (0.048 → 0.041, 0.058 → 0.048) and raises the test F1 at knob 0.5 by **0.004 (pe) and 0.009 (se)** over the
  distance set, and by 0.001-0.005 over its knob curve; PacBio +0.005 at 0.5, ONT unchanged. The knob curve then does
  harm (pe 0.9653 → 0.9555) and must be off when the feature is in. A one-line feature in protal and the trainer.
- **The denoising analogy reached the wrong population.** UNOISE and DADA2 judge a candidate against a more abundant
  parent. 83% of the false positives (454 of 546, training) have no present species in their genus at all; the parent
  is a species the database lacks, or a species of another genus. The 92 with a present congener are where the
  relatives features work, and where they conflict with real minor congeners.
- **A quarter of the false positives come from another genus** (159 of 764 from a present database species of another
  genus, at identity 0.99, false-positive rate 2.8-4.9% of such absent taxa against 0.3% beside a present congener;
  200 more from a held-out species of another genus). Marker-gene copies 99% identical across genera are not what
  marker genes do; the likely cause is a contaminated reference (a MAG's foreign contig) or a transferred gene, both
  detectable at build time with GTDB's taxonomy. The same taxa recur: 288 of the 764 false positives are 116 taxa
  called falsely in two or more samples.
- **The remaining third, clouds of 3 or more fragments beside a held-out congener, is the species boundary as
  protal's marker genes see it.** Their reads sit at identity 0.96-0.98 on conserved and fast genes alike, and the
  strains protal misses look the same (identity 0.963, excess share 0.5). GTDB's 95% ANI radius applies to the whole
  genome; the 168 marker genes compress a 6-10% genome divergence into 2-4%. Calibrating read identity per gene to
  genome ANI, and reporting a cloud below the species' radius as a novel species of the genus, is the one way the
  definition can be exploited for these; it needs an experiment with per-read gene data, not available offline.
- **Yesterday's calibrated calls (`--call-mode fdr`) are worse than the knob curve at r226** for all four read types
  (pe test 0.9573 against 0.9642). v5's database calls at them by default; `--call-mode curve` should be used.

## 1. What the r226 runs say about yesterday's ideas

Species held out / test at knob 0.5 / test at the knob curve, from the builds' training logs:

| model | v6 (`normalized+adjacency`, curve) | v5 (`+distance`, fdr) | v5 at its fdr target |
|---|---|---|---|
| pe | 0.9563 / 0.9542 / 0.9605 | 0.9625 / 0.9615 / 0.9642 | 0.9573 |
| se | 0.9508 / 0.9513 / — | 0.9579 / 0.9570 / 0.9593 | 0.9556 |
| pb | 0.9719 / 0.9714 / — | 0.9733 / 0.9714 / 0.9761 | 0.9711 |
| ont | 0.9690 / 0.9698 / — | 0.9716 / 0.9704 / 0.9696 | 0.9685 |

- The distance features gain 0.004-0.007 (pe), 0.006 (se), nothing for long reads: real, small. False positives with
  species held out: pe 608 → 546, se 663 → 596. The full relatives set retrained on v5's tables
  (`local/v5_relatives`) gets to 433 pe false positives and F1 0.9689 with species held out, yet
  [the v5/v6 report](../2026-10-03-r226-v5-v6-training/README.md) finds it no better than the distance set on the test set.
- The fdr calls lose 0.001-0.007 against the curve for every model. protal uses them by default when a model carries
  them, so v5's shipped database calls below its own curve.
- 78% of v5's remaining false positives have fewer than 10 fragments (336 of 433 in the relatives retraining, 502 of
  608 for v6).

## 2. Anatomy of the false positives

Paired-end, v5's model (`normalized+adjacency+distance`) at knob 0.5: 546 false positives with species held out
(56,761 rows, 11,529 present) and 218 on the test set (23,641 rows, 4,591 present), 764 in all
([`fp_anatomy_output.txt`](fp_anatomy_output.txt), [`fp_sources_output.txt`](fp_sources_output.txt)).

### By fragments and by source

| fragments | false positives | share |
|---|---:|---:|
| 1 | 399 | 52% |
| 2 | 96 | 13% |
| 3-9 | 122 | 16% |
| 10-99 | 100 | 13% |
| ≥ 100 | 47 | 6% |

The source of a false positive's reads, by the deepest rank it shares with the closest species simulated in the sample
(`meta_relative_rank`) and whether that species is in the database (`meta_novel_level` empty) or held out:

| the closest simulated species is | false positives | share | of such absent rows | false-positive rate | median identity |
|---|---:|---:|---:|---:|---:|
| a congener the database lacks | 351 | 46% | 21,564 | 1.6% | 0.973 |
| a congener the database has (spill-over from a present species) | 54 | 7% | 19,468 | **0.3%** | 0.986 |
| another genus of the family, in the database | 100 | 13% | 3,542 | 2.8% | 0.988 |
| another genus of the family, held out | 88 | 12% | 12,670 | 0.7% | 0.990 |
| another family of the order, in the database | 44 | 6% | 904 | **4.9%** | 0.990 |
| another family of the order, held out | 38 | 5% | 2,502 | 1.5% | 0.990 |
| class, phylum or domain only | 89 | 12% | 3,632 | 2.5% | 0.986 |

- **Spill-over from a present congener, the case the denoisers model, is the smallest class** and has the lowest
  false-positive rate: the forest already handles it (MAPQ, `congener_fit_share`, now the distance features).
- In the training set 454 of 546 false positives (83%) have **no present species in their genus**; 452 are the top
  species of their genus in the sample, with a median of one species of the genus holding records. There is no
  abundant parent for a skew feature to compare with, and no visible spray of low-identity reads over the siblings.
- **The false-positive rate is highest for taxa whose closest sample species is in another genus or family**, and
  those reads are at identity 0.99 (section 2.3).

### 2.1 Single fragments: perfect reads on the wrong species

399 false positives have one fragment. Quartiles of the 1-fragment rows (training and test together):

| | true, called (1,445) | true, missed (255) | false, called (399) | absent, not called (36,286) |
|---|---|---|---|---|
| identity | 0.983 / 0.990 / 0.997 | 0.960 / 0.974 / 0.986 | 0.980 / 0.990 / 0.995 | 0.910 / 0.927 / 0.950 |
| 1 − `excess_median` | 0.996 / 1 / 1 | 0.975 / 0.988 / 0.996 | 0.991 / 0.998 / 1 | 0.920 / 0.939 / 0.965 |
| `mean_mapq` | 15 / 92 / 107 | 7 / 26 / 90 | 14 / 87 / 104 | 8 / 78 / 89 |
| `lu_per_kb` | 45 / 80 / 111 | 0 / 12 / 38 | 30 / 60 / 100 | 0 / 0 / 0 |
| sample read pairs | 5k / 5k / 20k | 5k / 10k / 50k | **100k / 500k / 500k** | 500k / 1M / 5M |
| forest score | 0.83 / 0.93 / 0.97 | 0.06 / 0.20 / 0.36 | 0.65 / 0.80 / 0.92 | 0.002 / 0.002 / 0.003 |

- The false singletons are reads at identity 0.97 or more (348 of 399), with unique k-mers, at MAPQ 87: by every
  per-read measure they are a species' own reads. Their source is a species the database lacks (187 of 399, 78 of
  them congeners) or a present species of another genus (159).
- What separates them from true singletons is the sample: true single-read species occur in shallow samples, false
  ones in deep ones. At identity ≥ 0.99 in samples of 100,000 pairs or more there are 33 present single-fragment taxa
  and 1,119 absent ones; the forest calls 18 of the present and 163 of the absent. It scores a deep sample's singleton
  at 0.80 because it cannot see the depth; the knob curve is what rescues these (at the curve 81 of the test set's 125
  single-fragment false positives remain).
- A read-level rule does nothing: vetoing identity below 0.95 among taxa of fewer than 3 fragments removes 2 false
  positives and 5 true calls; below 0.97, 40 and 70. The identity-and-congener singleton rule of yesterday vetoed
  5,961 training rows, none of them present, but the forest called none of them anyway.

### 2.2 Three fragments or more: the species boundary on marker genes

269 false positives have 3 or more fragments, 243 of them beside a congener the database lacks. Quartiles against the
true calls and the misses with 3 or more fragments:

| | true, called (13,288) | true, missed (155) | false, called (269) |
|---|---|---|---|
| identity | 0.980 / 0.989 / 0.994 | 0.954 / 0.963 / 0.973 | 0.963 / 0.971 / 0.979 |
| `excess_high_share` | 0.004 / 0.03 / 0.08 | 0.30 / 0.50 / 0.62 | 0.24 / 0.38 / 0.50 |
| `low_mapq_share` | 0 / 0.04 / 0.26 | 0.26 / 0.67 / 0.87 | 0.04 / 0.31 / 0.70 |
| `congener_fit_share` | 0 / 0.01 / 0.13 | 0.14 / 0.39 / 0.62 | 0.01 / 0.20 / 0.44 |
| `conserved_hit_share` | 0.26 / 0.35 / 0.43 | 0.14 / 0.29 / 0.38 | 0.20 / 0.33 / 0.50 |
| `hit_gene_fraction` | 0.12 / 0.37 / 0.78 | 0.04 / 0.09 / 0.30 | 0.05 / 0.13 / 0.37 |

- The false clouds and the missed strains are one population: reads at identity 0.96-0.98, a third of them more
  divergent than their qualities explain, a third at low MAPQ. The forest's threshold runs through it; moving it
  trades one error for the other. This is why no classifier change gained anything in the F1 reports.
- The clouds are **not** concentrated on conserved genes (`conserved_hit_share` 0.33 against 0.35): the source's
  marker genes are 2-4% from the taxon's across the gene set. GTDB congeners are usually 80-94% ANI apart over the
  genome, but protal's 168 genes are conserved by selection: the build's own log says the genes' divergence between
  congeners correlates 0.76 with their within-species conservation factors, and a conserved gene of factor 0.5 shows
  half the genome's divergence. On marker genes, a sister species at 92% ANI and a strain at 96% ANI both read as
  "2-4% away".

### 2.3 Another genus: contaminated or transferred gene copies

159 false positives (21%) have as their closest sample species a present database species of another genus (100),
family (44), class (9) or phylum (6), and their reads are at identity 0.99. For 88 of them exactly one present species
of the sample shares that rank, the candidate source. Examples
([`fp_sources_output.txt`](fp_sources_output.txt)):

| false positive | candidate source | samples with the source | taxon has records | called | identity |
|---|---|---:|---:|---:|---:|
| s__Paratissierella segnis | s__Tissierella_B sp036844705 (family) | 5 | 3 | 3 | 0.996 |
| s__CAKPSH01 sp022795675 | s__UMGS1781 sp004554705 (family) | 6 | 1 | 1 | 1.000 |
| s__Chlorobaculum sp035560975 | s__Chlorobium sp013334755 (family) | 4 | 1 | 1 | 1.000 |
| s__Colivicinus sp900546285 | s__Catenibacterium sp900540665 (order) | 3 | 1 | 1 | 1.000 |
| s__Atopostipes sp036752205 | s__Abiotrophia defectiva (order) | 5 | 1 | 1 | 0.993 |

- A read of a Tissierella_B strain that fits Paratissierella's gene copy at 99.6% and the taxon's own copy worse
  means one of the two references carries a copy from the other lineage: a contaminating contig in a MAG (GTDB admits
  up to 10% CheckM contamination, and most r226 species are MAGs), or a transferred gene. Either way the copy is
  incongruent with the taxonomy, and the database can see that at build time without any sample.
- The same recurrence shows for congeners of held-out species: 288 of the 764 false positives are 116 taxa called in
  two or more samples (s__Methanococcoides methylutens_A in 9 of its 11 samples with records, s__Thermococcus_A
  sp001507935 in 6 of 8, s__Nitrosopumilus catalinensis in 6 of 7). Where a taxon is called falsely whenever a given
  relative is present, the cause is in the references, not in the sample.

## 3. Why the denoising analogy under-delivered

1. **The parent is not there.** UNOISE's skew and DADA2's Poisson test need the more abundant sequence that produced
   the candidate. In 83% of protal's false positives no species of the taxon's genus is present; the parent is a
   species the database lacks (46%) or one of another genus (47%). The relatives features had 92 training false
   positives to work on, and there they fight real minor congeners, which real samples have many of.
2. **The singletons are a prior problem.** DADA2 never calls singletons; protal must (71% of true single-read taxa are
   called). What decides a singleton is how many absent candidates the sample has, which is its depth, which the
   forest cannot see (section 4).
3. **The noise is at the species boundary.** The clouds of 3 or more fragments are divergent reads of a related
   organism, biological not technical, and on marker genes they overlap the strains GTDB puts inside the species.
   Error models from base qualities (DADA2's λ) cannot separate them; `excess_*` already measures what qualities
   explain.

## 4. Experiment: the sample's depth, and the genus's size, as features

`sample_log_fragments`: log10 of the sample's fragments over all its taxa, the number protal already computes to read
the knob curve (`profiler::DepthKnobAt`, `random_forest_cmdline.sample_depths`); here added to every row of the v5
tables and trained with the trainer of `27423c6` (64 trees, 512 leaves for short reads, 128 for long, seed 1,
`--depth-knobs`, no clade holdouts). `genus_species`: log10 of the number of species of the taxon's genus in the
database (from `internal_taxonomy.dmp`), a per-species constant standing in for how speciose the neighbourhood is.
Test-set F1 at knob 0.5 / at the model's knob curve / at the best test threshold, false positives and misses at 0.5
([`depth_feature_output.txt`](depth_feature_output.txt)):

| read type | features | species held out | test 0.5 / curve / best | FP / FN at 0.5 | 1-fragment FP | log loss |
|---|---|---:|---|---|---:|---:|
| pe | `na` (= v6) | 0.9565 | 0.9544 / 0.9591 / 0.9553 | 240 / 181 | 125 | 0.0563 |
| pe | `na` + depth | 0.9740 | 0.9634 / 0.9638 / 0.9647 | 133 / 201 | 48 | 0.0428 |
| pe | `na` + genus | 0.9566 | 0.9546 / 0.9612 / 0.9558 | 240 / 180 | 122 | 0.0559 |
| pe | `nad` (= v5) | 0.9625 | 0.9615 / 0.9642 / 0.9624 | 218 / 139 | 108 | 0.0484 |
| pe | `nad` + depth | 0.9747 | **0.9653** / 0.9555 / 0.9664 | **122** / 194 | **39** | **0.0411** |
| pe | `nad` + depth + genus | 0.9742 | 0.9652 / 0.9627 / 0.9665 | 121 / 196 | 42 | 0.0416 |
| se | `na` | 0.9512 | 0.9499 / 0.9585 / 0.9518 | 260 / 193 | 151 | 0.0673 |
| se | `na` + depth | 0.9696 | 0.9615 / 0.9577 / 0.9640 | 131 / 211 | 51 | 0.0517 |
| se | `nad` | 0.9579 | 0.9570 / 0.9593 / 0.9578 | 231 / 158 | 135 | 0.0584 |
| se | `nad` + depth | 0.9712 | **0.9660** / 0.9640 / 0.9678 | **118** / 185 | **43** | **0.0475** |
| pb | `nad` | 0.9733 | 0.9714 / 0.9761 / 0.9732 | 38 / 60 | 21 | 0.0872 |
| pb | `nad` + depth | 0.9753 | 0.9760 / 0.9773 / 0.9776 | 26 / 56 | 12 | 0.0815 |
| ont | `nad` | 0.9716 | 0.9704 / 0.9696 / 0.9720 | 50 / 54 | 29 | 0.0864 |
| ont | `nad` + depth | 0.9737 | 0.9707 / 0.9685 / 0.9728 | 34 / 68 | 11 | 0.0805 |

(`na` = `normalized+adjacency`, `nad` = `normalized+adjacency+distance`; about ±0.003 of test F1 is noise.)

- **The depth feature halves the false positives at knob 0.5** (pe 218 → 122, se 231 → 118, pb 38 → 26, ont 50 → 34),
  mostly single fragments (pe 108 → 39), at the cost of a third more misses for short reads. Net: test F1 at a fixed
  knob +0.004 (pe), +0.009 (se), +0.005 (pb), 0 (ont) over the distance set, and the log loss down 11-19%. It also
  does what the knob curve did and more: +0.001 (pe) and +0.007 (se) over the distance set's curve, the number protal
  calls with today.
- **The knob curve must be off when the feature is in.** Fitted on top of it, the curve still finds points that gain
  on the training samples (held out 0.9747 → 0.9760) and loses on the test set (pe 0.9653 → 0.9555, 300 misses):
  the forest has already used the depth, and the curve corrects twice.
- **Genus size gains nothing** beyond noise (+0.002 at best, 0 on top of depth and distance). The forest's score
  does not need to know how speciose the genus is; dropped.
- **Caveat, the same as the curve's**: the forest cannot extrapolate past the deepest training sample (10M pairs,
  6 samples). Real samples go to 20-50M pairs; the design should have a point there (`DEPTH:SAMPLES` supports it),
  or the feature should be capped at the deepest trained depth, which keeps the curve's last-point behaviour.
- **Implementation**: `TaxonFeatures` writes the sample's log10 fragments (already known for the knob) as a feature of
  every taxon; `model_features.SAMPLE_FEATURES` in the default set; the trainer skips the knob curve, or fits it only
  where a point gains 0.002 (suggestion 1 of [the v4 report](../2026-10-03-r226-v4-v2-rerun/README.md)), when the
  feature is present; a model's header says which. A table of an older protal lacks the column, so the trainer must
  fall back as it does for the distance features.

## 5. GTDB's species definition, and what protal can take from it

**What GTDB does** ([methods](https://gtdb.ecogenomic.org/methods), [FAQ](https://gtdb.ecogenomic.org/faq), Parks et
al. 2020): one representative genome per species; a genome belongs to the species if its ANI to the representative is
within the species' **circumscription radius, 95% by default and up to 97% where that keeps an existing name**, with an
alignment fraction of at least 50% (65% before r207); species whose representatives are more than 97% ANI apart are
synonyms; ANI and AF by skani since r220. The per-species radius and the cluster's mean and minimum intra-species ANI
and AF are in the release's `sp_clusters_r226.tsv` (auxiliary files), which the build does not read; the metadata the
build reads gives only the representative and the taxonomy.

**What follows for protal's false positives:**

1. **The boundary is a genome-wide ANI; protal measures marker-gene identity.** The 168 genes are the most conserved
   part of the genome (the build's congener comparison: between-species divergence tracks the within-species
   conservation factors, correlation 0.76). A read at identity 0.975 on a gene of factor 0.5 implies a genome about
   95% ANI from the reference; on a gene of factor 1.2, about 98%. Pooled per taxon, `1 − (1 − identity + quality
   excess) / factor_g` over the reads is an **estimate of the organism's ANI to the reference**, in the units of the
   species definition. The median read identity, which the forest sees, is not: it mixes the genes. This is what
   could sharpen the 3-plus-fragment clouds against the missed strains (section 2.2), the one population where every
   present feature overlaps. It needs the per-read gene and identity, so an offline test on the tables is not
   possible; it is a feature in `Profiler.h` (the factors are loaded already) and an r226 run.
2. **The radius per species, from `sp_clusters`.** A species whose radius is below 95% has a described neighbour
   within 95% ANI; between 95 and 97% named sisters are inseparable on marker genes by construction. protal could
   store the radius and the nearest representative's ANI per species (the build already compares each species' genes
   with 4 congeners, `CompareCongeners`; a Mash sketch of the whole representatives would give the ANI to the nearest
   congener for all) and: (a) feed `nearest_congener_ani` to the forest, a prior on how likely a 2-4% cloud is a
   sister rather than a strain (genus size, the crude version, gained nothing, so the prior would have to be sharper
   than that); (b) report a cloud whose estimated ANI (item 1) falls below the radius as **"novel species of genus G,
   nearest X at ~ANI a"** instead of as X. That turns a false species call into a correct genus-level call with the
   information the user actually has, and it is what GTDB-Tk does for a genome outside every radius.
3. **Species complexes.** Pairs of representatives within 97% ANI (radius reduced, or names kept) spill into each
   other systematically. protal's `relative_distance` sees such a pair when both are in a sample; a build-time table
   of close pairs would let protal report them jointly when the reads do not separate them, instead of one false and
   one true call.
4. **The alignment fraction has an analogue** in `hit_gene_fraction` and `gene_presence_ratio`; nothing to add.
5. **The taxonomy against the references: an incongruence scan at build time.** For each gene copy, the nearest
   other species' copy (the k-mer index knows which species share each k-mer; `unique_kmers.tsv` is the per-species
   count of that pass). A copy whose nearest copy is in another genus or family at 97% or more identity, while its
   congeners' copies are far, is contamination or a transfer; mark it foreign for that species (protal has the
   foreign-gene mechanism: `--keep_foreign_genes`, `Taxon::ForeignGene`), so that reads on it are mapped but do not
   count as evidence of the species. Upper bound: the 159 false positives from another genus's present species (21%),
   and part of the 200 from another genus's held-out species where a database relative carries the same copy. The
   scan also finds GTDB's own contaminated references, worth reporting back.

## Suggestions, most useful first

1. **`sample_log_fragments` as a feature**, knob curve off with it, depths of 20-50M pairs in the design: +0.004 to
   +0.009 test F1 at a fixed knob, false positives halved, from this experiment. One day including an r226 run.
2. **`--call-mode curve`** (the default) for the next r226 build; v5's fdr calls are below its curve for every model.
3. **The gene-copy incongruence scan at build time** (item 5.5): finds the cross-genus sources, which are 21-47% of
   the false positives and recur per reference; also reports suspect GTDB references.
4. **The ANI estimate by gene conservation** (item 5.1) as a feature, and the **novel-species report** below the
   species' radius (item 5.2b) with `sp_clusters_r226.tsv` read by the build: for the 3-plus-fragment clouds, which
   no present feature separates from missed strains; an r226 run decides.
5. **A multi-sample view**: taxa called in many samples of a run always beside the same relative are references'
   artefacts; a run-level report (`unreported_species.tsv` has the shape) would flag them for users today.

**Not worth pursuing**: genus size as a prior (nothing); per-taxon identity or ANI rules against singletons (2 false
positives per 5 true calls); more abundance-skew features (the parent is absent in 83% of cases).

No protal behaviour changed, so the website needs no update.
