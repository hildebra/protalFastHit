# Against false positives: the sample's depth and the divergence features, suspect gene copies, defaults, a run-level check

- **Date**: 2026-10-03.
- **Branch**: `fp-anatomy` (worktree `../protal-fp`), from `audit-fixes` at `4cae2f4`.
- **What**: the user's four ideas evaluated against the code and the r226 v5 tables, and the agreed items of
  [the anatomy report](../2026-10-03-false-positive-anatomy/README.md) implemented: the sample's depth as a feature
  (default, no knob curve with it), four divergence features (by gene conservation, by codon position, lost mates),
  suspect gene copies found at build time and left out of the evidence, `--fdr` and the singleton rule off by default,
  knob points that must gain, `scripts/recurrent_calls.py`.
- **Round two** (the same day, below): fragments before the filters and the EM's fragments, failed candidates (the
  ZF tag and minimal unmapped records), species priors from GTDB, a mixed training design, prevalence across samples
  as a postprocessing step; and a latent bug the parity check found.
- **Validation**: 325 unit tests (16 new), 131 e2e tests (the calibrated-calls test rewritten), 24 trainer tests, 50
  mini-db tests (the pipeline test with its parity check); a build-and-train pipeline on the benchmark world (WSL
  `~/bench071/world`, 765 species, 135 unknown) with the new protal and defaults, and the feature groups ablated on its
  tables. Scripts here: [`build_fp.sh`](build_fp.sh), [`e2e_fp.sh`](e2e_fp.sh), [`bench_eval.sh`](bench_eval.sh);
  outputs [`bench_summary.txt`](bench_summary.txt), [`ablation_output.txt`](ablation_output.txt).

## The user's ideas, against the code

1. **"Is protal checking that the mate aligns too, where the gene leaves room for it?"** Partly, and not where it
   counts. The aligner aligns the mates on their own, pairs them when they land on one taxon and gene in opposite
   orientation (no position or insert check), rescues a missing mate on the guide's gene or its neighbours
   (`MateGuidance.h`), and gives split mates a consensus taxon. The profiler counts a pair with one mate aligned as
   one fragment and never looks at the mate's state. `linked_share` (both mates with a kept record on the taxon)
   is in the dump but excluded from the model since the tuning world. On the r226 v5 tables it carries the signal
   the user expects, weakly: 79% of true single-read species have both mates on the taxon, 59% of the false ones;
   as a feature it adds 0.001 of test F1 ([`linked_share.py`](../2026-10-03-false-positive-anatomy/linked_share.py),
   [`run_linked.sh`](../2026-10-03-false-positive-anatomy/run_linked.sh)). The sharper version, the user's, is now
   `mate_lost_share`: of the fragments whose mate was expected on the taxon (both mates with a record on it, or one
   kept record with room for the fragment inside its gene, the room judged against the length within which 95% of the
   sample's fragments with both mates on one gene lie), the share whose mate has no record on the taxon
   (`MateRoom`, `FragmentSpan`, `RecordEvidenceCollector::MateSpan`, `Taxon::MateLostShare`).
2. **"Does protal model that fast genes are hit at lower identity than conserved ones?"** No: the four `conserved_*`
   features compare depths and record counts by gene kind, and `identity` and `excess_median` mix the genes. Now
   `excess_scaled_median` (each read's divergence beyond its qualities over its gene's conservation factor, the median:
   the genome's divergence in the species definition's units, which the marker genes compress two- to threefold) and
   `excess_conserved_fast_ratio` (log2 of the divergence beyond errors on the conserved genes over the fast ones) do
   this. The same thing as suggestion 4 of the anatomy report ("ANI by conservation"); implemented once.
3. **"Can mismatches be split by codon position?"** Yes, cheaply. The database's genes are GTDB's marker gene calls
   copied unchanged (Prodigal CDS, coding strand, in frame), so reference position p (0-based) is codon position p mod
   3; the CIGAR has `X` for every mismatch. `third_position_share` walks each record's CIGAR from its position
   (`MismatchesByCodonPosition`): sequencing errors fall on the three positions alike, a strain's differences mostly
   on the third, a relative's less so (within species pN/pS exceeds between species dN/dS, so the direction may be
   either; the forest learns it). Cost: one more pass over each CIGAR; nothing for the forest.
4. **"Exclude genes at >99% identity in other genera from the database?"** Flagging beats deleting. `--build` now
   sketches every species' copy of each gene (128 bottom hashes of 12-mers) and compares each with the copies sharing
   one of its 8 smallest hashes; a copy within 0.02 of a copy of another genus, family, order, class, phylum or
   domain (`--suspect_copy_distance`), whose nearest congener's copy is 0.02 farther or missing, is suspect:
   contamination or a transfer (`GeneIncongruence.h`). A copy a whole clade shares nearly unchanged is not: its
   congeners are as near. The suspect copies go into the database (`suspect_copies.tsv`) and a run leaves records on
   them out as if the reads had not aligned (`--keep_suspect_copies` counts them), so the species gets no credit
   while the read is still absorbed; deleting both copies would send the reads to the next best copy, and deleting
   only one needs to know which is the contaminant, which the congeners tell only when there are some. Every near
   pair across genera (within 0.05) is reported in `gene_incongruence.tsv` with the verdicts. The likely loss: the
   flagged copies' share of all copies, which the build log states (`Suspect copies: N of M gene copies ...`); at
   r226 the within-genus report found 0.5% of species with an identical congener copy of a conserved gene, so the
   cross-genus share should be well below 1%. On the benchmark world (synthetic genomes, no contamination) the scan
   flags 7 of 83,978 copies (0.01%; its first rule 232, see below). The r226 number comes with the next build.

## What changed

### protal

- `sample_log_fragments` (`SampleEvidence`, `ApplySampleContext`): log10 of the sample's fragments over all its taxa,
  what `DepthKnobAt` reads; written for every taxon.
- `RecordEvidence` (every best record, before the filters): `scaled_excess` (excess over the gene's factor), per kind
  of gene the differences, aligned bases and expected errors in ppm (integers: chunks add up alike), mismatches by
  codon position, `mates_linked` and `mate_room`; the collector's `m_spans` (fragments with both mates on one gene)
  gives `MateSpan()`, the 95th percentile, or 500 with fewer than 50. `Taxon::SetRecordEvidence(records, mate_span)`
  reduces them to `ExcessScaledMedian`, `ExcessConservedFastRatio` (0 with fewer than 200 aligned bases on either
  kind), `ThirdPositionShare` (1/3 without mismatches) and `MateLostShare`; `TaxonFeatures` emits the five.
- `PrepareMAPQ`: records on suspect copies (`MicrobialProfile::SuspectCopy`, `GenomeLoader::IsSuspectCopy`) are
  dropped before anything counts them, and counted (`SuspectRecords`, a line per sample); then the mates are judged.
- `--build`: `WriteSuspectCopies` after the congener report (`gene_incongruence::Scan`, OpenMP over genes, results
  sorted: the same for any thread count); `suspect_copies.tsv` bundled, `gene_incongruence.tsv` beside the database;
  `--suspect_copy_distance` (0.02; 0: no scan). `LoadSuspectCopies` at start; `--keep_suspect_copies`.
- `gene_conservation::MixKmer`, `BottomSketch`, `SketchDistance`: the sketches `SampleContext.h` had, now shared.
- `--fdr`: off by default; a model's calibrated calls are used only with `--fdr F` (the load log says
  "not used (--fdr F would use them)"). `--singleton_congener` 0 by default (`context::kSingletonCongener` 0).

### Trainer and build script

- `model_features.py`: feature groups joined by `+` (`normalized`, `adjacency`, `relatives`, `distance`, `depth`,
  `divergence`; `feature_set_columns`, `feature_set_name` for argparse, `has_sample_depth`); `DEFAULT_FEATURE_SET`
  `normalized+adjacency+distance+depth+divergence`. Tables of older protals: `normalized+adjacency+distance` or
  `normalized+adjacency`.
- `random_forest_cmdline.py`: no knob curve when the depth is a feature (the report says why, the metrics say
  `skipped`); a knob point keeps `--knob` unless its best knob gains `DEPTH_KNOB_MIN_GAIN` = 0.002 on its window
  (suggestion 1 of [the v4 report](../2026-10-03-r226-v4-v2-rerun/README.md)); `--singleton-congener` 0.
- `build_gtdb_database.py`: `suspect_copies` in `build_metadata.tsv`, `gene_incongruence.tsv` in `model_logs/`,
  `model_<type>_depth_knobs` says when no curve was fitted.
- `scripts/recurrent_calls.py`: a run's `*.profile.log` files read for species called on a read or two in several
  samples always beside the same more abundant relative (its companion and the rank shared): references' artefacts
  the database did not flag.

### Tests and documentation

- `tests/test_FalsePositiveFeatures.cpp` (10 tests): codon positions, mate room and span, the error probability,
  the features on a profile, lost mates through the SAM path on 1 and 3 threads, suspect copies dropped on 1 and 2
  threads, lineages from the taxonomy, the scan on three genera on 1 and 3 threads, the table's round trip.
  `test_SampleContext.cpp` sets the singleton rule on where it tests it.
- `scripts/test_model_pmml.py`: the groups, the default, no curve with the depth feature, the gain rule, the singleton
  default; `tests/e2e/test_protal_e2e.py`: `FalseCallsTest.test_calls_only_with_fdr`;
  `scripts/mini_db/test_mini_db.py`: the default set and the suspect-copies metadata.
- `docs/running.md` (`--fdr`, `--singleton_congener`, `--keep_suspect_copies`, the run-level check),
  `docs/model-training.md` (the two feature groups, `--features`, the knob curve's gain rule and its absence, calls
  only with `--fdr`, the singleton rule off), `docs/building-a-database.md` (the scan, the files, `--features`,
  `--call-mode`), `docs/database-files.md` (`suspect_copies.tsv`).

## Validation on the benchmark world

Pipeline `new` ([`e2e_fp.sh`](e2e_fp.sh): `build_gtdb_database.py` with the branch's protal, scripts and defaults;
pe and se; 8 samples per point at 1,000 to 1M read pairs, 4 test samples per point at 500 to 1M; seed 1):
[`bench_summary.txt`](bench_summary.txt). Then the feature groups ablated on its tables, every model with
`--depth-knobs` (fitted only without the depth feature), scored on its test set ([`bench_eval.sh`](bench_eval.sh),
[`ablation_output.txt`](ablation_output.txt)):

Test-set F1 at knob 0.5 (what protal calls with for a model without a curve) / at the model's knob curve where one was
fitted / at the best test threshold; false positives and misses at 0.5; `nad` = `normalized+adjacency+distance`:

| read type | features | species held out | test 0.5 / curve / best | FP / FN at 0.5 | FP / FN at the curve | log loss |
|---|---|---:|---|---|---|---:|
| pe | `nad` | 0.9893 | 0.9795 / 0.9748 / 0.9814 | 26 / 119 | 18 / 159 | 0.0676 |
| pe | `nad` + divergence | 0.9903 | **0.9807** / 0.9818 / 0.9827 | 25 / 112 | 28 / 101 | 0.0620 |
| pe | `nad` + depth | 0.9895 | 0.9748 / — / 0.9813 | 21 / 156 | — | 0.0719 |
| pe | `nad` + depth + divergence (default) | 0.9913 | 0.9775 / — / **0.9835** | 20 / 139 | — | **0.0637** |
| se | `nad` | 0.9866 | 0.9732 / 0.9707 / 0.9765 | 37 / 148 | 26 / 175 | 0.0752 |
| se | `nad` + divergence | 0.9881 | **0.9734** / 0.9752 / 0.9786 | 35 / 148 | 46 / 126 | 0.0705 |
| se | `nad` + depth | 0.9888 | 0.9664 / — / 0.9781 | 25 / 204 | — | 0.0807 |
| se | `nad` + depth + divergence (default) | 0.9901 | 0.9676 / — / **0.9791** | 23 / 198 | — | 0.0779 |

(The pipeline's own models, `normalized+adjacency+distance+depth+divergence` without a curve, score 0.9775 pe and
0.9676 se on its test set, [`bench_summary.txt`](bench_summary.txt). About ±0.003 of F1 is noise here; the first
pipeline, with the scan's first rule, gave 0.9769 and 0.9685, [`bench_summary_first_rule.txt`](bench_summary_first_rule.txt),
and its ablation the same picture, [`ablation_first_rule_output.txt`](ablation_first_rule_output.txt).)

- **The divergence features gain a little here and rank better**: +0.001 (pe) and +0.000 (se) at knob 0.5, the log
  loss down 8% and 6%, misses down (pe 119 → 112) with the false positives unchanged. Positive, within noise on this
  world; `excess_scaled_median` is the forest's top feature on it (importance 0.12 for pe and se, ahead of the
  unique-k-mer rates and `identity`), `mate_lost_share` and `third_position_share` are not in the top eight.
- **The depth feature loses here**: −0.005 (pe) and −0.007 (se) at knob 0.5, through misses in the deep test samples
  (pe 119 → 156, 37 of them single fragments), while the false positives, few to begin with on this world (26 of
  9,068 test taxa), fall by 5. Its ranking is the best of the sets (best-threshold F1 0.9835 against 0.9814, the lowest
  log loss with the divergence features), so the forest learnt the depth right and the threshold wrong: the training
  design (abundances of σ 1.3) has fewer rare species in deep samples than the test design (σ 2.0), and the forest
  with the depth feature learnt that prior, which the knob curve, fitted on the same training samples, had learnt
  too (it also loses here: pe 0.9795 → 0.9748) but applied more coarsely. On the r226 v5 tables, with the same two
  designs, the depth feature gained +0.004 (pe) and +0.009 (se) at knob 0.5 and halved the false positives
  ([the anatomy report](../2026-10-03-false-positive-anatomy/README.md)): there the database has 8,000 species and
  the absent single-read taxa of a deep sample outnumber the present ones 60 to 1 at identity ≥ 0.99, here 765
  species and few false positives, so the same prior pays there and costs here. **The default follows the r226
  result, as the user decided, and the next r226 build decides**; if it loses there too, `--features
  normalized+adjacency+distance+divergence`, or a training design whose abundances match the test's (σ 2.0, or rare
  species added to the deep points), which would let the forest learn the right prior for both.
- **Suspect copies**: the scan's first rule (a copy of a singleton genus suspect whenever another genus's copy is
  within 0.02) flagged 232 of 83,978 copies (0.28%) of 125 species on this world, 227 of them in singleton genera
  (165 of its 276 genera), on its slowest genes: a slow gene a young family shares, not contamination. The rule now
  asks such a copy to lie inside the other genus's cluster by the margin, and flags 7 copies (0.01%) of 5 species,
  two genera of one family whose gene 1 copies are at 0.011 while their congeners' are at 0.035 (the synthetic
  taxonomy disagreeing with that gene). 69,429 near pairs across genera within 0.05 are reported either way
  (`gene_incongruence.tsv`). The scan took 7.7 s for 84k copies on 6 threads; at r226 (18M copies) expect minutes.
  This world has no contamination, so the scan's benefit is not measurable here; the first run's 232 dropped
  copies changed the pipeline's test F1 by −0.0006 (pe) and +0.0009 (se) against the second's 7, within noise: the
  cost of flagging 0.3% of copies is small, the gain at r226 is the open question.

## What is open

- **GTDB scale decides.** The depth feature was validated on the r226 v5 tables offline (+0.004 pe, +0.009 se test
  F1 at a fixed knob; the anatomy report); the divergence features and the suspect-copies scan were only testable on
  the benchmark world here, whose genomes are synthetic (no contamination, synthetic gene rates). The next r226 build
  with these defaults tells: compare its `model_logs/summary.txt` with v5's (pe 0.9615 / se 0.9570 test F1 at knob
  0.5), and read `build_metadata.tsv`'s `suspect_copies` and `gene_incongruence.tsv` for what the scan flagged.
  A depth point beyond 10M read pairs is worth adding to the design (the forest cannot extrapolate past it).
- Not implemented from the anatomy report: the novel-species report below the species' radius
  (`sp_clusters_r226.tsv`; another session found genus-level reporting of called taxa does not pay,
  [genus fallback](../2026-10-03-genus-fallback/README.md)) and species complexes for pairs within 97% ANI.
- Dropped: genus size as a prior, identity rules against singletons, the singleton rule by default, more skew
  features.

## Round two: the strategies of the F1 list, items 2, 3, 4, 5 and 7

Implemented after the misses were dissected (328 misses of the r226 v5 model: 276 strains, 208 with one or two
fragments, a median of 41% of their records below MAPQ 10 while the EM gave them nearly all; the list is in the
memory notes and in the reply of 2026-10-03). Same branch; the tests of round one plus
[`tests/test_ReadEvidence.cpp`](../../../tests/test_ReadEvidence.cpp) (5 tests).

### 2. Fragments before the filters, and the EM's fragments

`fragments` counts the records the MAPQ and length filters keep; a divergent strain's reads tie with a congener's
reference and fall to MAPQ 0 before they count. `fragments_all` counts the taxon's reads with a best record before the
filters (a pair or a long read once, `RecordEvidenceCollector::NoteFragment`), and `em_fragments` is the share the
abundance-weighted assignment leaves to the taxon times that (`em_own_share × fragments_all`): the fragments the strain
would have had. Spill-over keeps a low EM share and gains little.

### 3. What GTDB knows of a species before any read

The converter writes `species_priors.tsv` (`gtdb_to_protal_db.py`, `SpeciesPriors.h`): per species the marker genes
found in the representative and how many of them twice (the converter kept a genome's first copy and dropped the rest
without counting; a second copy of a single-copy marker is CheckM's contamination signature), the representative's
CheckM completeness and contamination from the metadata (`checkm2_*`, else `checkm_*`), and from GTDB's
`auxillary_files/sp_clusters_r226.tsv` (now fetched by `download_gtdb.py` when the release lists it) the cluster's ANI
circumscription radius, mean and minimum intra-species ANI and size. `--build` bundles it, a run loads it
(`LoadSpeciesPriors`), and `TaxonFeatures` emits `rep_duplicate_share`, `rep_completeness`, `rep_contamination`,
`cluster_ani_radius`, `cluster_mean_ani`, `cluster_min_ani`, `cluster_genomes_log10` (-1 unknown). On the synthetic
worlds every value is unknown (the simulated releases have no duplicated markers, 100/0 quality and no clusters file),
so the group does nothing there; the r226 build is where it can act, and `build_metadata.tsv`'s converter log line says
how many species have values.

### 4. Reads that seeded on a taxon but did not align to it

`SimpleAlignmentHandler` records the taxa of the anchors it tries (`Attempted()`); the read paths subtract the taxa of
the valid alignments (`FailedCandidates`, after the mate guidance for pairs; `LongReadAligner::FailedTaxa` for long
reads) and hand the rest to the output handlers, which write them as `ZF:Z:<taxid>,...` on the read's first record,
or on a minimal unmapped record (flag 4, `*` for the sequence) when the read aligned nowhere; reads without failed
candidates get nothing, as before. The reader parses the tag on mapped records and counts the unmapped records' tags
per taxon before skipping them (`SamReader::FailedCandidates`); the profiler adds both into
`RecordEvidence::failed_candidates`, and `failed_candidate_rate` is the share of the reads that seeded on the taxon
and failed, over those plus the reads with a record. A SAM of an older protal gives 0. The e2e tests that read SAM
records skip the unmapped ones, and the count of skipped unmapped records in the "edited SAM" test allows for protal's
own.

### 5. The training design

`simulate_metagenomes --pln_sigma 1.3,2.0`: several sigmas are given to a design point's samples in turn
(`ProfileDesignOptions::pln_sigmas`, `SigmaForSample`), `collect_training_data.py --abundance lognormal:1.3,2.0`
passes them through, and `build_gtdb_database.py` defaults to it (half the training samples with the former sigma
1.3, half with the test set's 2.0), holds out 30% of the species instead of 20%, and adds a design point of 30M read
pairs (`30000000:1`: three samples) so that a model with the sample's depth as a feature has seen the depth real
metagenomes reach. The test design is unchanged (σ 2.0, 500 to 5M pairs).

### 7. Prevalence across the samples of a run, as a postprocessing step

`scripts/prevalence_calls.py OUTPUT_DIR`: for every taxon and sample, the odds of the model's probability times the
odds ratio of the taxon's prevalence in the run's other samples (mean probability, absence counting 0, smoothed by two
pseudo-samples at the base rate) to the base rate (the model's training prior if given, else the run's mean), capped
at a factor 4, then calls at the knob; writes `prevalence_calls.tsv` and, with `--profiles`, adjusted `.profile` files.
On the 42 paired-end test samples of the benchmark pipeline, drawn independently, the full update (`--direction both`)
removed 24% of the calls (564 of 2,388; 585 without the cap), nearly all of them true, and gained none, as a prevalence
prior must where samples share no species. So the default is `--direction up`, which only boosts prevalent species
and changed nothing there; the full update is for the samples of one study, and needs a study with a truth to be
validated, which the simulator cannot provide until its samples share a species pool. (The probabilities in these
`.profile.log` files are the collection model's, which the pipeline profiles with; the mechanics, not the F1, are what
this tested.)

### A latent bug the parity check found

The pipeline's parity check crashed protal on one deep sample: `em_fragments` was 9×10⁻³¹¹ for a taxon beside a
hugely abundant congener (the EM leaves it next to nothing), the dump wrote the subnormal number, the trainer (whose
forests take any float) used it, and cPMML, which reads a feature back with `stod`, took the underflow as a missing
value and threw. `em_own_share` has carried such values since the EM exists; the relatives set, opt-in, would have
crashed the same way. `FeatureString` now writes a subnormal value as 0, in the dump and at scoring alike, so the
trainer and protal agree (test `FeatureStringsFlushSubnormalValuesToZero`).

### Validation

325 unit tests (6 new), 24 trainer tests, 131 e2e tests (those reading SAM records skip the unmapped ones), 50
mini-db tests. A third benchmark-world pipeline with the round-two binary and defaults ([`bench3.sh`](bench3.sh):
30% of species held out, abundances σ 1.3 and 2.0 in turn, the same 1,000 to 1M read-pair points and test design as
before; its parity check passes), [`bench_summary3.txt`](bench_summary3.txt), and the groups ablated on its tables
([`ablation3_output.txt`](ablation3_output.txt)). Its numbers are not comparable with the second pipeline's: with 30%
of the species held out the test samples score 8,504 taxa against 9,068, with more of their species missing from the
database. Test-set F1 at knob 0.5 / at the model's knob curve where one was fitted / at the best test threshold:

| read type | features | species held out | test 0.5 / curve / best | FP / FN at 0.5 | AP | log loss |
|---|---|---:|---|---|---:|---:|
| pe | `nad` | 0.9809 | 0.9733 / 0.9697 / 0.9755 | 36 / 132 | 0.9928 | 0.0856 |
| pe | `nad` + divergence | 0.9829 | **0.9750** / 0.9744 / 0.9765 | 29 / 128 | 0.9930 | 0.0831 |
| pe | `nad` + divergence + unfiltered | 0.9827 | 0.9744 / 0.9749 / 0.9765 | 33 / 128 | 0.9956 | 0.0647 |
| pe | `nad` + divergence + unfiltered + priors | 0.9841 | 0.9741 / 0.9741 / 0.9765 | 38 / 125 | **0.9958** | **0.0629** |
| pe | + depth (default) | 0.9843 | 0.9734 / — / 0.9755 | 31 / 136 | 0.9954 | 0.0684 |
| se | `nad` | 0.9738 | 0.9628 / 0.9576 / 0.9670 | 52 / 176 | 0.9916 | 0.0977 |
| se | `nad` + divergence | 0.9757 | 0.9663 / 0.9624 / 0.9705 | 37 / 169 | 0.9923 | 0.0961 |
| se | `nad` + divergence + unfiltered | 0.9782 | 0.9659 / 0.9662 / 0.9703 | 41 / 168 | 0.9939 | 0.0836 |
| se | `nad` + divergence + unfiltered + priors | 0.9788 | **0.9670** / 0.9662 / 0.9720 | 39 / 163 | **0.9944** | **0.0799** |
| se | + depth (default) | 0.9796 | 0.9653 / — / 0.9714 | 37 / 175 | 0.9947 | 0.0813 |

- **The unfiltered group ranks far better and calls the same.** The log loss falls by 22% (pe) and 13% (se) and the
  average precision rises by 0.0026 and 0.0016 on top of the divergence features, while the F1 at knob 0.5 stays
  within noise (−0.0006, −0.0004) and the misses at one fragment fall from 53 to 50 (pe). The information the MAPQ
  filter had thrown away is now in the scores; on this world the strains it should rescue are too few, or too far,
  for the threshold to move. The r226 misses, 159 of 328 with half or more of their records dropped by the filter, are
  the population it was built for.
- **The priors do nothing here**, as expected: the synthetic release has no duplicated markers, 100/0 quality for every
  genome and no clusters file, so every value is unknown. The changes against the row above are noise.
- **The depth feature's cost shrank to noise with the mixed training design**: −0.0007 (pe) and −0.0017 (se) at knob
  0.5 against the set without it, where the second pipeline, trained on σ 1.3 alone, lost 0.005 and 0.007. The
  training prior now spans the test's abundance distribution, which is what the design change was for; the r226
  gain stays to be shown there.
- **The failed candidates are plentiful**: a 200,000-pair sample of this world has about 38,000 unmapped records
  (19% of its fragments seeded on a taxon and aligned nowhere, most of them reads of the 30% of species the database
  lacks), each a 70-byte line; the SAM grows by a few percent compressed, and the reader skips them in the same pass.
- Nothing changed in the pipeline's own calls through the parity check once the subnormal values were fixed.

## Website

The website lists protal's options: `--fdr` is now off by default, `--singleton_congener` defaults to 0, and
`--keep_suspect_copies` and `--suspect_copy_distance` are new; databases built from now on carry
`suspect_copies.tsv`, and models trained from now on use the new default feature set (older databases' models keep
working: protal provides every feature a model asks for). The website is not updated.
