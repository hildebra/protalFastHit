# Against false positives: the sample's depth and the divergence features, suspect gene copies, defaults, a run-level check

- **Date**: 2026-10-03.
- **Branch**: `fp-anatomy` (worktree `../protal-fp`), from `audit-fixes` at `4cae2f4`.
- **What**: the user's four ideas evaluated against the code and the r226 v5 tables, and the agreed items of
  [the anatomy report](../2026-10-03-false-positive-anatomy/README.md) implemented: the sample's depth as a feature
  (default, no knob curve with it), four divergence features (by gene conservation, by codon position, lost mates),
  suspect gene copies found at build time and left out of the evidence, `--fdr` and the singleton rule off by default,
  knob points that must gain, `scripts/recurrent_calls.py`.
- **Validation**: 318 unit tests (10 new), 131 e2e tests (the calibrated-calls test rewritten), 24 trainer tests, 49
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

## Website

The website lists protal's options: `--fdr` is now off by default, `--singleton_congener` defaults to 0, and
`--keep_suspect_copies` and `--suspect_copy_distance` are new; databases built from now on carry
`suspect_copies.tsv`, and models trained from now on use the new default feature set (older databases' models keep
working: protal provides every feature a model asks for). The website is not updated.
