# The four model features and the depth knobs on the v0.7.1 benchmark

- **Date**: 2026-10-01.
- **Code**: 0.7.2-dev = `25d457e` (branch `audit-fixes`; still Version 0.7.1): 0.7.1 (`1c11a00`) plus the parallel
  profiling of `f8b1b6b`-`f5ee645` (byte-identical profiles, `../2026-10-01-multithreading-audit`), the four model
  features of `e680d91` (`excess_median`, `excess_high_share`, `conserved_fast_depth_ratio`, `conserved_hit_share`)
  and the depth knobs of `25d457e` (the trainer's `--depth-knobs`, which `build_gtdb_database.py` passed for pb and
  ont). Exported with `git archive`, built in Release in WSL.
- **Data**: the v0.7.1 benchmark (`../2026-10-01-v071-benchmark`, WSL `~/bench071`): its world, and 0.7.2-dev's own
  pipeline on its release (`V072dev`), the 0.7.1 pipeline's design and seed with 0.7.1's simulator, so the same
  training and test samples as `V071`; the same species held out of the training database (checked). Then the
  benchmark's samples (26 paired-end incl. 2 of 5M pairs, 24 single-end, 8 PacBio, 8 Nanopore) against its finished
  database (full) and its training database (missing: 290 species left out).
- **Commands**: `scripts/run.sh` (build, pipeline, runs `v072.*`, and `v072k.*`: the long-read samples again with
  `--knob 0.5`, i.e. the features without the depth knobs), `scripts/compare.py` (scores with the benchmark's
  `score.py`; paired differences per sample with a 95% bootstrap interval). Results: `results/summary.md`,
  `results/runs.tsv`, `results/depth_knobs.txt` (the trainer's sections on them).
- **Machine**: WSL2, 6 threads, shared with another session's builds.

## The four features: a gain on paired-end and Nanopore reads

Paired difference per sample, 0.7.2-dev (or, for long reads, 0.7.2-dev at `--knob 0.5`: the features alone) less
0.7.1, mean (95% interval):

| reads | full database: F1 | FP per sample | Bray-Curtis | missing database: F1 | FP per sample | Bray-Curtis |
|---|---|---|---|---|---|---|
| pe | **+0.0031** (+0.0009, +0.0054) | **−0.62** | −0.0016 | +0.0057 (+0.0000, +0.0136) | **−0.50** | −0.0032 |
| se | +0.0005 (−0.0030, +0.0037) | +0.12 | −0.0009 | +0.0036 (−0.0030, +0.0103) | 0.00 | −0.0035 |
| pb, features alone | +0.0020 (−0.0013, +0.0060) | −0.38 | 0.0000 | +0.0008 (−0.0028, +0.0043) | 0.00 | −0.0003 |
| ont, features alone | **+0.0071** (+0.0019, +0.0128) | **−1.12** | −0.0018 | +0.0057 (−0.0039, +0.0158) | −0.38 | −0.0013 |

As the offline refit predicted (+0.003 pe, +0.007 ONT, `../2026-10-01-f1-opportunities`), mostly through fewer false
positives: paired-end 2.12 to 1.50 per sample with the full database, Nanopore 1.88 to 0.75. No read type loses.
The pipeline's own estimates agree: F1 of species held out pe 0.9791 to 0.9816, se 0.9747 to 0.9790, ONT 0.9750 to
0.9834 (PacBio 0.9729 to 0.9722); on its independent test set pe 0.9682 to 0.9726, PacBio 0.9658 to 0.9723.

## The depth knobs: they do not carry over; PacBio loses

| reads | the knobs (bin: knob) | full: F1, knobs less `--knob 0.5` | FP per sample | missing: F1 | FP per sample |
|---|---|---|---|---|---|
| pb | 2: 0.13, 3: 0.63, 4: 0.74 | **−0.0074** (−0.0177, −0.0002) | +1.00 | **−0.0163** (−0.0276, −0.0086) | +2.38 |
| ont | 2: 0.18, 3: 0.49, 4: 0.55 | +0.0005 (+0.0000, +0.0015) | −0.25 | −0.0057 (−0.0164, +0.0020) | +1.12 |

On the pipeline's own test set they looked fine (ONT F1 0.9547 to 0.9660, PacBio 0.9723 to 0.9727) and on species
held out they gained (chosen there, so optimistic). On the benchmark's samples they cost PacBio up to 0.016, barely
help Nanopore with the full database (+0.0005) and cost it 0.006 with species missing (not significant). Why (`results/depth_knobs.txt` and the runs' logs, `Sample ...: N fragments, knob K`):

- **Few samples per bin, and the shallow ones have few absent taxa.** Bin 2 (fewer than 1,000 fragments) of the
  PacBio model rests on 5 training samples with 288 taxa, 256 of them present: with 32 absent taxa, the threshold
  with the highest F1 drops to 0.13. The benchmark's 3 Mb samples have more reads per sample and more absent taxa
  with reads, and at 0.13 they come through.
- **Hard bin edges.** The 3 Mb samples have 400 to 1,450 fragments, on both sides of the 1,000 edge: PacBio samples
  with 770, 920 and 998 fragments were called at 0.13, one with 1,010 at 0.63.
- The offline estimate of +0.015 (ONT) and +0.007 (PacBio) in `../2026-10-01-f1-opportunities` was on the 0.7.1
  pipeline's own test set, which has the training design's depths; it overstated what knobs by depth give on other
  samples.

**Consequence.** `build_gtdb_database.py` no longer passes `--depth-knobs` by default (`--depth-knob-read-types`,
empty by default); the trainer option and protal's use of knobs in a model stay. Before they are used again they
would need knobs that change smoothly with depth (interpolated between the training depths rather than per bin),
shrunk towards `--knob` by how many samples and absent taxa a depth has, and training designs with several samples
per depth that also hold absent taxa at low depth.

## Not covered

Real samples and a GTDB-scale database; 8 long-read samples per database here, so the long-read intervals are
wide. The conservation features rest on the simulator's gene rates, checked read by read in
`../2026-10-01-conservation-pattern`.
