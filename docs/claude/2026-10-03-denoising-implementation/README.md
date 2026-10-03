# Amplicon-denoising ideas in protal: congener-rich training, relatives features, read EM, the singleton rule and calls at a target share of false calls

- **Date**: 2026-10-03.
- **Branch**: `denoise` (worktree `../protal-denoise`), from `audit-fixes` at `6ca86ee`; uncommitted when written.
- **What**: suggestions 1-5 of [the amplicon-denoising report](../2026-10-02-amplicon-denoising/README.md), as the user
  asked, and the user's refinement while they were being written: model a low-abundance congener as part of the
  abundant one with the error rates its reads show, and from how its reads lie on the genes.
- **Validation**: unit, e2e and mini-database tests; two build-and-train pipelines on the benchmark world (WSL
  `~/bench071/world`: 765 species of a GTDB-like release, 135 more in no database, 3 genomes each, congeneric species
  3-12% apart), one with the new defaults and one with the old design, cross-scored (below). Scripts and outputs here.

## Summary

All five are in protal, the simulator and the scripts, tested; on the benchmark world only part of them pays:

- **Congener groups** in the training and test samples (on by default): within noise for the old features; they make
  the test sets show what congeners do.
- **Relatives features and the read EM**: the model ranks better (AP, log loss, best-threshold F1), but test F1 at
  protal's knobs falls by up to 0.005 and minor congeners are missed 1.6-2.3 times as often. Their four distance
  features alone are the best paired-end set (+0.001 to +0.006). Opt-in (`--features normalized+adjacency+relatives` or
  `+distance`), to be judged on an r226 build, where false positives dominate.
- **Singleton rule**: as first specified it removed true calls of minor congeners. With the user's condition on the
  read (vetoed only if the EM gives the read mostly to the abundant congener, or it is more than 5% from the reference)
  it removes only false calls. On by default.
- **Calls at a target share of false calls**: 0.0005-0.004 below the knob curve. protal uses them when a model has them;
  the build trains them only with `--call-mode fdr`.

## What changed

### 1. Congener groups in the training and test samples

- `simulate_metagenomes --congener_groups SHARE:MIN-MAX` (`ProfileDesignOptions::congener_share/min/max`,
  `parse_congener_groups`, `CommunityProfileDesigner::design_profile`): after the `--genus` and `--taxon` demands, genera
  in random order among those with MIN species or more, each with MIN to MAX of its species, until about SHARE of the
  sample's species are in such groups. Without it nothing is drawn there, so older commands make the same samples.
- `collect_training_data.py --congeners` takes `SHARE:MIN-MAX` (passed on as `--congener_groups`), or N as before (one
  genus of N per design point), or 0.
- `build_gtdb_database.py --congeners` defaults to `0.25:2-5` (was 0), for the training data and the test set;
  `build_metadata.tsv` records it in `classifier_training_design`.

### 2. Features from the other taxa of a sample

`MicrobialProfile::ApplySampleContext` (`src/Profiling/SampleContext.h`, `Profiler.h`), run at the end of
`ApplyRecordEvidence`, once a sample's reads are in, gives every taxon a `SampleEvidence`; `TaxonFeatures` writes:

| feature | |
|---|---|
| `genus_skew`, `family_skew`, `genus_share`, `genus_spill` | the rank features of the offline experiment: against the genus's most abundant other species, the family's of another genus, the genus's share, and seen over expected spill-over at 10⁻³ per congener fragment and 10⁻⁴ per fragment of the family's other genera |
| `relative_skew`, `relative_distance`, `relative_spill` | against the congeners with more fragments (at most 4) by the distance d of the references: median Mash distance (k = 12, bottom-64 sketches per gene, at least 10 shared genes, else 1) and a spill rate 0.01 × 10^(−d/0.05) per fragment; the likeliest source's skew and distance, and seen over expected summed spill-over |
| `relative_close_share` | the user's "distribution in the low congener": of the taxon's fragments on the genes it shares with the likeliest source, the share on the half of them where the two references are most alike, over that half's share of their length. ~1 for a species' own reads, up to 2 for spilled-over ones |
| `em_own_share`, `em_kept_own_share` | suggestion 3 with the user's "observed error rates": an EM over every best record's alternatives (ZA, within 5 edits), each read shared among the taxa it fits in proportion to their records times (θ / (1 − θ))^(edits more), θ the taxon's own differences per aligned base in the sample (0.002-0.2); the share of the taxon's records (all, and those kept by the MAPQ and length filters) it keeps |
| `genus_top_fragments` | the singleton rule's input (not in the feature sets) |

The distances are computed only for the pairs of congeners that share a sample, from the references' genes, and kept
for the run (`CongenerDistances`, capped); they need no database file, so existing databases work. Everything is
summed and visited in taxid and key order, so a sample's features are the same on any number of threads (the build
test's parity check passes). `model_features.RELATIVE_FEATURES` (`--features normalized+adjacency+relatives`) and its
four distance features alone, `DISTANCE_FEATURES` (`normalized+adjacency+distance`); neither is the default set (below).

### 4. The singleton rule

A taxon of one fragment beside a congener of 100 fragments or more is never reported if its read looks like the
congener's: its EM share below 0.5, or its identity below 0.95 (`Taxon::Vetoed`,
`TaxonFilterForest::Calls`, used by every output: profiles, gene logs, the dump's `prediction`, statistics, strain MSAs,
`unreported_species.tsv`). `protal --singleton_congener N` (0: none); `random_forest_cmdline.py --singleton-congener`
counts it in every call (from `fragments`, `genus_top_fragments`, `em_own_share` and `identity`). The condition on the
read came from the validation below: without it the rule removed true calls of minor congeners.

### 5. Calls at a target share of false calls

- Trainer `--fdr-calls`: an isotonic fit of presence on the scores of species held out (64 quantile points and 0, 1), the
  share of present rows (the prior), and the target (0.005-0.5) with the highest F1 on species held out; in the model's
  header as `protal_calibration`, `protal_prior`, `protal_fdr`; reported on the test set next to the knob curve and 0.5.
- protal (`context::FalseCallKnob`, `SampleAdjusted`): each sample's scores made probabilities by the curve, adjusted to
  the sample's share of present candidates (Saerens et al. 2002, with 20 pseudo-candidates at the prior: without them a
  sample of three likely taxa reached a share of 1 and called everything even at a target of 10⁻⁶, which the e2e test
  found), and the highest-scoring taxa called while the mean of their 1 − probability is at most the target. Used by
  default ahead of the knob curve when a model has them; `--fdr F` another target, `--fdr 0` off; not with `--knob`.
- `build_gtdb_database.py --call-mode fdr` passes `--fdr-calls`; its default, `curve`, does not (below).

## Tests

- C++: 308 unit tests pass (12 new in `tests/test_SampleContext.cpp` and `tests/test_CommunityDesign.cpp`: sketches,
  distances, spill rate, edit ratio, EM shares, prior adjustment, the false-call knob, the header, a profile of two
  congeners and a species of another genus with the features, the veto and the close share; congener groups and their
  option).
- e2e: 131 pass (3 new: calls at the model's target, `--fdr` and `--knob`, malformed calibrations).
- `scripts/test_model_pmml.py`: 23 pass (the extension round trip, the trainer's prior adjustment, calls at a target
  checked against the rule sample by sample, the calibration, the singleton rule, `--fdr-calls` in the model).
- `scripts/mini_db/test_mini_db.py`: 49 pass, among them `GtdbBuildTest` (the whole pipeline with the new defaults, its
  parity check included).

## Validation on the benchmark world

### Setup

- **Pipelines** (`bench.sh`, the branch's protal and scripts, `build_gtdb_database.py`, pe and se, seed 1, 6 threads, 13.5
  minutes each):
  - **new**: congener groups 0.25:2-5;
  - **base**: `--congeners 0`, uniform draws as before.
  - Design: 8 samples per point at 1,000, 20,000, 200,000 and 1M read pairs (4); test sets of 4 samples at 500,
    10,000, 100,000 and 1M (2).
  - The same species held out (the same training database). Their summaries: `pipeline_*_summary.txt`.
- **Cross-scoring** (`eval2.sh`):
  - **Features:** `na` (`normalized+adjacency`) and `na+relatives`, plus the groups `na+em`, `na+distance` and `na+rank`
    (`ablation_features.patch`).
  - **Training and testing:** every model trained on each training table and scored on both test sets (uniform:
    "base"; with congener groups: "new").
  - **Trainer settings:** `--depth-knobs --fdr-calls`, 64 trees, 512 leaves, the refined singleton rule.
  - **Calls:** at 0.5, at the model's knob curve, and at its target share of false calls (`eval2_summary.py`, output
    `eval_summary.txt`).
  - **Minor congeners:** present species beside a present congener of ten times their fragments or more (331 and
    397-407 of them in the two test sets).
- The test sets have another design than the training data (abundance σ 2.0 against 1.3, 10-300 species, more
  strains). With ~3,500 present taxa each, about ±0.003 of F1 is noise.

### Results

Test F1, paired-end (pe) and single-end (se), on the uniform / congener-rich test set:

| model | trained on | at 0.5 | at the knob curve | at the target share of false calls | AP | best threshold |
|---|---|---|---|---|---|---|
| pe `na` | uniform | 0.9772 / 0.9784 | 0.9762 / 0.9781 | 0.9751 / 0.9765 | 0.9943 / 0.9951 | 0.9790 / 0.9811 |
| pe `na` | congener groups | 0.9757 / 0.9784 | 0.9760 / 0.9788 | 0.9719 / 0.9766 | 0.9945 / 0.9947 | 0.9792 / 0.9801 |
| pe `na+relatives` | congener groups | 0.9753 / 0.9754 | 0.9743 / 0.9742 | 0.9748 / 0.9757 | 0.9952 / 0.9964 | 0.9801 / 0.9812 |
| pe `na+distance` | congener groups | **0.9771 / 0.9791** | **0.9794 / 0.9805** | **0.9778 / 0.9794** | 0.9944 / 0.9962 | 0.9800 / 0.9816 |
| pe `na+em` | congener groups | 0.9770 / 0.9780 | 0.9773 / 0.9786 | 0.9722 / 0.9747 | 0.9948 / 0.9955 | 0.9814 / 0.9812 |
| pe `na+rank` | congener groups | 0.9764 / 0.9788 | 0.9749 / 0.9770 | 0.9785 / 0.9786 | 0.9947 / 0.9959 | 0.9801 / 0.9811 |
| se `na` | uniform | 0.9682 / 0.9740 | 0.9700 / 0.9736 | 0.9674 / 0.9714 | 0.9928 / 0.9951 | 0.9704 / 0.9743 |
| se `na` | congener groups | 0.9703 / 0.9728 | 0.9672 / 0.9708 | 0.9656 / 0.9703 | 0.9926 / 0.9948 | 0.9722 / 0.9772 |
| se `na+relatives` | congener groups | 0.9702 / 0.9690 | 0.9657 / 0.9665 | 0.9641 / 0.9659 | 0.9941 / 0.9957 | 0.9739 / 0.9774 |
| se `na+distance` | congener groups | 0.9708 / 0.9727 | 0.9679 / 0.9683 | 0.9677 / 0.9706 | 0.9935 / 0.9956 | 0.9733 / 0.9759 |

Minor congeners missed at 0.5 (uniform / congener-rich test set), trained on congener groups: pe `na` 20 / 19,
`na+relatives` 35 / 44, `na+distance` 25 / 26; se `na` 28 / 29, `na+relatives` 46 / 65, `na+distance` 38 / 47. At the
curve and the target the same order.

- **Congener groups (1)** change the models with the old features by ±0.003, within noise. They are kept on by default:
  real samples have congeners, and the relatives features need them.
- **The relatives features (2, 3) carry information but do not pay at protal's knobs here.**
  - They raise the cross-validated F1 (pe 0.9877 → 0.9891, se 0.9852 → 0.9877), the log loss, AP and the test sets'
    best-threshold F1.
  - Yet at 0.5, at the curve and at the target, test F1 falls (by up to 0.005) and 1.6-2.3 times as many minor
    congeners are missed: the thresholds chosen on the training design do not carry to the test design.
  - Alone, the four distance features are the best paired-end set at every knob (+0.001 to +0.006), and single-end
    within noise, still with 1.3-1.6 times the minor-congener misses.
  - `em_own_share` is among the forest's top five features, but the EM group alone helps no knob.
  - On this world false positives are few (20-45 per test set), while at GTDB r226 they were the larger error and the
    rank features gained 0.009-0.014 offline. So the set is not the default; an r226 build with
    `--features normalized+adjacency+relatives` (or `+distance`) decides.
- **The singleton rule (4)** without the read's condition, as first implemented, was wrong here. It vetoed 4-13 present
  taxa per test set and removed 2-11 true calls against 2-4 false ones (`singleton_count_only.txt`). In training it had
  vetoed 1,035 and 1,099 rows with none present: there abundances are less uneven (σ 1.3), so a minor congener rarely
  has one fragment beside a hundred. The present ones read like their own reference (identity ~0.975, EM share ~0.98,
  close share 0), the vetoed absent ones like the congener's (0.93, 0.45, 2.2; `singleton_features.txt`). With the
  condition (EM share below 0.5, or identity below 0.95) the rule vetoes 364-402 rows with 0-1 present and removes
  0-4 false calls and no true one (`singleton_refined.txt`). The models with the relatives features call none of these
  rows themselves.
- **Calls at a target share of false calls (5)** were 0.0005-0.004 below the knob curve for both read types.
  - On test samples deeper than any trained (models trained on samples of up to 200,000 read pairs, tested at 1M), pe
    was 0.003 below and se 0.002 above.
  - The target chosen on species held out (0.005-0.015) moved between models. The F1 over the targets is flat there, so
    small changes of the training rows move it.
  - So `build_gtdb_database.py` keeps `--call-mode curve` by default. protal uses calibrated calls whenever a model
    carries them.

## Defaults after the validation

| | default | opt-in |
|---|---|---|
| congener groups in training and test samples | `build_gtdb_database.py --congeners 0.25:2-5` | `--congeners 0` for none |
| relatives features | not in the default set (`normalized+adjacency`) | `--features normalized+adjacency+relatives` or `+distance` |
| abundance-weighted EM | computed for every taxon (in the dump; used by the singleton rule) | as a feature in `+relatives` |
| singleton rule | on, with the read's condition (`--singleton_congener 100`) | `--singleton_congener 0` |
| calls at a target share of false calls | protal uses them when a model has them | `build_gtdb_database.py --call-mode fdr`; `random_forest_cmdline.py --fdr-calls`; protal `--fdr F` |

## Not covered

- GTDB scale: the r226 tables of the earlier report lack the dump's new columns, so a new r226 build is needed to test
  the relatives features and the calibrated calls where false positives dominate. A sensible run: `--features
  normalized+adjacency+distance --call-mode fdr` against the defaults.
- Long reads (pb, ont): the features are computed for them too, but they were not trained or tested here.
- Real samples.

## Website

The website lists protal's options: `--fdr` and `--singleton_congener` are new, and the singleton rule changes calls of
databases' existing models by default (it is a rule, not a model feature). The website is not updated.
