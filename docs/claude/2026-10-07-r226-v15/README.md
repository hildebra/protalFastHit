# The r226 v15 build: the sample's complexity, and what the models learn about species

**Data.** The r226 v15 build (SLURM 24027266, q512n10, 2026-10-07 06:28-11:10, 84 threads; console log in the request,
`local/protal_r226_v15_logs.tgz` and `_training.tgz`, extracted to `local/v15`). The logs archive arrived cut off
inside `model_logs/error_reads/pb/test`: everything tar would have written after it is missing, among it
`build_metadata.tsv`, the se and pb `scenario_predictions`, the pe and pb `calls.tsv.gz`, `error_reads.log` and the
classifier logs. The training archive is complete. The ont error reads and part of pb's are in WSL `~/v15`.

**Which commit.** `build_metadata.tsv` is lost. The build started at 06:28, 13 minutes after `ce85bd7` (the complexity
features `096567f` and the v14 report) and before the next commits (`049f7a1`, 07:46), so v15 is `ce85bd7`: boosting at
250 rounds and 0.1 (`15006c2`), the default features `normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity`
(55, v14's 52 plus `sample_log_taxa`, `sample_low_identity`, `sample_identity`), `--evaluation basic`, scenario weight
0.25 and one global knob, the error reads (`652ed53`). Its reads are v14's: ART for Illumina, `hifi_reads.py` for
PacBio and pbsim3 for Nanopore, with pbsim3's error-free Q93 reads. Since then `1b76405` and `acd4163` (long reads in
`simulate_metagenomes`, without those) and `c4cf9bc` (Illumina without ART) changed all of them: the next build's rows
will not pair with v15's.

**The same rows as v14.** [`compare_tables.py`](compare_tables.py) ([output](compare_tables.txt)): every training and
test table of v15 has v14's rows in v14's order (sample, taxon, truth), and every shared feature is identical but
`excess_scaled_median` (at most 1.5e-4 apart). v15 adds the three complexity columns. So every difference between v14
and v15 is the model's: the features it is given, 250 rounds at 0.1 instead of 500 at 0.05, and the knob it chose.

Scripts here, run in WSL (4 pinned cores, niced; the shared machine was busy) with the scikit-learn venv of the v13
report and the trainer of the working tree:

    python3 compare_tables.py --builds v14=local/v14,v15=local/v15 > compare_tables.txt
    python3 compare_builds.py --builds v14=local/v14,v15=local/v15 > compare_builds.txt
    python3 ../2026-10-06-r226-v14/samples_held_out.py --build local/v15 > samples_held_out_v15.txt
    python3 ablate.py --build local/v15 --out ~/v15/ablate --read-types pe,ont,se,pb -t 4    # run_ablate.sh: pe, ont, se's v15
    python3 ablate.py ... --read-types pb --variants v15,no_complexity,no_ref                 # run_queue.sh: no_ref, pb,
    python3 test_species_held_out.py --build local/v15 --ablate ~/v15/ablate --out ~/v15/tsho #   test rows by fold models
    python3 ablate_compare.py --dir ~/v15/ablate --build local/v15 > ablate_compare.txt      # and --fixed-knob 0.5, 0.8,
                                                                                             # --pairs no_ref:v15, pb
    python3 knob_by_familiarity.py --build local/v15 > knob_by_familiarity.txt               # and v13, v14
    python3 familiarity_variants.py --dir ~/v15/ablate > familiarity_variants.txt
    python3 design_vs_test.py --build local/v15 --read-type pe > design_vs_test_pe.txt
    python3 error_taxa_summary.py ~/v15/model_logs/error_reads/ont --heldout local/v15/heldout_species.txt > error_taxa_ont.txt
    python3 own_seed_check.py --heldout local/v15/heldout_species.txt <taxa.tsv files> > own_seed_check.txt

The refits' knobs and times are in [`ablate_knobs.tsv`](ablate_knobs.tsv); their row scores stay in WSL `~/v15/ablate`.

`ablate.py`'s `v15` variant reproduces the build's model exactly (pe and ont design test scores: max difference 0,
knobs 0.80 and 0.5). The ablations ran for pe, ont and PacBio; se's (20 minutes a variant on the busy machine) were
stopped after its `v15` refit: se has no gut or shallow-soil scenario, and the v14 report already showed that its knob's
move to 0.5 comes with 250 rounds at 0.1, not with the features.

## Summary

- **The sample's complexity features did what they were added for.** With whole samples held out, pe's shallow soil
  went from 0.908 to 0.936 (95% interval +0.002 to +0.073 over samples) and its shallowest sample from 0.758 to 0.919;
  pe's "samples held out" F1 at 0.5 rose from 0.9475 to 0.9526 (misses 3,526 → 2,308), pb's and ont's by 0.002. With
  species held out nothing moved (all read types within 0.0005). Isolated by refitting v15's table without them:
  pe shallow soil with samples held out +0.023, ont and pb gut +0.004 to +0.006 (intervals above 0, also with
  species held out), soil hold-out samples +0.001 (pe) to +0.007 (ont, interval spans 0).
- **They cost the design's test set nothing that the data can show.** pe 0.9649 → 0.9645 at the knob, but that is the
  knob: at a fixed 0.5 the features gain +0.003 (interval +0.0002 to +0.006), at 0.8 -0.001 (spans 0). ont's
  design test fell 0.9759 → 0.9723 (interval -0.009 to +0.001): -0.0016 from 250 rounds at 0.1, -0.0024 from the
  features, neither distinguishable from noise on 42 samples.
- **Soil itself barely changed** (hold-out samples: pe soil +0.001, shallow soil -0.001, se soil -0.001, pb +0.000 and
  +0.001, ont +0.004 and -0.002). The rows are v14's, so the soil errors are v14's classes: divergent strains,
  1-2-fragment taxa, near-identical relatives of novel species. The features aimed at the last (branch `fp-features`,
  `11cdf6b`) are not in this build.
- **The knob is the least stable part.** pe 0.82 → 0.80, se 0.73 → 0.5 (the flip the v14 report foresaw for 250 rounds:
  se soil with species held out -0.004, design test +0.001); refits of near-identical models choose 0.73-0.80 for pe.
  The design's test set prefers 0.5: per sample, 0.8 saves it as many false positives as the training design (0.8 and
  1.0) but costs it twice the present taxa (1.4 and 0.7). 90% of those are strains, and the test design has more of
  them (77% of present taxa another genome than the database's representative, against 68.5%):
  [the knob](#the-knob-strains).
- **The models barely learn the species.** The test rows of species seen present in training are best called at
  0.07-0.20, those seen only absent at 0.88-0.95, but scored by fold models that never saw the species they are within
  0.002 of the final model's F1 (pe, ont, PacBio): the difference is mostly in the reads, not memorised. The `ref`
  features, constant within a species, are a real gain too: +0.001 to +0.002 with species held out at a fixed knob:
  [below](#do-the-models-learn-the-species).
- **Sample-level features still extrapolate at the edges.** The pe samples worst with samples held out are now the ones
  with the lowest `sample_identity` in their scenario (shallow s4 0.903 against 0.929 with species held out; soil s2),
  and no training sample lies between the design's 4,600 taxa and soil's 9,800.
- **The error reads:** the logs archive was cut off inside them; pe's extraction failed (its log is in the lost part);
  and `error_reads.py` counts held-out species as "unseen" (the training database's `genome2tiid.tsv` lists every
  species), which inflates them 6-26x in soil. Without those, the unseen species are depth-limited, with a small
  aligner loss in Nanopore soil: [below](#the-error-reads).

## What changed from v14, as the builds report it

[`compare_builds.py`](compare_builds.py) ([output](compare_builds.txt)) scores both builds' own predictions on the
same rows, each at its knob, with a paired bootstrap over samples (95% interval of v15 - v14):

| | v14 | v15 | v15 - v14 (95%) |
|---|---|---|---|
| pe design test (knob 0.82 / 0.80) | 0.9649 | 0.9645 | -0.0003 (-0.0025, +0.0020) |
| pe soil hold-out | 0.9646 | 0.9655 | +0.0009 (-0.0014, +0.0044) |
| pe shallow soil hold-out | 0.9408 | 0.9396 | -0.0012 (-0.0033, -0.0001) |
| pe shallow soil, samples held out | 0.9076 | 0.9355 | +0.0279 (+0.0015, +0.0733) |
| pe soil, samples held out | 0.9610 | 0.9639 | +0.0029 (+0.0013, +0.0050) |
| pe gut, species held out | 0.9929 | 0.9951 | +0.0022 (+0.0009, +0.0035) |
| se design test (knob 0.73 / 0.5) | 0.9612 | 0.9622 | +0.0010 (-0.0015, +0.0032) |
| se soil, species held out | 0.9648 | 0.9611 | -0.0037 (-0.0047, -0.0028) |
| se soil hold-out (build summary) | 0.9645 | 0.9635 | |
| pb design test | 0.9791 | 0.9785 | -0.0006 (-0.0073, +0.0051) |
| pb gut, species held out | 0.9907 | 0.9942 | +0.0035 (+0.0015, +0.0052) |
| pb soil / shallow, samples held out | 0.9564 / 0.9516 | 0.9588 / 0.9539 | +0.0024 / +0.0023 (both span 0) |
| ont design test | 0.9759 | 0.9723 | -0.0036 (-0.0087, +0.0009) |
| ont gut hold-out | 0.9878 | 0.9937 | +0.0058 (+0.0024, +0.0102) |
| ont soil hold-out | 0.9439 | 0.9483 | +0.0044 (-0.0003, +0.0139) |
| ont shallow soil, samples held out | 0.9435 | 0.9471 | +0.0037 (+0.0009, +0.0073) |

The scenario sets have 3 (hold-out) or 6 (hold-in) samples, so their intervals are wide and a difference of one
community moves them; the design sets have 41-309 samples. The se and pb hold-out samples' predictions are in the lost
part of the archive (their pooled numbers are in the build's summary). With species held out every read type is
within 0.0005 of v14 overall (pe 0.9519 both, se 0.9628 → 0.9633, pb 0.9595 → 0.9596, ont 0.9542 → 0.9543).

## The complexity features, isolated

v15 differs from v14 in two things besides the knob: the three features and 250 rounds at 0.1. [`ablate.py`](ablate.py)
refits v15's table without the features (= v14's features at v15's settings; the v14 report's `gbm250` refit of v14's
table agrees within 0.001), [`ablate_compare.py`](ablate_compare.py) compares ([own knobs](ablate_compare.txt),
[at 0.5](ablate_compare_knob0.5.txt), [at 0.8](ablate_compare_knob0.8.txt)). v15 minus "without", 95% interval:

| | pe, own knobs (0.80 / 0.75) | pe, both at 0.5 | ont (0.5) | PacBio (0.5) |
|---|---|---|---|---|
| design test | -0.0030 (-0.0056, -0.0004) | +0.0028 (+0.0002, +0.0057) | -0.0024 (-0.0067, +0.0018) | -0.0013 (-0.0080, +0.0045) |
| soil hold-out | +0.0014 (-0.0017, +0.0042) | +0.0004 (-0.0011, +0.0020) | +0.0067 (-0.0004, +0.0193) | +0.0001 (-0.0006, +0.0014) |
| shallow soil hold-out | +0.0008 (-0.0018, +0.0028) | +0.0028 (+0.0014, +0.0037) | -0.0002 (-0.0016, +0.0007) | +0.0011 (-0.0007, +0.0028) |
| gut hold-out | +0.0004 | 0.0000 | +0.0038 (+0.0012, +0.0089) | +0.0021 (+0.0012, +0.0026) |
| gut, species held out | +0.0020 (+0.0007, +0.0035) | +0.0018 (-0.0002, +0.0040) | +0.0059 (+0.0036, +0.0092) | +0.0040 (+0.0013, +0.0060) |
| soil, samples held out | +0.0022 (+0.0006, +0.0039) | +0.0016 (-0.0022, +0.0052) | +0.0020 (-0.0017, +0.0074) | +0.0032 (+0.0003, +0.0078) |
| shallow soil, samples held out | +0.0228 (+0.0007, +0.0609) | +0.0146 (-0.0009, +0.0461) | +0.0036 (-0.0000, +0.0082) | +0.0005 (-0.0029, +0.0031) |

- The features' purpose, a sample whose depth lies outside its scenario's training samples, holds: pe shallow soil
  with samples held out +0.015 to +0.026 whatever the knob. Per sample ([`samples_held_out_v15.txt`](samples_held_out_v15.txt)
  against the v14 report's): s6 (4.86 log10 fragments) 0.758 → 0.919, s1 0.904 → 0.943, s2 0.925 → 0.934; s3 and s5
  unchanged; s4 0.901 → 0.903.
- Gut gains with species held out as well (ont +0.006, PacBio +0.004, pe +0.002, intervals above 0); why was not
  looked into.
- pe's design test at the knobs (-0.003) is the knob: without the features the refit chooses 0.75, with them 0.80. At
  one knob the features gain (0.5) or do nothing (0.8). ont's -0.0024 spans 0; boosting at 250 rounds and 0.1 had
  already cost ont -0.0016 in the v14 report.

Without `sample_log_fragments` (keeping the complexity; `no_depth` in the tables): pe soil and shallow-soil hold-out
+0.001 and +0.002 (intervals above 0), samples held out +0.003 (spans 0), but the design test -0.003 at 0.5 (interval
-0.0058 to -0.0004), and gut -0.001 with species held out; ont: soil hold-out +0.002, design +0.002 (spans 0), gut
-0.001. The depth still carries information for the design's samples; keep it.

## Where sample-level features still extrapolate

The six hold-in samples per scenario are all the soil-like samples the model has. Held out, a sample whose complexity
lies at the edge of its scenario's is scored by what the model learned elsewhere: in v15 the two worst pe samples with
samples held out are the ones with the lowest median identity of their scenario (`sample_identity`):

| pe sample | log10 fragments | log10 taxa | low-identity share | sample_identity | F1 samples / species held out | FN |
|---|---|---|---|---|---|---|
| shallow soil s4 | 5.27 | 4.15 | 0.0227 | 0.9847 (lowest of 9) | 0.903 / 0.929 | 534 / 287 |
| soil s2 | 5.83 | 4.43 | 0.0272 | 0.9855 (lowest of 9) | 0.954 / 0.960 | 308 / 231 |
| shallow soil s6 | 4.86 | 3.99 | 0.0178 | 0.9897 | 0.919 / 0.924 | 354 / 278 |

The design's samples span 0.5-3.7 log10 taxa (up to 4,600 taxa with reads, at 10M pairs), the soils' 4.0-4.6 (9,800 to
35,000): nothing in training lies between, where a moderately complex real sample (or a soil sample of 1-2M pairs)
falls. The model then has to interpolate between two clusters with one sample-level feature per few samples.

## Do the models learn the species?

The trainer's report warns that the `ref` features, constant within a species, "can name a species". The test rows
split by what the training table holds of their species ([`knob_by_familiarity.py`](knob_by_familiarity.py),
[v15](knob_by_familiarity.txt), [v13 and v14](knob_by_familiarity_older.txt)) look like it: pe's design test rows of
species present in training are best called at 0.20 (F1 0.987), of species only absent in training at 0.92 (0.911),
of species not in training at 0.08 (0.988); likewise every read type, and v13's forest without `ref` (0.37, 0.85,
0.47). Boosting sharpened it (v15's thresholds are further apart).

But the species seen only as absent are the congeners of held-out species, which the test samples' held-out species
spill reads onto again: their reads look like presence. [`test_species_held_out.py`](test_species_held_out.py) scores
the test rows a second time by the trainer's five species folds (each row by the fold model that never saw its species):

| | pe design test | pe scenario hold-out | ont design test | ont scenario hold-out |
|---|---|---|---|---|
| final model (as built) | 0.9645 | 0.9551 | 0.9723 | 0.9492 |
| fold models that never saw the species | 0.9642 | 0.9540 | 0.9737 | 0.9478 |
| present in training: final / fold models | 0.9718 / 0.9715 | 0.9685 / 0.9683 | | |
| only absent in training: final / fold models | 0.8982 / 0.8974 | 0.8519 / 0.8451 | | |

([`test_species_held_out.tsv`](test_species_held_out.tsv), each variant and familiarity; PacBio: design test 0.9785 /
0.9778, scenario hold-out 0.9594 / 0.9573.) Having seen a species is worth at most 0.002 of F1 overall: the test set is
not optimistic for that reason, and most of the gap between the thresholds is in the evidence. Within the rows of
species seen only absent the final model is ahead by 0.001-0.007 (pe) and 0.017-0.028 (PacBio, 114 fewer false
positives on them in the scenarios): it has learned a little that those species tend to be absent. In a real sample those
species' held-out congeners are in the shipped database and take their own reads, so that lesson is at worst a small
bias against calling them. Without `ref` (`no_ref`, [own knobs](ablate_compare_ref.txt),
[at 0.5](ablate_compare_ref_knob0.5.txt), [PacBio](ablate_compare_pb.txt)) pe loses with species held out (soil
-0.0012, shallow soil -0.0024 at 0.5, intervals below 0), on the design test (-0.0024 at 0.5) and with samples held out
(-0.004); ont 0 to -0.002 (shallow soil with species held out -0.0020, interval below 0); PacBio -0.002 on the soil
hold-out samples and -0.001 with species held out (intervals below 0), 0 on the design test. The v14 report's caution
about `ref` is answered: its gain carries over to species the model never saw.

## The knob: strains

pe's knob of 0.80 is chosen with species held out over the training rows (design rows weight 1, scenario rows 0.25),
where the design's rows prefer it too (0.9712 against 0.9687 at 0.5). The design's test set prefers 0.5 (0.9698 against
0.9645), also when scored by fold models (0.9682 against 0.9642). [`design_vs_test.py`](design_vs_test.py)
([output](design_vs_test_pe.txt)): at the depths both designs have (50k and 200k pairs) the knob gains 0.003 on the
training rows and loses 0.012-0.013 on the test rows. Over all depths, 0.8 against 0.5 saves the training design 305
false positives and costs it 212 present taxa (309 samples: 1.0 and 0.7 per sample), the test design 64 and 106 (78
samples: 0.8 and 1.4). Of the present taxa lost, 90% are strains, another genome than the representative, in both sets;
the test design's present taxa are more often such strains (77.1% against 68.5%; in-silico strains 25.6% against
22.2%), which explains part of the doubled loss (they may also be more divergent; not looked into). The knob trades
strains for false positives at the training design's rate. Real samples are mostly not the GTDB representative genome,
so they are closer to the test design than to the training design here.

Knobs chosen by refits of near-identical models: pe 0.73 (without `ref`), 0.75 (without complexity), 0.79 (without
depth; and the v14 report's `gbm250` on v14's table), 0.80 (v15), 0.82 (v14). With species held out the training F1
varies by 0.001 over 0.70-0.85 (`trained_model.thresholds.tsv`: 0.9576, 0.9579, 0.9584, 0.9574); the design test moves
by 0.003 across that range. se moved from 0.73 to 0.5 by the 0.002 rule.

## The error reads

- **The archive.** `model_logs/error_reads` (se 4.3 GB, PacBio 2.4 GB, Nanopore 2.2 GB of SAMs) sits inside
  `model_logs`, so tar packed it before the files after it; the download was cut off inside PacBio's and everything
  behind it was lost (see Data). The se error reads did not arrive at all. What to archive instead:
  [2026-10-07-share-logs](../2026-10-07-share-logs/README.md).
- **pe failed** ("taking the reads of the pe model's errors failed (1)", after 5 minutes); `error_reads.log`, which
  says why, is in the lost part.
- **"Unseen" includes the held-out species.** `error_reads.py` takes the training database's species from
  `training_db/genome2tiid.tsv`; `gtdb_to_protal_db.py --from_db --exclude_species` copies that file unchanged (it keeps
  every species in the taxonomy files and drops only the held-out ones' genes). A held-out species, novel to the
  training database by design, is then "unseen" although no read of it can align to itself. [`error_taxa_summary.py`](error_taxa_summary.py)
  ([Nanopore](error_taxa_ont.txt), [PacBio, partial](error_taxa_pb_partial.txt)) splits them off with
  `heldout_species.txt`: per Nanopore soil sample 5,810-6,160 of the 6,040-6,470 unseen species are held out, per
  shallow-soil sample 5,850-5,970 of 6,980-7,260; per design sample about half.
- **The rest are mostly too rare.** Of the unseen species the training database has, per Nanopore sample: soil 233-310,
  shallow soil 1,131-1,295 (against 145-200 model FN), design 36-61 (against 1-2 FN). Most have own reads that left a
  record but aligned nowhere. [`own_seed_check.py`](own_seed_check.py) ([output](own_seed_check.txt); the three
  hold-out samples of each soil, 12 design samples each of Nanopore and PacBio) reads their records: only 6-11%
  (shallow soil), 11-29% (soil) and 0-11% (design; PacBio 0-2%) had an own read that seeded on their own genes; the
  rest seeded on other taxa by a few chance k-mers, as long reads from outside the marker genes do. The species was
  missed for want of a read on its markers. What is left for the aligner: 51-92 species per Nanopore shallow-soil
  sample (16-42 in soil), each with one read that seeded on its own gene and failed to align.

## The build

4 h 42 min (v14 5 h 09 min): the training simulations (4:05, from 0:13) are the critical path; collection 3:09,
training 2 min 32 s (v14 16 min: 250 rounds at 0.1, folds side by side), parity 24 s, relatives 2:31, error reads 12:27
(+5 min for pe's failure), `--add_model` 2 s (v14 2:49). `training_db` took 1:02 (v14 0:34) beside `protal_db` in the
background, off the critical path. The scratch peak of 474.2 GB (v14 252.7 GB) is not the error reads: the simulations
run ahead of profiling until 30 GB are left (`--keep-free`), and this node had 527 GB free (v14's 263 GB); the least
free space the collector logged was 130 GB.

## What to do next in training

1. **Keep the complexity features** (default since `096567f`); they fix the failure they were made for and help gut.
2. **Fill the gap between the design and the soils.** A scenario of moderate complexity (1,000-5,000 species) and soil
   at lower depths (0.1-0.5x), and more samples per scenario (6 → 10-12; the scratch peak follows the free space, so
   the node's disk is not the limit): sample-level features then interpolate instead of telling clusters apart, and an
   edge sample (s4) has neighbours.
3. **Revisit the knob.** It is chosen from a flat curve (0.73-0.82 among near-identical pe models, se flipping between
   0.5 and 0.73) and on rows with fewer strains than the test set; at 0.5 pe's test set gains 0.005. Options, cheapest
   first: a larger margin for leaving 0.5 (0.005 instead of 0.002) or a knob only if the test set (scored by fold
   models, so not optimistic) agrees; more strains in the training design (its present taxa 68.5% another genome,
   the test design's 77%); or a knob per sample complexity, now that the model knows it (v14 found no gain per taxa
   count, before the features).
4. **Test the false-positive features** (`fp-features`, `11cdf6b`, unmerged): soil's remaining false positives are the
   near-identical relatives they aim at, and nothing in v15 addressed them.
5. **Expect every read type's samples to change.** The next build at the current commit simulates PacBio and Nanopore
   in `simulate_metagenomes` without pbsim3's error-free reads (`1b76405`, `acd4163`) and Illumina reads without ART
   (`c4cf9bc`): none of its rows will pair with v15's. Compare it on species and samples held out and the test sets'
   F1, not row by row; a refit of v15's model on its tables (`ablate.py`'s `v15` variant) tells the model's share
   from the reads'.
6. **The error reads:** exclude the held-out species from "unseen" (`heldout_species.txt`, or the species with genes in
   `reference.map`), find pe's failure (`error_reads.log` of v15), and share a lean archive without the error reads
   ([2026-10-07-share-logs](../2026-10-07-share-logs/README.md)).

## Follow-up: pe's error reads, the knob, and what was implemented

**Why pe's error reads failed.** v15's `error_reads.log` (sent afterwards): after "the contigs of 54483 genomes in
109 s" the process pool broke (`concurrent.futures.process.BrokenProcessPool: A process in the process pool was
terminated abruptly`), a worker killed from outside: out of memory. se, pb and Nanopore ran to the end (se 327 s,
4.3 GB of SAMs). [`fragment_memory.py`](fragment_memory.py) ([output](fragment_memory.txt)) builds what `extract()`
held per tracked fragment: 1.2-2.3 kB in v15's layout (a set of reasons, and every aligned and seeded-on taxon in two
more sets). se's soil samples tracked 13.6M fragments each (the log: 81.7M in 6 training samples), so 17-32 GB a
sample; the 84 workers started on the largest SAMs first, pe's 18 soil and shallow-soil samples at once.
[`heldout_fragments.py`](heldout_fragments.py) ([Nanopore](heldout_fragments_ont.txt)): 85-99% of the tracked
fragments were own reads of held-out species counted as unseen. The fixes:

- `d5ffc20` (another session, the share-logs report): `--heldout`, so held-out species are no unseen ones; the SAMs
  capped to FP and FN.
- `2543408`: a fragment keeps its reasons as one tuple of shared strings in the dict pass 1 made, and of its taxa
  only the sample's error taxa (the taxa tables look up no other: the 7 earlier tests' tables unchanged): 273-300 B
  (measured as above), 4.5-8x less; the workers start while their estimates (300 MB and twice the SAM's size) fit in
  `--memory` (60% of the least of the machine's memory, `SLURM_MEM_PER_NODE` and the cgroup limit); a worker killed
  anyway loses no other sample (those beside it run again one at a time, one killed twice is named and left out).

**The knob is stable for one model.** [`knob_stability.py`](knob_stability.py) ([output](knob_stability.txt)) on
the refits of v15's tables: pe's best threshold over 500 bootstrap resamples of the training samples is 0.76-0.80
(5-95%), and 0.80 beats 0.5 in all of them; with the test set weighted as the training rows (design 1, scenarios
0.25) 0.80 wins there too (0.9591 against 0.9566; only the design's test rows prefer 0.5). The spread of 0.73-0.82
above is between models, each calibrated its own way, not noise within one. se's best is 0.70 (0.63-0.75), but the
0.002 rule kept 0.5 (+0.0015) and the test set prefers 0.5 too; Nanopore 0.60 (0.56-0.62) would gain a little on both;
PacBio 0.5. So the knob's real cost is the trade between the design's samples (strains: 0.5) and the soils' (false
positives: ~0.77), which no single knob settles; more strains in training put the strains' cost into the choice.

**Implemented** (branch `training-next` on `audit-fixes` `6a4b16c`, the merge of `fp-features`):

| Commit | What |
|---|---|
| `6a4b16c` | `fp-features` merged into `audit-fixes` (the 15 false-positive features in the default set, 70 features; the end-to-end sample thin with one pair of Mockella alpha); verified on the same merge onto `c2061fc` (`b5377d8`): unit, end-to-end (139), accuracy, script and trainer tests; the GTDB build test's rerun failure at its parity check is the one `6aee667` fixes |
| `2543408` | `error_reads.py`: compact fragments, a memory budget, killed workers lose no other sample (above) |
| `fb6c42d` | scenarios: each sample's depth 1/8-2x (`depth_range`), its species log-uniform and stratified from the scenario's range on a stream of its own (`simulate_metagenomes --species_per_sample N1,N2,...`), gut 150-1,000, soil 3,000-11,000, evenness sigma 1.0-2.5, a `moderate` scenario (1,000-5,000 species, 30% lacking), 10 hold-in and 4 hold-out samples |
| `c922d0f` | the global knob: the bootstrap median, kept if it beats `--knob` in 95% of the resamples and not worse on the test set; the training design's strains 0.5,0.2 (the test set's) |

What the next build changes at once: the reads of every type (`c4cf9bc`, `1b76405`, `acd4163`), 15 features, the
scenarios' samples, the strains and the knob rule. Its comparison with v15 is by held-out estimates; to tell the
parts apart, refit v15's models on its tables (`ablate.py`) and its models on v15's (for the shared features).
