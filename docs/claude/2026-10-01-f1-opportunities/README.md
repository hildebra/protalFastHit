# Where protal 0.7.1's F1 goes, and what could still raise it

- **Date**: 2026-10-01.
- **Question**: what precision protal 0.7.1 reaches and by which mechanisms, where its F1 is lost, and which
  changes not yet made could raise it, tested where they can be.
- **Code**: protal 0.7.1 (`1c11a00`); its presence models refitted offline as its trainer fits them (64 trees,
  at most 128 leaves, balanced classes, the 28 normalized features of `scripts/model_features.py`), with the
  helpers of [the alignment-features report](../2026-10-01-alignment-features/README.md) (`exp_lib.py`).
- **Data** (WSL `~/bench071`, from [the v0.7.1 benchmark](../2026-10-01-v071-benchmark/README.md)): its 368
  runs on the operon world (`error_budget.py`); and the 0.7.1 pipeline's training tables (36 pe, 36 se, 12 pb and
  12 ont samples) and independent test set (18, 18, 6, 6 samples, another design and seed), with their profiles
  and SAMs, for the refits; a second training collection of the same design with seed 101
  (`scripts/more_training.sh`, `~/bench071/V071_training2`).
- **Scores**: F1 at the knob 0.5. Refits are scored by cross-validation with species held out (5 folds, the
  trainer's scheme) and on the test set by a forest fitted on all training rows, seeds 1-3; a test difference
  carries a 95% interval from a paired bootstrap over test samples.
- **Run**: `scripts/error_budget.py ~/bench071`, `test_errors.py`, `thick_misses.py`, `features_exp.py`,
  `own_cluster.py`, `qual_calib.py`, `combo.py`, `more_data_exp.py` (all with the `protal-db-build` Python); the
  tables are in `results/`.

## Summary

0.7.1's precision is high: 0.97-0.99 for every read type, about 2 false positives per sample, under 1% of the
reported abundance on them. Its F1 is limited by recall, and almost all of what limits recall is evidence the
classifier does not get:

1. **No reads.** At 1,000 read pairs 42-53% of the present species have no read at all; F1 0.63-0.71 against a
   ceiling of 0.69-0.75 for a perfect classifier. In real genomes the markers are about 3% of a genome against 46%
   in this world, so real samples give about 14 times fewer marker reads per species at the same depth.
2. **A close relative the database lacks.** Present species with a missing congener in the sample are missed two
   to three times as often (paired-end 8.2% against 3.5%), and they are nearly all of the misses with more than 10
   fragments (18 of 20). Their reads are 65-98% the relative's, at identities no window separates.
3. **Distant strains**: 59 of the benchmark's 79 paired-end misses with reads are strains 2-4% from their
   reference, against 52% of the present species (10,000 pairs, 2x150).

On the classifier side, nothing large is left. More trees, more leaves and gradient boosting gain nothing for
short reads. Neither do thresholds, the congener context in the sample, the gene-neighbour shares, or features of
a taxon's own read cluster. Features of the sample's depth lose up to 0.026: they do not carry over to depths the
training lacks. Small gains that held in cross-validation and on the test set:

| change | pe | ont | others |
|---|---|---|---|
| divergence beyond the base qualities (new feature) | +0.0021 (−0.0008, +0.0052) | +0.0045 (−0.0028, +0.0136) | not run |
| conservation pattern of the hit genes (new feature) | +0.0017 (−0.0007, +0.0044) | +0.0021 (−0.0019, +0.0073) | se +0.0009, pb +0.0017 |
| both | **+0.0028 (−0.0001, +0.0058)** | **+0.0068 (−0.0016, +0.0195)** | |
| knob per depth bin, chosen on the training calls | +0.0009 | **+0.0146 (−0.0021, +0.0447)** | pb +0.0074, se −0.0013 |
| twice the training samples | +0.0021 (−0.0006, +0.0049) | +0.0074 (−0.0008, +0.0189) | se −0.0004, pb 0 |

Test F1 change against 0.7.1's model, 95% interval. None is significant alone. The two features also raise
cross-validated F1, except the conservation pattern for ONT (0.9737 to 0.9735).

## How 0.7.1 gets its precision

- **Unique k-mers.** The index flags k-mers found in one species only (and those within one difference of another
  species', "long unique"); a taxon's rate of them per aligned kb (`lu_per_kb`, `lsu_per_kb`, the per-gene rates)
  is high for its own reads and low for a relative's.
- **Read identity**: `identity`, `top_identity`, and the depth from reads within 0.08 of the best.
- **Coverage pattern**: genes hit against what the fragments predict (`gene_presence_ratio`), `depth_cv`.
- **Alleles**: a relative's reads add variant and multiallelic sites, and allele-frequency profiles (`RAF*`).
- **Other candidates** (0.7.1): the reads' MAPQ (the score gap to the next candidate), and the shares that a
  congener or another genus fits as well (the ZA tag).
- **Read assignment** (0.7.1): a consensus keeps a pair's mates and a long read's genes on one taxon; a sure mate
  guides the other.
- **Training to reject relatives**: the training database leaves out whole clades and single species (290 of 765
  here), so the model sees relatives of missing species and learns to reject them.

## Where the F1 is lost

The benchmark's runs of 0.7.1 (`results/error_budget.md`), paired-end, full database, summed over each point's
samples:

| read pairs | present, called | present, no reads | present, seen, not called | false positives | F1 | best threshold: F1 | ceiling |
|---|---|---|---|---|---|---|---|
| 1,000 (2x150) | 284 | 218 | 14 | 2 | 0.708 | 0.13: 0.712 | 0.732 |
| 1,000 (2x100) | 242 | 249 | 34 | 4 | 0.628 | 0.08: 0.667 | 0.689 |
| 10,000 (2x150) | 700 | 44 | 16 | 9 | 0.953 | 0.35: 0.955 | 0.970 |
| 500,000 (2x150) | 339 | 0 | 1 | 9 | 0.985 | 0.94: 0.996 | 1.000 |
| 500,000 (2x100) | 790 | 0 | 3 | 20 | 0.986 | 0.74: 0.996 | 1.000 |
| 5,000,000 | 360 | 0 | 0 | 6 | 0.992 | 0.71: 0.999 | 1.000 |

The best threshold is an oracle, chosen on the same samples; the ceiling calls every present species with reads
and nothing else. At low depth the reads are the limit, at high depth the false positives. Pooled over depths,
the best single threshold per read type is within 0.004 of 0.5. A depth-aware threshold looks large here
(0.986 to 0.996 at 500,000 pairs), but it does not carry over for short reads: chosen on the training calls, it
gains +0.0009 (pe) and −0.0013 (se) on the test set. The test set's per-point optima range from 0.08 to 0.94 (pe
and se) with no pattern in depth.

**The misses** (the 0.7.1 test set, `results/test_errors.md`). Present taxa missed with 1 fragment: 17.5%; with
2-3: 7.1%; with more than 100: 2.0% (14 of 698). 69 of the 76 paired-end misses were simulated from another
strain than the reference.

| present taxa | pe missed | se | pb | ont |
|---|---|---|---|---|
| no congener of theirs missing from the database | 3.5% | 4.5% | 5.7% | 4.8% |
| a congener in the sample that the database lacks | 8.2% | 8.9% | 6.2% | 13.6% |
| misses with more than 10 fragments in the second group | 18 of 20 | 20 of 24 | 4 of 5 | 8 of 8 |

The 14 paired-end misses with more than 100 fragments (`results/thick_misses_pe.md`): 13 have a missing
congener in the sample, and 65-98% of their reads were simulated from another species (`own_cluster.py`, from the
read names). *Salixnebacter limoum* has 2,320 fragments, 96% of them the relative's. Of the reads within 0.04 of
the misses' best identity, 33-96% are borrowed: the relative is that close. The model sees an absent species living
on borrowed reads, and here it is wrong only because a few of the reads are the species' own.

**The false positives** are relatives too: 10 of the 19 paired-end test false positives are at genus level, the
rest at family to class; most rest on 1-3 fragments (the same as on the V2 test set,
[report](../2026-10-01-false-positives-v2-test/README.md)).

## What was tried

Test F1 change against 0.7.1's model (interval); cross-validated F1 in `results/features_exp*.md`,
`own_cluster*.md`, `qual_calib.md`, `combo.md`, `more_data_exp.md`:

| change | pe | se | pb | ont |
|---|---|---|---|---|
| sample's depth: log10 fragments, taxa with reads, the taxon's share | −0.0167 (−0.039, +0.003) | −0.0262 (−0.061, +0.002) | +0.0024 | −0.0032 |
| congeners in the sample: their number, the deepest one's depth ratio, top of genus, their share | −0.0014 | −0.0037 (−0.008, −0.0001) | +0.0008 | +0.0010 |
| conservation pattern of the hit genes | +0.0017 | +0.0009 | +0.0017 | +0.0021 |
| gene-neighbour shares (dumped by 0.7.1, not trained) | −0.0009 | −0.0006 | +0.0020 | +0.0019 |
| own read cluster within 0.01 / 0.02 / 0.04 of the best | −0.0019 / +0.0004 / −0.0021 | −0.0005 / +0.0001 / −0.0009 | | |
| identity shape: 98th percentile − median, share at 0.99 or more | +0.0005 | −0.0007 | | |
| divergence beyond the base qualities | +0.0021 | | | +0.0045 |
| knob chosen on the training calls (one) | −0.0013 | −0.0023 | +0.0021 | +0.0046 |
| knob per depth bin, chosen on the training calls | +0.0009 | −0.0013 | +0.0074 | +0.0146 |
| 256 trees / 512 leaves | 0.0000 / 0.0000 | −0.0004 / −0.0003 | −0.0009 / 0.0000 | −0.0005 / 0.0000 |
| gradient boosting (HistGradientBoosting) | −0.0011 | −0.0029 | +0.0020 | +0.0044 |
| 25% / 50% / 75% of the training samples | −0.0042 / −0.0043 / −0.0024 | −0.0003 / −0.0015 / −0.0004 | −0.0104 / −0.0112 / −0.0009 | −0.0080 / −0.0141 / +0.0013 |
| twice the training samples (+ seed 101) | +0.0021 | −0.0004 | 0.0000 | +0.0074 |

- **The conservation pattern**: the log ratio of the median depth of a taxon's conserved genes (factor below 1
  in `gene_conservation.tsv`) to its fast ones, and the conserved share of its hit genes. A missing relative makes
  the fast genes deeper ([gene-scaled margin report](../2026-10-01-gene-scaled-margin/README.md)).
- **Divergence beyond the base qualities**: per read, its differences per aligned base less the mean error
  probability of its bases (10^(−Q/10)); per taxon the median, and the share of reads above 0.02. Present
  taxa's reads exceed their predicted errors by 0.010-0.013, absent taxa's by 0.051-0.066. It separates a
  genome's divergence from the read's errors, which differ from read to read most in long reads.
- **Sample depth** gains in cross-validation (pe 0.9793 to 0.9842) and loses on the test set (0.9690 to 0.9520):
  the training samples have 1,000, 20,000 and 200,000 pairs, the test 500, 10,000 and 500,000, and the forest's
  depth-specific splits do not extrapolate.
- **Twice the data**: a second collection alone scores −0.003 to +0.006 against the first, so a training set's
  draw carries about that much noise; the pooled one gains for pe and ONT and not for se and PacBio.

## Opportunities, ranked

1. **Divergence beyond the base qualities and the conservation pattern as model features** (pe +0.003, ONT
   +0.007 together). Both are cheap in protal: the base qualities are read for SNP calling already, so a
   per-taxon median of each read's excess is one accumulator in the profiler; the conservation pattern needs the
   per-gene depths and `gene_conservation.tsv`, which every 0.7.1 database has. Then a pipeline run to retrain,
   and the GTDB r226 test set to confirm.
2. **A knob per depth for long reads.** For long reads the best threshold depends on depth: on the benchmark 0.15
   at 3 Mb and 0.63 at 90 Mb for Nanopore, 0.26 and 0.52 for PacBio (an oracle); thresholds per depth bin chosen
   on the training calls gained +0.015 (ONT) and +0.007 (PacBio) on the test set. The pipeline already writes each model's F1-optimal threshold;
   storing it with the model (a `knob` per read type in the database) and using it as the default is small. For
   short reads it does not pay.
3. **More marker sequence per genome** (the "no reads" class). In real samples this class is about 14 times larger
   than here (3% against 46% marker share). More genes per species (more universal single-copy genes, or
   species-specific ones), or any read rescue from non-marker regions, raises low-abundance recall directly; no
   classifier can. Testable first on a world with real genome sizes (`simulate_gtdb_release.py --genome_length
   4000000`), which shows the real error budget.
4. **Training data**: more samples gain a little for pe and ONT. A denser depth design (the collector's default
   1,000-500,000 in five steps rather than this run's three) would also let depth-aware thresholds and features
   be learned without extrapolating. The r226 build's design decides how much is left.
5. **Distant strains**: 75% of the benchmark's paired-end misses with reads (all depths) are strains 2-4% from
   their reference (52% of the present species at 10,000 pairs). The stress test found such reads aligning less often (0.344 against 0.405 of the mates of strains
   5% or more away); more seeds or a lower identity floor for reads with unique k-mers of the taxon would give them
   more of their own reads. Measure the aligned share by strain distance first.
6. **Cohort information** (not testable on independent communities): species recur across the samples of a study.
   A lower knob for species called confidently in other samples of the same run could rescue low-depth misses;
   it needs a cohort-structured simulation to measure, and it trades on the cohort's design.
7. **The missing-relative misses are an information limit.** Their reads are mostly a relative's, within 0.01
   of their best identity; the remedy is a more complete database, not a feature. Read reassignment (an EM
   over the ZA candidates) would help where both relatives are in the database (7 of the 32 V2 false
   positives), not here.

Not worth pursuing on this evidence: model capacity, gradient boosting for short reads, a sample-depth feature
(harmful without denser training depths), congener context, the gene-neighbour shares, own-cluster features.

## Caveats

One simulated world, whose genomes are mostly markers, whose relatives all differ at every gene by rate, and
whose read errors are ART's and pbsim3's. The gains found are of the size of a training set's own noise (±0.003);
the GTDB r226 build's larger test set, on real genomes, is where they should be confirmed before shipping.
