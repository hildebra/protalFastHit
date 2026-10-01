# Where a missing relative's reads land, by how fast the gene evolves

- **Date**: 2026-10-01.
- **Question** (the user's): the model features of e680d91 were documented with "a relative's reads pile up on the
  fast genes". Biologically one expects the opposite: a species the database lacks differs least from its
  congeners on the conserved genes, so its reads should align there. Is the statement right, or an artefact of
  the simulations, and can real genomes tell?
- **Data**: the v0.7.1 benchmark world (`docs/claude/2026-10-01-v071-benchmark`, WSL `~/bench071`): 900
  simulated species, genes of different conservation (`simulate_gtdb_release.py --gene_rates categories`), and
  the 0.7.1 pipeline's 36 paired-end training samples (`V071/training/points/rl*`: 2x100, 2x150 and 2x250 at
  1,000, 20,000 and 200,000 pairs, 4 samples each), aligned by protal 0.7.1 to the training database
  (`V071/training_db`), which leaves out 290 species, 119 of them as single species whose genus stays. No real
  genomes: none are on this machine (every genome here is synthetic, `GCF_999...`).
- **Scripts**: `scripts/trace_relatives.py` (read by read: a read's name starts with its source genome, its
  primary record names the taxon and gene it aligned to; the gene's true rate from the simulator's
  `gene_rates.tsv`); `scripts/feature_by_class.py` (the two features as protal computes them from the profiles'
  gene logs, by class of taxon). Results in `~/bench071/results_conservation/`.

## The statement was wrong

**Read by read, a missing relative's reads align best on the conserved genes.** Per unit of their source's
coverage, relative to what a species' own reads give on the same gene (R: gene length and alignability cancel),
medians over the genes of each class of the simulator's rate:

| gene rate | genes | R, all primary records | on a congener | MAPQ < 4 | own reads MAPQ < 4 | R of the records the profiler keeps (MAPQ >= 4) | kept records on the most-hit taxon | relative's identity | own identity |
|---|---|---|---|---|---|---|---|---|---|
| < 0.7 | 58 | 1.01 | 0.77 | 0.70 | 0.11 | 0.38 | 0.89 | 0.958 | 0.983 |
| 0.7-1 | 33 | 0.77 | 0.92 | 0.46 | 0.03 | 0.43 | 0.83 | 0.936 | 0.977 |
| 1-1.4 | 41 | 0.54 | 0.98 | 0.35 | 0.01 | 0.36 | 0.81 | 0.928 | 0.972 |
| >= 1.4 | 36 | 0.19 | 1.00 | 0.23 | 0.01 | 0.14 | 0.84 | 0.919 | 0.964 |

Spearman correlation of a gene's rate with R over the 168 genes: −0.975. On the conserved genes the relative's
reads align as often as a species' own; on the fastest a fifth as often. But 70% of them on the conserved genes
fit several congeners equally (MAPQ below 4) and the profiler drops them, against 23% on the fastest genes; the
species' own reads on conserved genes are ambiguous too (11%, against 1% on fast genes) when congeners are in the
database. After the filter the relative's reads count about evenly up to rate 1.4 and less on the fastest genes;
they do not pile up on the fast genes either way.

**As the model sees it** (the features of e680d91, from the profiles' gene logs; factor below 1 in the
database's `gene_conservation.tsv` = conserved):

| taxa | rows | with 5+ hit genes | conserved / fast depth, log2, median (5+ genes) | conserved share of the hit genes (5+ genes) |
|---|---|---|---|---|
| absent, congener of a species held out of the database in the sample | 1618 | 840 | +0.33 | 0.44 |
| absent, other | 3794 | 1030 | +0.10 | 0.83 |
| present, a congener held out in the sample | 536 | 407 | −0.02 | 0.46 |
| present, no congener held out | 1835 | 1495 | +0.01 | 0.47 |

A taxon that holds only a missing relative's reads has its conserved genes deeper than its fast ones; present
taxa have them even. The feature carries what one expects biologically; the forest does not depend on the sign,
it learns whichever way the training data go.

**Where the wrong sentence came from.** The gene-scaled margin report (`../2026-10-01-gene-scaled-margin`) found
that a *present* species outnumbered by a held-out congener had conserved genes at 0.71 times the depth of its
fast ones, in a world of 12 genera of 10 species each. That is a ratio after the profiler's filters, where on
conserved genes the relative's reads (and the species' own) fit many congeners and are dropped; it says nothing
about where the reads align. Here, in a world of other design (genera of 1 to 38 species, 173 of them with one),
the same group is at −0.02 (0.99 times); why the two worlds differ was not traced.
The feature documentation turned it into a statement about where a relative's reads land, without measuring it.
Corrected in `docs/model-training.md`, `scripts/model_features.py`, the comments of `Profiler.h`, and noted in
the two earlier reports.

## What the simulation assumes that real genomes may not

- **One rate per gene at every timescale.** The simulator multiplies the divergence of every branch, within and
  between species, by the gene's rate, so a gene conserved within species is conserved between congeners in
  proportion. In real genomes the two correlate but need not be proportional (selection differs between the
  timescales, fast sites saturate, some genes recombine between species).
- **Uniform sites within a gene.** Substitutions fall anywhere (no stop codons, some codon indels); real genes
  have conserved stretches and variable ones, and divergence concentrated at third codon positions, which changes
  how many seeds a relative's read keeps and how many congeners it fits equally.
- **No recombination or gene transfer between congeners**: a marker copied from one species to another would be
  identical between them, i.e. ambiguous, whatever its rate.
- **Congener distances and counts** (3-12% at a gene of factor 1 here; 1 to 38 species per genus) are the
  world's design; the MAPQ share above depends on them directly.

So the direction of R should hold on real genomes (a relative aligns where it differs least), but how much of it
the MAPQ filter removes, and how often conserved genes are identical between congeners, needs real data.

## Follow-up: the drop before the MAPQ filter as a feature (f359d49)

The user's reading of the table above: the drop from conserved to fast genes is clear before the MAPQ filter and
fuzzy after it, so the records before the filter should show protal a congener it lacks. `f359d49` adds
`conserved_fast_record_ratio` (log2 of the depth of all of a taxon's best records, MAPQ 0 included, on its genes of
factor below 1 over that on its other genes) to the normalized features, and the three checks below. Validated as
0.7.3-dev with its own pipeline on the benchmark world, against 0.7.2-dev (`25d457e`, the same pipeline design and
seed; `../2026-10-01-features-depth-knobs/scripts/run.sh` with `VERSION=0.7.3dev TAG=v073 REV=f359d49`, and
`compare.py`; `results/benchmark_v073.md`).

**The signal is there** (the trainer's new section, `results/pipeline_feature_classes.txt`; median, quartiles):

| reads | present, no congener held out | present, beside a held-out congener | absent, congener of a held-out species | absent, other |
|---|---|---|---|---|
| pe | +0.03 (−0.19, +0.30) | +0.05 (−0.24, +0.36) | **+1.33** (+0.62, +2.26) | +3.19 |
| se | +0.01 | +0.03 | **+1.34** | +2.96 |
| pb | +0.10 | +0.14 | **+1.97** | +5.18 |
| ont | +0.06 | +0.10 | **+0.93** | +4.51 |

after the filter (`conserved_fast_depth_ratio`) every class has median 0. Its importance in the forests: 0.015
(ONT) to 0.046 (se).

**It does not raise F1.** 0.7.3-dev less 0.7.2-dev per benchmark sample (long reads: against 0.7.2-dev at `--knob
0.5`, as 0.7.3-dev has no depth knobs), mean (95% interval):

| reads | full database | missing database |
|---|---|---|
| pe | −0.0010 (−0.0031, +0.0007) | −0.0002 (−0.0026, +0.0021) |
| se | −0.0034 (−0.0063, −0.0005) | −0.0026 (−0.0067, +0.0003) |
| pb | −0.0013 (−0.0052, +0.0013) | −0.0006 (−0.0025, +0.0013) |
| ont | −0.0029 (−0.0092, +0.0011); FP +0.88 per sample | +0.0024 (−0.0002, +0.0051) |

The pipeline's own estimates move as little (species held out pe 0.9816 to 0.9833, se 0.9790 to 0.9805; test set pe
0.9726 to 0.9712, se 0.9666 to 0.9683). Refitting a forest with one feature more changes its trees, which alone
moves these scores by a few thousandths; only single-end with the full database is outside that, and against it.
Why the signal does not pay, presumably (not traced taxon by taxon): the absent taxa it marks are mostly ones the
forest already tells apart by its other features (`excess_median` alone is 0.048 for them against 0.013 for present
taxa), and where protal still errs, a
present species beside a congener it lacks, it hardly moves (+0.05 against +0.03): that species' own reads cover its
genes evenly and outnumber the relative's. A signal of an unknown congener, rather than a feature of the taxa that
carry its reads, would be a sample-level report: unreported taxa of one genus with a high
`conserved_fast_record_ratio` and `excess_median` point to a species of that genus the database lacks (not built).

**The checks on the benchmark world** (the pipeline's `model_logs/`): `gene_congeners.tsv`, 1,634 pairs of species
of 111 genera: the between-species factors correlate 1.00 (Spearman, 168 genes) with the within-species ones, conserved
genes 0.62, fast ones 1.42; the nearest congener's copy identical in 0.1% of species on conserved genes, within 1% in
3.4%, never on fast genes. The correlation of 1 is the simulator's assumption (one rate per gene at every timescale),
so it confirms the check, not the biology; a GTDB build will measure it. `relatives_by_gene_conservation.txt`
(by the database's factors instead of the simulator's rates) repeats the trace above: R 1.02, 0.79, 0.58 and 0.20 by
class, 72% of the relative's records below MAPQ 4 on the most conserved genes.

PacBio's `excess_median` is about −0.97 for every class: pbsim3's PacBio reads here carry base qualities near 0, so
the expected error is near 1. The differences between classes remain (−0.969 present, −0.941 absent beside a held-out
congener), but the feature means something else than for the other read types.

## How real genomes could tell (implemented in f359d49)

1. **At `--build`, from the database itself.** The build already compares each gene's copies within species
   (`gene_conservation.tsv`, from `full_reference.fna`). The same k-mer distances between each species'
   representative copy and those of its congeners (`reference.fna`, genus from the taxonomy) would give, per
   gene, a between-species factor, its correlation with the within-species factor, and the share of species whose
   copy is (nearly) identical to a congener's: how often a gene is ambiguous. Cheap; a GTDB build gives it for all
   genera.
2. **During the training collection** (`build_gtdb_database.py`): the held-out species are real GTDB genomes and
   their reads land on real congeners; `trace_relatives.py` with the database's factors in place of the
   simulator's rates would write the table above for real genomes into `model_logs/` on every build.
3. **In the trainer's report**: the two features by class of taxon (`feature_by_class.py`), from the training
   table's `meta_novel_level`/`meta_relative_rank`, so each model's report shows what the feature carries in its
   training data.

## Reproduce

```
python3 docs/claude/2026-10-01-conservation-pattern/scripts/trace_relatives.py ~/bench071
~/micromamba/envs/protal-db-build/bin/python docs/claude/2026-10-01-conservation-pattern/scripts/feature_by_class.py ~/bench071
```

`feature_by_class.py` reads the 0.7.1 database's `gene_conservation.tsv` unpacked by
`docs/claude/2026-10-01-f1-opportunities/scripts/features_exp.py` (`protal --unpack_db`).
