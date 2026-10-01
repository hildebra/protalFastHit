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

## How real genomes could tell (not implemented)

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
