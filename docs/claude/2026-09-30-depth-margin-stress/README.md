# Depth identity margin: fixed or anchored on the gene median, on a stress world

- **Date**: 2026-09-30.
- **Question**: the [0.7 vs 0.6.0a benchmark](../2026-09-30-v07-vs-v06/README.md) found that
  `--depth_identity_margin 0.04` undercounts strains that differ from the reference, and that 0.08, or a
  threshold anchored on the median over genes of each gene's median read identity, removes the bias; on
  its world (strains at most 4% from the reference) no rule beat a fixed 0.08. Which rule holds when
  strains reach the species boundary, when a minor divergent strain sits beside a dominant close one, and
  when relatives of species the database lacks are in the sample?
- **Code**: 0.7 = `audit-fixes` at `4b21427`, built in WSL (`~/fix-build`); the rules are an experimental
  build of it (`../2026-09-30-v07-vs-v06/scripts/depth_rule_patch.py`, the rule from
  `$PROTAL_DEPTH_RULE`, not in the repository), which profiles 0.7's SAMs again (`--profile_only`).
- **Machine**: WSL2, Intel Core Ultra 7 258V, 6 threads, 23 GB.
- **Data** (`scripts/world.sh`, `scripts/design.py`; in `~/stress`):
  - 120 GTDB-like species (`gtdb_like_lineages.py --species 120 --archaea 0.1 --seed 5`), 5 genomes each
    (`simulate_gtdb_release.py --genomes_per_species 5 --genome_length 200000 --strain_divergence
    0.002-0.03 --species_divergence 0.015-0.06 --seed 11`): every genome, the representative included,
    0.2–3% from its species' ancestor, so a strain is 0.4–6% from the reference (each species' most
    distant strain 1.5–6%, median 4.0%); congeneric species 3–12% apart at the markers.
  - Two databases built by 0.7: every species (`db_full`), and without a quarter of them (`db_missing`:
    30 species, each with a congener that stays in the database, `results/heldout.txt`).
  - 12 samples (`results/design.tsv`, simulated with `--from_manifest`): 6 at 2x100 (ART HS20, 300±40)
    and 6 at 2x150 (HS25, 350±50), 60 species each, 15 of each kind: the representative alone; one other
    strain alone (the strains in turn); the representative with the most distant strain as a 10, 20 or
    30% minor (`close_major`); the reverse (`far_major`). Species depth lognormal, median 12x, 2–80x.
- **Run**: `scripts/run.sh` (simulate, align once per database with `--no_profile`, profile each rule with
  `--knob 0 --no_strains`), `scripts/score_stress.py` (`results/scores.md`), `scripts/congener_diag.py`,
  `scripts/aligned_share.py`.
- **Scores**: with `--knob 0` every taxon with reads is reported, so the scores are of abundance, not of
  detection. Bray-Curtis between the true relative abundances of the species present (less the held-out
  ones against `db_missing`) and the reported ones renormalised over them; the median log2(reported /
  true) by kind of species.

## Result

No rule wins both cases. With every species in the database, the gene median and no filter at all are
best, and a fixed 0.08 below the top is slightly worse. With species missing, a fixed 0.08 is best by a
clear margin at 2x150, because the gene median follows an abundant relative's reads down and lets them
in. The current 0.04 is the worst everywhere, 2–5 times the best rule's error.

| Rule | all species, 2x100 | all species, 2x150 | species missing, 2x100 | species missing, 2x150 |
|---|---|---|---|---|
| no filter (`--depth_identity_margin 1`) | 0.040 | 0.020 | 0.048 | 0.066 |
| 98th percentile − 0.04 (0.7) | 0.171 | 0.102 | 0.162 | 0.101 |
| 98th percentile − 0.08 | 0.045 | 0.023 | 0.046 | **0.049** |
| 98th percentile − 0.12 | 0.040 | 0.020 | 0.048 | 0.066 |
| gene median − 0.04 | 0.046 | 0.034 | 0.052 | 0.053 |
| gene median − 0.06 | 0.041 | 0.024 | 0.046 | 0.056 |
| gene median − 0.08 | 0.040 | 0.020 | 0.046 | 0.063 |
| genes' 20th percentile − 0.06 | 0.040 | 0.020 | 0.046 | 0.065 |
| min(gene median − 0.06, 98th percentile − 0.08) | 0.041 | 0.021 | 0.046 | 0.058 |
| 98th percentile − (0.04 + 2 × (1 − gene median)) | 0.040 | 0.033 | 0.046 | 0.061 |

Bray-Curtis, mean of 6 samples. The share of reported abundance on taxa that are not present is 0.6–1.4%
with every species in the database and 7.6–9.3% with species missing, for every rule but 0.04 (0.3–0.7%
and 4.7–4.9%).

By kind of species, median log2(reported / true), all species in the database:

| Rule | representative alone | strain alone, <2% | 2–4% | ≥4% | close major, distant minor | distant major, close minor |
|---|---|---|---|---|---|---|
| no filter | +0.08 | +0.05 | 0.00 | −0.20 | +0.04 | −0.02 |
| 98th percentile − 0.04 | +0.44 | +0.17 | −0.25 | −0.81 | +0.25 | −0.61 |
| 98th percentile − 0.08 | +0.10 | +0.07 | 0.00 | −0.20 | +0.04 | −0.05 |
| gene median − 0.06 | +0.09 | +0.07 | +0.01 | −0.18 | +0.01 | −0.01 |
| gene median − 0.08 | +0.08 | +0.05 | 0.00 | −0.20 | +0.03 | −0.02 |

- **Mixtures**: a distant minor strain beside the representative is not cut by the gene median (+0.01 at
  0.06), so the risk raised in the benchmark did not show: at 10–30% of a species, the minor strain's
  reads that fall below the threshold are too small a part of its depth to show.
- **Relatives**: 86% of the gene median's extra error over a fixed 0.08 (species missing, 2x150) is on
  species with a held-out congener in the sample. *Vinaococcus venuris* (12.5x of its own beside a 70x
  held-out congener) is reported at 2.0× its true share with 0.08 below the top, 2.8× with the gene
  median and 2.7× with no filter. When a relative's reads outnumber a species' own on most genes, the
  gene median sits at the relative's identity and its threshold admits them all. In this simulation a
  congener's reads land on every gene, since all markers diverge at one rate; in real genomes
  relatives' reads would concentrate on the conserved genes, where a median over genes resists them
  better. That cannot be tested on this world.
- **Distant strains lose reads before any rule**: strains 4% or more from the reference are 13% low
  even with no filter. Their reads align less often: of the simulated mates of strains 5% or more away,
  0.344 align (primary records) against 0.405 of the representative's at 2x100, and 0.401 against 0.425
  at 2x150 (`aligned_share.py`). That is seeding and `-a` (0.9 by default), not the depth rule.

## The strain MSAs take reads by the same margin

Until this change the MSA's rows took the reads within `--depth_identity_margin` of the best ones too
(round 4, step N). The strain audit's runs A (one strain per species, 1–50x) and B (two strains per
species, minor 2–50%), aligned and profiled by 0.7 at margins 0.04 and 0.08 (`scripts/strain_check.sh`,
scored by the benchmark's copy of the audit's `evaluate.py`; run B by `scripts/mix_summary.py`):

| | 0.04 | 0.08 |
|---|---|---|
| run A, 10x: called cells of strain rows after qcmsa | 0.941 | 0.921 |
| run A, 10x: IUPAC codes in the raw MSA | 0.01% | 0.05% |
| run A, 1x: false alternative calls per million | 499 | 584 |
| run B: mixture sites called with both alleles, minor 20% / 30% / 50% | 0.71 / 0.92 / 0.98 | 0.72 / 0.94 / 0.99 |
| run B: wrong calls per million where the two strains agree, raw MSA | 46–62 | 512–582 |
| run B: the same after qcmsa | 10–23 | 24–38 |

The wider margin lets a relative's reads (the world's congeners are ~7% apart) into the strain rows:
10 times the false calls in the raw MSAs, about twice after qcmsa, which also removes 2% more cells.
So the two thresholds are separate from here on: `--depth_identity_margin` 0.08 for the abundance, and
a new `--msa_identity_margin` 0.04 for the MSA rows. With the new defaults, the MSAs of runs A and B are
byte-identical to those of margin 0.04 (48 files), and only the profiles change.

## Recommendation

Make 0.08 below the 98th percentile the default for the abundance (done: `--depth_identity_margin`,
with `--msa_identity_margin` keeping the MSAs at 0.04). It is within 0.003–0.005 of the best rule with every
species in the database and the best with species missing, which is the case the margin exists for.
The gene median is the better anchor for divergent strains alone, but a relative that outnumbers a
species pulls it down; before it replaces the fixed margin, it needs a guard against that (e.g. the gene
median only where the species' 98th percentile and gene median agree within a few percent), and a test on
real genomes, whose genes differ in conservation. The undercount of strains 4% or more from the
reference is at the alignment, and would be the next thing to look at for strains near the species
boundary.

## Follow-up: genes that differ in conservation

In the first world every marker evolves at one rate, so a relative's reads land on every gene, which is
why the gene median followed them. Real genes differ: ribosomal proteins are among the most conserved,
replication, repair and metabolism genes the least. A second world (`scripts/run_gcat.sh`, `~/stress2`;
binaries of `1b6161d`) is the first one with `simulate_gtdb_release.py --gene_rates categories` (new): each
marker evolves at its category's rate (ribosomal proteins 0.4, translation and transcription 0.8, tRNA
synthetases and modification 1.1, the rest 1.4) times a gamma draw of CV 0.35, scaled to a mean of 1, at
every level from the domain down to the strains (`results/gene_rates_true.tsv`). The design, the held-out
species and the samples are drawn as in the first world.

Two rules scale by the gene, with a factor r per gene (mean 1):

- per-gene margin: a gene's reads count if their identity is at least the 98th percentile − (M0 + r × K),
  wide on fast genes, narrow on conserved ones (`gtop:M0:K`);
- gene-scaled median (the gene median, scaled per gene): D = the median over genes of (1 − the gene's
  median read identity) / r, the strain's divergence from the reference measured on every gene; a gene's
  reads count if their identity is at least 1 − r × (D + K) − M (`gmedc:M:K`).

The factors come from the simulator's true rates, or are estimated from the database's genes as a build
could (`scripts/estimate_gene_rates.py`): for each species and gene, the median k-mer distance of the
other genomes' copies (`full_reference.fna`) to the representative's, over that species' median across its
genes; the median of that over species. Per gene, the estimate correlates 0.991 with the true rates;
averaged per category, 0.70 (the categories' own spread): ribosomal 0.37, translation 0.78, tRNA 1.24,
other 1.24 (`results/gene_rates_estimate.txt`).

Bray-Curtis, mean of 6 samples (`results/scores_gene_rates.md`):

| Rule | all species, 2x100 | 2x150 | species missing, 2x100 | 2x150 | sum |
|---|---|---|---|---|---|
| no filter | 0.049 | 0.026 | 0.052 | 0.037 | 0.165 |
| 98th percentile − 0.04 | 0.187 | 0.137 | 0.172 | 0.126 | 0.621 |
| 98th percentile − 0.08 (the default since `1b6161d`) | 0.058 | 0.033 | 0.057 | 0.034 | 0.182 |
| gene median − 0.06 | 0.052 | 0.030 | 0.052 | 0.034 | 0.168 |
| gene median − 0.08 | 0.049 | 0.026 | 0.051 | 0.035 | 0.161 |
| min(gene median − 0.06, 98th percentile − 0.08) | 0.053 | 0.028 | 0.053 | 0.034 | 0.168 |
| per-gene margin 0.03 + r × 0.05, estimated per gene | 0.056 | 0.029 | 0.056 | 0.032 | 0.173 |
| per-gene margin 0.03 + r × 0.05, true rates | 0.056 | 0.029 | 0.057 | 0.031 | 0.173 |
| gene-scaled median, M 0.03, K 0.03, estimated per gene | 0.050 | 0.028 | 0.051 | 0.032 | 0.160 |
| the same, per category | 0.050 | 0.028 | 0.050 | 0.033 | 0.160 |
| the same, true rates | 0.050 | 0.028 | 0.051 | 0.032 | 0.160 |

By kind (median log2(reported / true), all species): strains 2–4% from the reference get −0.010 with the
gene-scaled median, −0.015 with the gene median − 0.08, −0.019 with no filter and −0.032 with 0.08 below
the top; strains 4% or more away −0.19 to −0.24 with every rule (the reads lost at alignment); against
`db_missing`, the congeners of missing species +0.018 to +0.021 with the gene-median rules and no filter,
+0.046 with 0.08 below the top.

- **With genes that differ, the fixed margin is the worst of the reasonable rules.** Reads of a strain
  on its fast genes fall further below the top than 0.08, so it undercounts strains and, in relative
  terms, overcounts the reference-genome species; no filter does better (0.165 against 0.182).
- **The relatives' pull is gone.** A relative's reads now align mostly on the conserved genes, so the
  gene median, and a depth that is a median over genes, follow the species' own reads: even no filter
  scores 0.037 with species missing at 2x150 (0.066 in the first world).
- **Scaling the gene median by the gene helps a little** where relatives are present (0.032 against
  0.035 for the gene median − 0.08 at 2x150) and for strains 2–4% away; the factors estimated from the
  database's genes do as well as the true rates, and category averages nearly as well.
- **Over both worlds** (sum of the eight scores): the gene-scaled median 0.327 (its first-world scores are
  those of the gene median − 0.06, since its factors are all 1 there), the gene median − 0.08 0.330,
  no filter 0.339, 98th percentile − 0.08 0.344. The differences between the first three are small
  (0.002–0.008 per score); the fixed margin wins only in the first world with species missing.

**Revised recommendation.** Where genes differ in conservation, as real ones do, the gene median should be
the anchor: the gene median − 0.08 needs nothing but the reads and scores within 1% of the best rule
over both worlds; the gene-scaled median adds per-gene factors that a build can estimate from
`full_reference.fna` (0.99 correlation with the true rates here) for a further ~1%. The fixed 0.08 below
the top (the current default) should give way to one of them once they are checked on real genomes.
