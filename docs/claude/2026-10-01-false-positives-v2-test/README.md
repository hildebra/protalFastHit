# False positives of the presence models, traced to their reads

- **Date**: 2026-10-01.
- **Data**: the independent test sets of the `V2` build of
  [Models per read type and a test set](../2026-09-30-read-type-models.md) (`~/tune/V2` in WSL): the
  GTDB-like tuning world (765 species, congeners 3-12% apart at the markers, strains 0.2-2% from the
  representative), 48 test samples (pe and se 18 each at 500, 10,000 and 500,000 read pairs; pb and ont
  6 each at 0.15, 3 and 90 Mb), 10-300 species per sample, profiled against the training database
  (`training_db`, which lacks the species of `heldout_species.txt`: single species and whole genera,
  families, orders, classes and phyla).
- **Code**: the profiles, dumps and models are those of the `V2` run (branch `audit-fixes` at `5385719`).
  This analysis reads them only; it ran at `80c317f`.
- **Run**: `python3 trace_false_positives.py ~/tune/V2 .` (pandas, numpy, `zstdcat`), output in
  [`output.txt`](output.txt), tables [`fp_calls.tsv`](fp_calls.tsv) (every false positive and its closest
  simulated relatives), [`fp_reads.tsv`](fp_reads.tsv) (their reads by source genome),
  [`fp_main_source.tsv`](fp_main_source.tsv).

A false positive is an absent taxon that the read type's model calls present (p ≥ 0.5, the knob). Its
reads are traced in the sample's SAM: the reads aligned to the taxon's genes, grouped by the genome they
were simulated from (paired- and single-end read names carry the genome's accession, long-read names its
index in the manifest), with the identity of each alignment (matches over aligned columns).

## How many

| read type | samples | absent taxa seen | false positives | per sample | of absent taxa |
|---|---|---|---|---|---|
| pe | 18 | 2,164 | 11 | 0.61 | 0.51% |
| se | 18 | 1,923 | 14 | 0.78 | 0.73% |
| pb | 6 | 156 | 2 | 0.33 | 1.28% |
| ont | 6 | 341 | 5 | 0.83 | 1.47% |
| all | 48 | 4,584 | 32 | 0.67 | 0.70% |

The 32 calls are 28 events: 4 are the same taxon in the same sample called from its paired-end and its
single-end reads. None is at the shallowest paired-end depth (500 pairs); per depth the rates are in
`output.txt`. Raising the knob trades them for missed species: at 0.7 there are 13 false positives and 328
false negatives instead of 32 and 225; the F1-optimal thresholds of the four models on this test set
are 0.47, 0.34, 0.48 and 0.14, below the knob, not above it.

## Where their reads come from

Every false positive is a database species whose reads came from a close relative in the sample:

| main source of the reads | false positives |
|---|---|
| a species held out of the database, in the same genus (or in 2 cases a sister genus) as the called taxon | 21 |
| a species held out with its genus or family, sharing a family or an order with the called taxon | 4 |
| a species the database has (a congener of the called taxon, present in the sample) | 7 |
| no relative (reads from unrelated genomes) | 0 |

- **A missing species' reads land on its congener (21).** The classic case: the database lacks the
  sample's species, and its reads align to the closest species the database has, at 0.90-0.96
  identity. For example *Halthivibacter licoa* (pe, 500,000 pairs): 156 of its 171 reads come from the
  three genomes of *H. terzuis*, held out of the database. *Kawoplasma romihalensis* (pe, 10,000 pairs):
  8 reads, all from two held-out congeners.
- **Deeper holdouts give few (4).** Species held out with their whole genus or family hardly produce
  false positives: their reads reach the database's nearest species, which shares only a family or an
  order with them, at 0.91-0.96 identity, and 1-5 of them do. This matches the rates by relation: 1.06%
  of absent taxa whose closest simulated species is a held-out congener are called, 0.69% at family
  level, 0 at order or above.
- **Congeners that both are in the database (7).** A present species' reads that align better to a
  congener than to their own reference: 1-5 reads at 0.93-0.98 identity (*Vipehalella nefiensis* from
  *V. vihalis*, *Lirotervibrio neyaltulis* from three *Lirotervibrio* in the sample). The sample's strain
  is up to 2% from its own reference, so for some genes another species' reference is as close. 0.55% of
  absent taxa next to a present congener are called.

Some calls pile up stray reads from several relatives at once: *Literbacter nagenus* got 5 reads from 5
genomes of 4 genera.

## Why the model calls them

**The evidence is thin.** The median false positive rests on 3 reads (18 of 32 on at most 3, only 3 on
20 or more). The model calls taxa from very few reads because real low-abundance species look like
that: of the present taxa with a single read hit, 71% are called; of the absent ones, 1%.

| read hits on the taxon | absent taxa | called (FP) | FP rate | present taxa | called | sensitivity |
|---|---|---|---|---|---|---|
| 1 | 1,232 | 12 | 0.97% | 483 | 343 | 71.0% |
| 2-3 | 798 | 12 | 1.50% | 668 | 618 | 92.5% |
| 4-8 | 828 | 5 | 0.60% | 512 | 494 | 96.5% |
| 9-48 | 1,108 | 1 | 0.09% | 667 | 657 | 98.5% |
| 49+ | 618 | 2 | 0.32% | 1,869 | 1,862 | 99.6% |

**Identity separates them only partly.** With up to 8 hits, the called false positives align at 0.952
mean identity (top read 0.963), the true positives at 0.970 (0.984): the relatives' reads are a little
worse, but a quarter of the true positives are at 0.960 or below (strains 2% from their reference,
sequencing errors), where most false positives are. With long reads the gap closes further: the ONT
false positives align at 0.86-0.95, and a quarter of the ONT true positives with few reads align at
0.938 or below (median 0.957).

**The few thick ones are borderline.** The three false positives with 20 or more reads (*H. licoa* 171,
*Kadococcus serpraixum* 81, *H. belvennais* 21, all from an abundant held-out congener at 0.93-0.94
identity) are called at p 0.53-0.59. Here the evidence is not thin but consistently a little too
divergent: the true positives with 20 or more reads align at a median of 0.972 (pe) and 0.975 (se),
and 95% of them above 0.953 and 0.956.

## What this means

- The false positives are what the training design is built to produce and the model is trained to
  reject: relatives of species the database lacks, and congeners. The model keeps them to about 0.7 per
  sample and 0.7% of the absent taxa it sees, at 1-3 reads mostly.
- This world's congeners are 3-12% apart at the markers by construction; how the rates carry over to
  GTDB r226's real congeners shows in the r226 build's report (section "False positives and false
  negatives by taxonomic rank"), and the same trace runs on its test set (`trace_false_positives.py
  OUTDIR`), as long as the test points' SAMs are kept.
- To lower them without losing the single-read sensitivity, two directions follow from the traces, both
  untested: a feature for the identity a taxon's reads should have at its depth (many reads, all at
  0.93-0.94, is a relative, not a strain), which helps the thick cases; and for thin calls, whether the
  taxon's reads also fit a congener the sample has (the 7 cases from present congeners).

The rank each source shares with the called taxon is in `fp_reads.tsv` (`shared_rank`), the true
positives' identity by read type and read hits in `output.txt`.
