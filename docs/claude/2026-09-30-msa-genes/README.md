# Which genes reach the strain MSAs

Date: 2026-09-30. Branch `strain-fixes`, commit `a1369eb` (the strain fixes L–P of the
[strain audit](../2026-09-29-strain-audit/README.md), and the long-unique threshold below). WSL
Ubuntu 24.04, 8 cores, 6 threads.

The question: are a species' genes without unique k-mers still left out of its MSAs, and can
the ones that go in be trusted?

## How protal chooses the genes

`SelectGenesForTaxon` (`src/RunProtal.h`):

1. It takes the genes with reads in the species' MSA samples. Before step O it took the "hittable" genes, those with unique k-mers.
2. For a species with relatives in the database it then drops the genes without long unique k-mers: a relative may share such a gene unchanged, and its reads would show as a second strain.

A long unique k-mer is one unique to its taxon whose 15-mer core also holds other k-mers (a flex block). Whether the filter applies depends on the share of the species' genes that have them.

## Finding 11, fixed

Step O's choice by reads did not bring back *Dummya solo*'s missing genes: the long-unique filter dropped them again.
- **Old rule:** the filter applied once half of a species' genes had long uniques.
- ***Dummya*:** it has no relative in the database, and has long uniques in exactly half of its genes (59 of 118).
- **Species with congeners:** long uniques in 98–99% of their genes.

`a1369eb` applies the filter from 90%. On the audit's accuracy run A (`run_P2.sh`, label P2 against P = `ac0f773`):

| | P | P2 |
|---|---|---|
| *Dummya* MSA genes | 59 | 118 |
| called at 50x | 0.615 | 0.999 |
| called at 2x | 0.515 | 0.839 |
| accuracy | 0.99999 | 1.00000 |
| false alternative calls | unchanged | unchanged |

No other species' MSA genes changed.

## Genes without unique k-mers in larger databases

`gene_uniques.py` on three unpacked databases (`unpack_dbs.sh`):

| Database | Species | Genes | Without any unique k-mer | With long uniques | Dropped by the long-unique filter |
|---|---|---|---|---|---|
| audit world (`~/audit5/world`) | 8 | 878 | 0 | 80.5% | 4 |
| tuning world H2 (`~/tune/H2/protal_db`, GTDB-like, with congeners) | 765 | 85,929 | 3 | 100% | 4 |
| `~/protal-perf/db900` | 900 | 101,081 | 4 | 100% | 5 |

In the two larger databases every gene has long uniques: nearly every 15-mer core is shared, so every unique k-mer counts as long.
- **What the filter drops:** it then applies to every species, and drops only the genes with no unique k-mer at all.
- **How many:** a handful. Both indices say "checked single-entry uniques", so their uniqueness is real (the check exists since `78f7a92`; before it, every gene was hittable).
- **GTDB r226:** `gene_uniques.py` on its `unique_kmers.tsv` gives the number there. Species that share marker genes unchanged will have more.

## What reaches the MSAs

**Data.** 48 samples of the tuning world's `test_congeners` set:
- read lengths 100, 150 and 250 bp, at 30,000 and 100,000 pairs;
- 11–40 species per sample, congeners together;
- depth about 1–6x.

**Run.** Profiled with the H2 database and `--map` (`genes_run.sh`); qcmsa run on each species as protal runs it (`genes_qcmsa.sh`, `--reapply-hcov 1000`); counted with `msa_genes.py`.

234 species got an MSA, with 24,877 genes in their genomes:

| | Genes | Share |
|---|---|---|
| with reads in the MSA's samples | 24,659 | 99.1% |
| in the raw MSAs | 24,608 | 98.9% |
| after qcmsa | 19,935 | 80.1% |

The missing genes, by cause:

| Cause | Genes |
|---|---|
| no reads in any of the samples | 218 |
| long-unique filter | 1 |
| no position with enough reads | 50 |
| qcmsa coverage gate | 4,505 |
| qcmsa multi-allelic filter | 1 |
| no filtered MSA at all | 167 |

The 167 are two species whose samples, bar at most one, have fewer than 1,000 written bases (`--reapply-hcov`).

**Uniqueness does not decide which genes go in.** Genes with 1–9 unique k-mers reach the filtered MSAs in 64% of cases (11 genes), those with 10–99 in 78%, those with 100 or more in 81%.

The coverage gate's losses are almost all in species with few samples (`gate_split.py`). A gene needs 2 samples that each write 30% of it:

| MSA samples | Species | Genes the gate removes | Median sample depth |
|---|---|---|---|
| 2 | 120 | 31% | 1.4x |
| 3 | 60 | 9% | 1.5x |
| 4 | 23 | 2% | 2.1x |
| 5+ | 31 | 0% | 1.5x |

With two samples, either one failing drops the gene; that gene says nothing about how the two samples relate.

## Can the genes with few unique k-mers be trusted?

In this world a strain's substitutions are spread evenly over its genes. So in every raw MSA row, each gene's rate of bases that differ from the reference row should be near the row's median (`cell_trust.py`). A congener's reads in a gene would raise the rate, or add IUPAC codes. A cell is "suspect" when it holds 3 more differences than the row's median rate predicts, and more than twice that rate.

| Unique k-mers of the gene | Congener in the sample | Cells | Rate / row median | Suspect | IUPAC per kb |
|---|---|---|---|---|---|
| 1–9 | no | 6 | 0.73 | 0 | 0 |
| 1–9 | yes | 11 | 0.66 | 0 | 0 |
| 10–49 | no | 2,672 | 0.96 | 0.64% | 0.011 |
| 10–49 | yes | 2,103 | 0.95 | 0.62% | 0.093 |
| 50–99 | no | 10,315 | 1.00 | 0.38% | 0.011 |
| 50–99 | yes | 5,619 | 1.00 | 0.52% | 0.023 |
| 100+ | no | 26,746 | 1.02 | 0.17% | 0.009 |
| 100+ | yes | 13,976 | 1.02 | 0.33% | 0.020 |

Genes with few unique k-mers carry no excess differences, with a congener present or not.
- **Suspect cells:** 0.6% for genes with 10–49 unique k-mers, against 0.2–0.3% for 100+.
- **IUPAC codes:** the genes with 10–49 unique k-mers hold 0.09 per kb when a congener is present, against 0.01 when not. This is the only trace of a congener's reads, about one code per 11 kb, and qcmsa's cell and gene filters see such codes.

The reads that step N keeps (within `--depth_identity_margin` of the species' best) exclude congeners of a few percent divergence.

## Conclusions

- Genes without unique k-mers are no longer left out in any number: 1 of 24,877 genes here. Genes with few unique k-mers go in, and their cells are as good as the others'.
- What remains out is genes that no sample can support (no reads), and genes qcmsa's coverage gate removes in species with 2–3 samples. That gate is by design: a gene needs 2 samples.
- Not tested: genes without any unique k-mer, which these databases barely have. At GTDB scale the long-unique filter drops them all, since there every gene has long uniques. Whether a relative that shares such a gene is present differs per sample. If GTDB has many such genes, a rule per sample could let them in where no species sharing the gene is reported. That needs the species sharing each gene, which the database does not store.

## Data

In WSL (not in git):
- `~/genes_study/{world,tuneH2,db900}`: the unpacked databases;
- `~/genes_study/cong`: the run, its map, the merged manifest and the logs;
- `~/audit5/accuracy/prot_*_P2`: the P2 evaluation.

The scripts are in `scripts/`; the sample reads are the tuning world's (`~/tune/test_congeners/points/*/sim`).
