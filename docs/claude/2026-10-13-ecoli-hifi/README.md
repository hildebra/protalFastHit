# Why HiFi reads of *E. coli* K-12 are called *Escherichia sp004211955*

2026-10-13, on the local mock test's data ([2026-10-13-mock-local-test](../2026-10-13-mock-local-test/README.md)):
- **Database:** the reduced r226 database of 761 species, built at `4dfc549`, kept in WSL as `~/Mocks/reduced_db_r226`.
- **Alignments:** the SAMs of run `0.8.3_4dfc549`, 1 Gb per long-read sample.
- **Features:** re-profiled from those SAMs with `--profile_only --profile_truth` by protal at `d865681`. Its protal
  code is that of the database's build plus the output-folder change; the probabilities match the run's to the digit.

## The question

In both HiFi samples of the test, *E. coli* is not called, and the congener `s__Escherichia sp004211955` is called
instead:

| Sample | Strains | *E. coli* depth | *E. coli* P |
|---|---|---|---|
| ATCC MSA-1003 | *E. coli* K-12 MG1655, 10% of cells | 46.5× | 0.122229 |
| Zymo D6331 | five *E. coli* strains | 25.1× | 0.122229 |

The same species is called from ONT reads (P 0.59) and from paired-end reads (P 0.80 in the HKU sample).

**The genomes.** r226's *E. coli* representative is GCF_003697165.2, the type strain ATCC 11775 (phylogroup B2).
K-12 is in r226 as a non-representative *E. coli* genome; the reduced database has 5 other *E. coli* strains as
alleles, none of them K-12. *Escherichia sp004211955* is *Escherichia* sp. E1V33 (2 genomes in r226).

## The reads

**K-12's reads in the ATCC sample** ([sam_taxa.py](sam_taxa.py),
[results/sam_taxa_atcc_hifi.txt](results/sam_taxa_atcc_hifi.txt)). That mock holds no other *Escherichia*, so every
*Escherichia* record is a K-12 read.

| Taxon | Reads | Genes | Mean identity (10th percentile) | Mean MAPQ | Records within 5 edits of another taxon |
|---|---|---|---|---|---|
| *E. coli* | 2,835 | 108 | 0.984 (0.970) | 11 | 2,472 of 5,559 |
| *sp004211955* | 470 | 51 | 0.989 (0.975) | 8 | 585 of 925 |
| *E. fergusonii* | 458 | 56 | 0.983 (0.948) | 6 | 632 of 870 |

**Three patterns** ([lost_reads.py](lost_reads.py), [results/lost_reads.txt](results/lost_reads.txt)):

1. **K-12 is not close to the representative.** Its reads are 98.4% identical on average, so no gene-long HiFi
   alignment (1-4 kb) is perfect. The best reads' identity (`top_identity`, the 98th percentile) is 0.9957, where
   every other species in the sample has 1.0.
2. **Some of K-12's genes are nearer a congener than its own representative.** On the neighbouring genes of the
   *secE-nusG-rplKAJL-rpoB* locus (database genes 3 PF00466, 68 TIGR00922, 72 TIGR00964, 92 TIGR01169, 99 TIGR01632,
   103 TIGR02013), K-12's reads have no *E. coli* record at all. They align to *sp004211955* at 99.3%.
   - **Not missing genes:** the representative's copies are in GTDB's marker files and in the database (all 120 genes,
     with unique k-mers); K-12's copies of this locus are simply nearer E1V33's.
   - **Elsewhere:** genes 10, 18, 38, 97, 109, 111 and 113 lean to *E. fergusonii* or *sp004211955* in the same way.
   - **Whole reads follow:** a 10-kb read over such a locus goes whole to the congener ("locus": 237 reads).
3. **The read consensus has a cycle.** `ConsensusOfRead` (src/Alignment/ReadConsensus.h) compares two taxa on the
   genes where both have a candidate. A taxon that scores less on them is beaten, and the consensus is a taxon nobody
   beats. When every taxon is beaten by one, it takes the one with candidates on the most genes.
   - **The example, read 13443:**
     - *E. coli* beats *sp004211955* on their 6 common genes (by 1-3 edits on 5 of them).
     - *E. coli* is beaten by *E. fergusonii* by one edit on gene 13 in a 3-gene comparison.
     - *E. fergusonii* is beaten by *sp004211955*.
     - Nobody is unbeaten, so *sp004211955*, with candidates on all 8 genes, takes the read. Its genes are written
       with `ZR:i:1` and MAPQ 22.
   - **Count:** 101 reads lost this way ("cycle"), where *E. coli* fit at least as well on most genes and better on one
     at least.
   - **Another genus:** 36 of them went to *Kosakonia sp000814905*.
   - **Rest of the loss:** 288 more reads are mixed.

**The loss, per sample.** The congeners' calls are made of these lost reads:

| Sample | Kept | Lost | Locus | Cycle | Mixed |
|---|---|---|---|---|---|
| ATCC HiFi | 2,791 | 626 (18%) | 237 | 101 | 288 |
| D6331 HiFi | 1,703 | 402 (19%) | 199 | 58 | 145 |
| D6331 ONT | 2,270 | 485 (18%) | 233 | 85 | 167 |

## The model

[explain_call.py](explain_call.py) changes each of the taxon's features in turn to every reference taxon's value in
the same sample, and scores the result with the database's HiFi model (`model_PB.xml`, through
`scripts/model_pmml.py`; it gives protal's 0.122229).
- **Outputs:** [results/explain_ecoli_pb.txt](results/explain_ecoli_pb.txt) and
  [results/explain_sp004211955_pb.txt](results/explain_sp004211955_pb.txt).
- **Side by side:** [compare_rows.py](compare_rows.py) puts the features next to each other
  ([results/feature_comparison.tsv](results/feature_comparison.tsv)).

**What *E. coli* would need** (against the sample's 17 true positives):

| Feature | *E. coli* | True positives (median) | P if set to a true positive's value |
|---|---|---|---|
| `top_identity` | 0.9957 | 1.0 | +0.57 on average (to ~0.69: called) |
| `excess_scaled_median` | 0.014 | 0.0045 | +0.31 |
| `lu_per_kb` (unique k-mer hits per kb) | 28 | 139 | +0.26 |
| `lsu_per_kb` | 18 | 131 | +0.17 |
| `split_fragment_share` | 0.014 | 0 | +0.12 |

**The congener mirrors it.** *sp004211955* is called (P 0.55) with no unique k-mer hits at all (`lu_per_kb` 0) and 13%
of its genes. Its locus reads give it `top_identity` 0.998. Setting that to a true negative's value (median 0.963)
drops it to about 0.08.

**One feature decides both calls:** `top_identity`, measured against the representative only.
- **Known as a HiFi feature:** [features.md](../../features.md) has it among the top three of the PacBio models,
  "where a present species' best reads reach 1.0 and a relative's do not".
- **What it gets wrong here:** a strain 1.6% from its representative is taken for a relative, and the congener whose
  copies happen to match a locus of that strain for the present species.
- **What the database already knew:** *E. coli*'s strain alleles explain 84% of the reads' differences
  (`allele_explained_share` 0.84 against a true-positive median of 0.26; `allele_identity_gain` +0.016, so the reads
  are about 99.6% identical to their best alleles). The model hardly uses that; its identity features ignore the
  alleles.

**Why the other read types call it:**

| Read type | *E. coli* `top_identity` | Sample's true positives | P |
|---|---|---|---|
| ONT | 0.994 | 0.9975 (read errors) | 0.59 |
| Paired-end | 1.0 (some 150-bp reads match perfectly at 98% identity) | 1.0 | 0.80 |

## What this means at r226

The cause is not specific to the reduced database. At r226 *E. coli* has the same representative, and K-12's best HiFi
reads will match it at the same 0.9957.

What may differ at r226:
- **The HiFi model:** trained on many more real strains (with the build defaults since 2026-10-11, real strains are
  ~54% of the simulated strains), so it may weigh `top_identity` less.
- **The alleles:** up to 16 per gene, drawn from thousands of *E. coli* genomes, may include K-12-like ones.
- **The locus and the cycles:** these depend on the congeners, and r226 has all of them (*E. albertii*,
  *E. marmotae*, *E. ruysiae*, E1V33 and others). They will probably be the same or worse.

To check on the cluster: the ATCC MSA-1003 HiFi run (SRR11606871) and the Zymo D6331 HiFi runs at r226
(`download_mocks.py --sets atcc1003_hifi,d6331_hifi`).

## What could fix it

These are proposals; none is implemented. Each changes the features the models are trained on, so each needs a
rebuild and a comparison on simulated test sets and on these mocks.

1. **Read consensus without the cycle fallback.** When every taxon is beaten, take the one whose largest defeat is
   smallest (minimax over the pairwise margins on common genes), not the one with the most genes. A defeat by one
   edit on one gene would then no longer hand a read to a taxon that lost by many.
   - **Cost:** small, in ReadConsensus.h, and testable on the read level by itself.
   - **Effect here:** at most the 101 "cycle" reads per 1 Gb, including the reads sent to another genus.
2. **Identity against the best allele.** Compute `identity` and `top_identity` from each record's best allele where
   the species has alleles: the alignment already scores alleles (`--no_allele_scores` turns that off). A strain the
   alleles cover would then look like its species, and a congener would not borrow its locus.
   - **Biggest lever:** it changes the feature that decides both calls.
3. **More divergent real strains in the HiFi training data.** The model should see present species whose best reads
   stay below 1.0, with the allele features telling them from relatives.
4. **The locus effect itself** (a strain's genomic region nearer another species' copies) is biology. The allele index
   (`--index_alleles`, off by default) would give such reads a candidate of their own species, if an allele of that
   locus is in the database; its earlier F1 gain was too small for its cost.

## Files

| File | What it does |
|---|---|
| [dump_features.sh](dump_features.sh) | re-profiles the four samples' SAMs with `--profile_truth` |
| [truth_taxids.py](truth_taxids.py) | the truth files for that, from the benchmark's truth tables |
| [sam_taxa.py](sam_taxa.py) | per taxon and gene: records, identity, MAPQ, tags; mixed reads |
| [lost_reads.py](lost_reads.py) | a species' reads lost to congeners, by cause |
| [explain_call.py](explain_call.py) | which features change a model's probability for a taxon |
| [compare_rows.py](compare_rows.py) | taxa's features side by side |
| [results/](results/) | their outputs on this data |

Working folders `~/ecoli_dbg` (the unpacked database, the feature dumps) and `~/nofolders` (the build at `d865681`)
were removed afterwards. `dump_features.sh` with a protal build makes them again in a few minutes.
