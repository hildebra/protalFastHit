# Why a real D6331 sample missed *Akkermansia muciniphila* with 1,111 reads, and spike genes left out of the evidence

2026-10-15. **The question** (the user's, after the [F1 audit](../2026-10-15-f1-audit-since-0.6.0a.md)):
*Akkermansia muciniphila* had 1,111 fragments on 38 genes at 2.3× in a real D6331 mock sample and was not called
(P 0.484). How did that happen, and is there a robust change in settings that recovers such cases?

**Short answer.**
- **The library is not clean shotgun data.** In that sample (Dohlman et al., SRR38158090, "D6331 DNA only"), every
  species of the standard has one to three marker genes holding 48-97% of its reads. *Akkermansia*'s 2,177 records:
  2,108 on one ~1.4 kb stretch of *secA* (TIGR00963) at identity 1.0, and 1-8 on each of 37 other genes. That is a PCR
  product or fragment carried over into the library (or a spike-in), not the genome's coverage. The same mock in two
  other labs' libraries (HKU, Auckland) has no such piles.
- **The coverage features read the piles as a relative's reads.** A relative the database lacks puts its reads on a
  few genes, and the model has learnt that `gene_presence_ratio` far below 1 and a large `depth_cv` mean "absent".
  Every present species of the sample scored lower: *Bifidobacterium adolescentis* with 4,719 fragments P 0.69,
  *Bacteroides fragilis* (4,505) 0.55, *Fusobacterium animalis* (6,248) 0.59, where species with that much evidence score
  near 1 in simulations. *Akkermansia* (0.48) and *E. coli* (0.38) fell below the knob.
- **No setting fixes it robustly.** A knob of 0.45 calls *Akkermansia* but not *E. coli*, and adds a false positive in
  another sample. The model kind, the gates and the fragment count do not touch the coverage features.
- **A sample-level spike mode fixes it:** before a sample is profiled, protal counts each gene's records, finds each
  taxon's spike genes and checks whether the sample shows the pattern of carried-over products (spikes on several of its
  most abundant taxa, a share of its reads). Only then is the sample in spike mode, and its spike genes are left out of
  the evidence (as suspect gene copies are; `--spike_genes auto`, the default). The Dohlman sample is the only one of
  ten mock samples and 62 simulated test samples in spike mode. There *Akkermansia* (P 0.48 → 0.80) and *E. coli*
  (0.38 → 0.55) are called and every present species scores the same or higher; in every other sample not one probability moves
  ([section 4](#4-spike-mode-and-what-it-changes)). Without the sample check (`--spike_genes on`), the per-taxon rule
  also takes absent relatives' piles, and three such congeners were then called in two long-read samples.

## Setup

- **Code:** the mock runs of [2026-10-14-mock-sparse-gates](../2026-10-14-mock-sparse-gates/README.md) (protal 0.8.4 at
  `5369c99`, run `0.8.4_5369c99_dbgates-gbm`: the database's gates, pe C 4 stride; the alignments do not depend on the
  model kind), profiled again with `--profile_truth` by the same binary (`~/c4test/protal/build/protal`) for every
  taxon's features; spike mode was built on `dd01eb1` in the worktree branch `f1-spike-genes`, and every SAM profiled
  by that binary with `--spike_genes auto`, `on` and `off` (`off` is protal as before; the profiles of the first
  version, without the sample check, were made with its `--keep_spike_genes` and are cited as such).
- **Database:** the reduced r226 database of the mock benchmark (761 species, `~/c4test/db/protal_db`, built at `5369c99`
  with pe 150 bp, HiFi and ONT models, forests by default). Its training table `~/c4test/db/work/training/training_data.tsv`
  and its simulated test samples' SAMs (`~/c4test/db/work/test/points/*/protal/alignments`: 30 paired-end samples
  with truth, 32 long-read samples without).
- **Machine:** WSL on 4 cores (`taskset -c 0-3`), shared and loaded; no timing below is a benchmark.

Scripts here, in the order used:

| Script | What it does |
|---|---|
| [dump_features.sh](dump_features.sh) | profiles six pe mock SAMs again with `--profile_truth` (every taxon's features) |
| [feature_shift.py](feature_shift.py) | where the mock samples' present taxa lie among the training table's present taxa, feature by feature |
| [swap_features.py](swap_features.py) | a taxon's probability with one feature (or a group) at a time set to a reference: its row in another sample, or the training median |
| [duplicates.py](duplicates.py) | duplicate fragments per taxon of a SAM |
| [gene_spikes.py](gene_spikes.py) | taxa whose reads pile on one or two genes |
| [spike_rule.py](spike_rule.py), [spike_rule_all.sh](spike_rule_all.sh) | the rule protal applies, as a prototype over every mock and simulated SAM |
| [spike_compare.sh](spike_compare.sh), [spike_calls.py](spike_calls.py) | every SAM profiled in each spike mode, and the calls of two modes compared |
| [spike_check_long_sims.sh](spike_check_long_sims.sh) | protal's sample check on the simulated long-read test samples (no truth there) |

## 1. Every present species of that sample is depressed, not only *Akkermansia*

The model's probability for the 12 species of the standard in the Dohlman sample against the HKU sample (same standard,
another lab), forest model at the database's gate (`dump_features.sh`; the column "spike mode" from
[section 4](#4-spike-mode-and-what-it-changes)):

| Species | Dohlman fragments | P | P, spike mode | HKU fragments | P |
|---|---|---|---|---|---|
| *Bifidobacterium adolescentis* | 4,719 | 0.688 | 0.984 | 4,204 | 1.000 |
| *Limosilactobacillus fermentum* | 6,295 | 0.797 | 1.000 | 2,473 | 1.000 |
| *Faecalibacterium prausnitzii_E* | 5,024 | 0.641 | 0.656 | 6,086 | 0.891 |
| *Bacteroides fragilis* | 4,505 | 0.547 | 0.750 | 5,753 | 0.953 |
| *Fusobacterium animalis* | 6,248 | 0.594 | 0.719 | 2,422 | 0.922 |
| *Clostridioides difficile* | 738 | 0.609 | 0.609 | 499 | 0.938 |
| *Escherichia coli* | 713 | **0.375** | 0.547 | 2,669 | 0.781 |
| *Akkermansia muciniphila* | 1,111 | **0.484** | 0.797 | 1,083 | 0.875 |

What the features show (`feature_shift.py`, against the 1,014 present training rows of 200-20,000 fragments): in the
Dohlman sample every present species has `depth_cv` above the training's 99th percentile, and `gene_presence_ratio`
and `hit_gene_fraction` near its bottom; the other samples' present species do not. Setting *Akkermansia*'s coverage
features (gene_presence_ratio, breadth_ratio, hit_gene_fraction, depth_cv, gene_dispersion) to its HKU values raises
its probability from 0.484 to 0.859, `gene_presence_ratio` alone to 0.781; every other group of features moves it by
0.03 at most (`swap_features.py`).

Two other things the comparison showed, not the cause here:
- **A database trained on 150 bp reads only.** In its training rows `top_identity_bound` is 0.975 for every present
  taxon (the Wilson bound of a perfect 150 bp read), so the feature is a constant there and real reads of other lengths
  fall outside it. The GTDB builds now train on 75-250 bp reads.
- **The strain alleles' features** (`allele_explained_share`, `polymorphic_known_share`) are 0 in the training's present
  rows (representatives and in-silico strains) and ~0.37 in the real samples, the audit's E1 in action: the finished
  database's alleles come from every genome.

## 2. The piles: one gene per species, at identity 1.0

Not PCR duplicates: *Akkermansia*'s 1,112 fragments are 97% distinct (`duplicates.py`). Its records per gene in the
Dohlman sample: 2,108 on gene 71 (*secA*, TIGR00963), then 8, 6, 5, 4, 3, ... (37 genes); in HKU 20-58 on each of 114
genes. The 2,108 lie between positions 1,900 and 3,300 of the gene, NM 0, MAPQ 111.

`gene_spikes.py` over every mock SAM: in the Dohlman sample the top gene holds 97.5% of *L. fermentum*'s records, 92%
of *F. prausnitzii*'s, 88% of *E. coli*'s, 86% of *B. adolescentis*', 58% of *F. animalis*', 56% of *C. difficile*'s and
48% (plus 47% on a second gene) of *B. fragilis*'; absent *Veillonella dispar* has 18,152 records on one gene and
12,683 on another. In no other sample does a present species have such a gene; absent relatives do, with smaller piles
(their reads gather on the genes where they differ least from the present congener).

## 3. Settings

The calls of the Dohlman sample and the others near the knob (`dump_features.sh` output; present taxa below 0.6, absent
ones above 0.3):

- **Knob:** at 0.45, *Akkermansia* (0.484) and Treichel's *B. ovatus* (0.469) and hCom2's *P. vulgatus* (0.469) are called,
  and Treichel's absent *Roseburia_C* sp900315735 (0.469); *E. coli* (0.375) is not. A lower knob trades errors and moves
  every sample.
- **Model kind, gates, fragment count:** they leave the coverage features as they are. The GBM model of the same
  database missed *Veillonella rogosae* at 20% of the cells in the earlier report.

## 4. Spike mode and what it changes

**The rule** (`src/Profiling/SpikeGenes.h`; `--spike_genes auto|on|off`, default `auto`):
- Before a sample is profiled, protal counts each gene's records in the SAM: primary and supplementary records with a
  reference on a gene with unique k-mers (the records the profile counts), a pass over the lines on the sample's threads.
- **Per taxon** (`SpikesOf`): the background is its records on its genes other than its 3 with the most records per
  base, over those genes' length (every hittable gene of the taxon counts, hit or not). Those top genes are spikes, in
  order, while each holds at least 50 records and 50 times what the background gives a gene of its length, the taxon
  keeping 6 or more other genes with records; together they must hold at least half of the taxon's records.
- **Per sample** (`Judge`): the sample is in spike mode if at least 2 of its 10 taxa with the most records have spike
  genes and all spikes together hold 5% of its records or more. A library with carried-over products puts them on many
  species, the abundant ones among them; a relative's reads pile on a few genes of an absent congener, a thin taxon.
- **In spike mode** the spike genes' records are left out of the evidence as if the reads had not aligned (as suspect
  gene copies'). The log gives each sample's check where a taxon has spike genes, and `<profile>.spike_genes.tsv`
  lists the genes left out. `on` takes every sample's spikes, `off` counts nothing.

**The sample check on every SAM** (protal's log line per sample with `--spike_genes auto`; `spike_compare.sh` and
`spike_check_long_sims.sh`):

| Sample | Taxa with spike genes | Records on spikes | Of the 10 taxa with most records, with spikes | Spike mode |
|---|---|---|---|---|
| D6331 DNA (Dohlman) | 13 | 79.0% | 8 | yes: 20 genes, 80,244 records left out |
| D6331 HKU, Auckland, Treichel, hCom2 | 0 | 0 | 0 | no |
| hCom1a | 1 | 0.2% | 0 | no |
| ATCC HiFi, D6331 HiFi, D6331 ONT | 1-3 | 0.6-2.0% | 0 | no |
| strain mix ONT (*E. coli*, *K. pneumoniae*) | 1 | 0.5% | 1 | no |
| 30 simulated test samples (pe 150 and 250 bp, 10k-1M pairs) | 0-1 | 0-0.3% | 0 | no |
| 32 simulated long-read test samples (HiFi, ONT; 3 and 30 Mb) | 0 | 0 | 0 | no |

Outside the Dohlman sample every taxon with spike genes is absent: *Veillonella* congeners in D6331 HiFi and ONT,
*Cereibacter_A johrii* in ATCC HiFi, *E. fergusonii* in the ONT strain mix, *Hungatella effluvii* in hCom1a and
*Hydrogenobaculum* sp002878215 in one simulated sample (the `on` profiles' `spike_genes.tsv` against the truth). No
present species' gene was taken anywhere but in the Dohlman sample.

**How the thresholds were chosen** (`spike_rule_all.sh`, records per gene of every SAM):
- A fold of 10 also took absent relatives' piles in the clean libraries (x10-x88) and a few genes of present species in
  the simulated samples (x11-x82, 30-179 records); 50 with the share rule leaves both alone.
- The share rule (the spikes hold most of the taxon's records) is what separates artefacts from a relative's piles in
  short reads: the clean samples' absent relatives have at most 39% of their records on their top genes, the Dohlman
  species 56-97% on one to three.
- 100 records at least left *C. difficile*'s third pile (79 records, x168) in place, and the species, now with fewer
  fragments and still uneven coverage, fell below the knob; 50 takes it.
- In long-read samples the rule took absent *Veillonella* and *Escherichia* congeners' piles only, and in the D6331 HiFi
  sample two of them were then called (P 0.34 → 0.53 and 0.22 → 0.64). The first version kept long reads out; the
  sample check replaces that limit, since it tells the two patterns apart whatever the read type.
- Records count only where the profile counts them: on a gene with unique k-mers (others are skipped as no evidence). The
  prototype took every gene of `reference.map`, and in the D6331 sample one more pile of absent *V. dispar*, 12,683
  records on a gene without unique k-mers that the profile never counted.

**Calls in each mode** (`spike_compare.sh`, then `spike_calls.py <mode> off`; the ten mock samples and the 30 simulated
test samples with truth, the database's default forest models; TP/FP/FN):

| Mode | Dohlman | The other nine mock samples | Mock F1 | 30 simulated samples | Simulated F1 |
|---|---|---|---|---|---|
| `off` (protal as before) | 10/0/2 | 189/61/13 | 0.8397 | 671/15/45 | 0.9572 |
| `auto` (default) | **12/0/0** | 189/61/13, not one probability moved | **0.8445** | 671/15/45, none moved | 0.9572 |
| `on` | 12/0/0 | 189/64/13 | 0.8392 | 671/15/45 | 0.9572 |

- **`auto`:** in the Dohlman sample 24 taxa's probabilities move. Every present species scores the same or higher (the
  table of [section 1](#1-every-present-species-of-that-sample-is-depressed-not-only-akkermansia)), *Akkermansia* and
  *E. coli* cross the knob, and *C. difficile* stays at 0.609 (the first version, at 100 records, lost it).
- **`on`:** the same gain, and three absent congeners called where a relative's reads pile up: in the D6331 HiFi sample
  *Veillonella dispar* (0.344 → 0.531) and *V.* sp900549845 (0.219 → 0.641), in the ONT strain mix *Escherichia
  fergusonii* (0.312 → 0.547). Spike mode for every sample costs more than it gains.
- **Time:** the extra pass over the SAM took 0.1-0.4 s of the profiling step's 0.4-1.9 s per sample here (about a
  fifth; the alignment, not timed here, takes far longer). `off` skips it.

## 5. What this does not settle

- One contaminated library among the ten: the thresholds are tuned on what these samples show. The next r226 build
  trains and tests with `auto` on (the collector runs protal's defaults); its simulated samples have no such piles, so
  no sample there should be in spike mode, but the build's test set, not this report, has to show it.
- An absent congener of the spiked sample moved too: *Escherichia* sp004211955 from 0.359 to 0.469 under `auto`, below
  the knob but closer (not traced; the congener features compare a taxon with its relatives' evidence, which the
  removed pile changes).
- What the piles are. They look like amplicons of single-copy genes (one stretch of one gene per species, identity
  1.0); the study's methods were not checked.
- A long-read library with real artefacts of this kind is judged the same way, but none was at hand to check it.
