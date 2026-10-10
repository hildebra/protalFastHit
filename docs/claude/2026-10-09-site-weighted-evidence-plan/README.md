# Site-weighted evidence: per-column conservation from the family, in the features and the alignment (a plan)

**Status.** A plan, nothing implemented. Written 2026-10-09 (evening) on top of `6ae06fc` plus the uncommitted
changes of [the ancestry true-positive test](2026-10-09-ancestry-true-positive-test/README.md), section 10.

**The request (the user's).** Score a read's mismatches by how conserved the column is, with the weight growing as
the odds (99% conservation ten times 90%); estimate conservation at the family level for every genus (a reference
from at most 10 genera, aligned with WFA2), averaged over genera rather than across phyla; polarise two-member
genera the same way; consider the amino-acid (codon) version; say what it costs at r226 scale; and cost the
alternative of putting deep alleles into the index.

## 1. What protal has, and what the plan builds on

| today | where | reused for |
|---|---|---|
| the ancestry cache compares a hit copy with up to 17 congeners along shared 12-mers, once per run, and keeps per position the congeners' votes | `AncestrySites.h` (`Cache`, `Consensus`) | the run-time side: the per-copy weight view lives beside the sites |
| per-gene conservation factors from the full reference (Mash distances of the species' other copies) | `GeneConservation.h` | the fallback weight for a copy without a family estimate |
| the congener-gap build step aligns every copy against its congeners' copies with WFA2, gene by gene, all threads, ~140 MB of sequences per gene at r226, ~10 min of the training build | `Build.h` `WriteCongenerGaps`, `CongenerGapsTable.h` | the build side: the same loop, the same aligner, two more passes |
| codon positions are counted from the copy's start (`MismatchesByCodonPosition`): the copies are gene calls in frame | `Profiler.h` | the codon variant's frame |
| the profiler knows each species' genus and family (`GenusOf`, `FamilyOf`) | `Profiler.h` | grouping |
| the site shift of unsure reads counts every differing position alike | `AlignmentStrategy.h` `SettleBySites`, `StrainAlleles.h` `ShiftOf` | where the weights enter the alignment |
| r226: 16.6 M copies, ~1 kb each, 3.1 G index entries, ~35 GB resident | [memory audit](2026-10-03-memory-audit/README.md) | the budgets below |

## 2. The weight

For a column with conservation `p` (the share of comparisons that agree with the consensus), the odds that a
mismatch there is a real substitution of a close relative, against a sequencing error or a distant read's
difference, scale with `1 − p`. The per-mismatch weight is therefore

    w = −log(1 − p̂)        p̂ = (k + 1) / (n + 2)

with `k` agreeing of `n` comparisons (the pseudocounts cap what few sequences can claim: two congeners never
assert more than 0.75). Summed over a read's mismatches this is a log-likelihood ratio; the model sees the
product of the odds, which is the ten-fold step per nine the user asked for. Stored quantised to 4 bits in half
nats (0..7.5 nats, i.e. up to 99.9% conservation), which is all the sample sizes support.

Two time scales, both informative, estimated from the same alignments and stored as one weight each:
- **`w_within`**, from the species of each genus against their genus's reference copy, averaged over the genera
  of the family with equal weight per genus (the user's choice: no genus dominates, and no phylum-wide consensus is
  needed). This is the rate at which a strain's own mutations hit a column: what separates a close relative from
  errors and outsiders.
- **`w_among`**, from the genus references against the family reference (at most 10 genera). This is the rate at
  long range: the homoplasy risk of a derived site, and the reliability of a difference between twin species.

A phylum-level estimate is not needed for either and would be weak (the user's point): the alignments across 60-70%
divergence are unreliable, and the consensus means little. If it is ever wanted, the family references can be
chained upward (family references aligned to an order's), averaging `p̂` per column across families, at no
per-copy cost.

## 3. The build step: `column_weights` (per family and gene)

Added to the training and the final build after the congener gaps, in the same gene-by-gene loop (the gene's
copies are in memory then):

1. **The family reference.** For each family and gene, the genus references are the gene's copy of one species per
   genus, chosen by a hash of the genus name (deterministic), for up to 10 genera chosen the same way. One of them,
   again by hash, is the family reference: the coordinate system of the columns. The other genus references are
   aligned to it with WFA2 (as `StrainAllelesBuild.h Diff` does, both ends partly free, up to 0.3 divergence; a
   genus whose reference covers less than 0.7 of the family reference is dropped from the family estimate and
   keeps its own genus-level columns). Per column: the agreeing genus references of those compared → `w_among`.
2. **Within each genus.** Every species' copy of the genus is aligned to the genus reference (all of them up to
   24, else the nearest by sketch and a hashed sample, as the congener gaps do). Per genus-reference column: the
   share agreeing. The genus reference's alignment to the family reference maps the genus columns onto the family
   columns; the family column's `w_within` is the mean over its genera of their shares, with the pseudocounts
   applied per genus.
3. **Polarising two-member genera.** With the family reference and the other genus references in hand, the
   consensus of `AncestrySites.h` gets outgroups for free: where fewer than three congeners are compared, a
   position is the species' derived state when its congener's base matches the family consensus. This replaces
   the 0.7.9 fallback (which takes the nearest congener's differences, half of them the congener's own derived
   states) and removes the H3-fallback loss of the test (a novel species agreeing at 0.5 instead of near 0).
4. **The table `column_weights.bin`** (stored in `database.protal`, `--add_tables` for an existing database):
   per (family, gene) the family reference's length and two 4-bit weights per column; per copy its alignment to
   the family reference as run-length pairs (copy position, column; a handful per copy, more only across indels).
   Species without a family estimate (a family of one genus of one species) map to their genus's columns, or to
   nothing: then the gene's conservation factor stands in as a constant weight.

**Cost at r226** (16.6 M copies, ~4,500 families, 120 genes):
- within-genus alignments: one per copy, ≤0.1 divergence, ~0.3 ms each → ~1.5 CPU-hours, ~3 min on 32 threads
  (the congener gaps align up to 24 per copy and take ~10 min);
- family alignments: ≤10 per family and gene, up to 0.3 divergence, ~2-5 ms each → 5.4 M alignments, 3-8
  CPU-hours, ~10 min on 32 threads;
- storage: 4,500 × 120 × ~1 kb × 1 byte (two 4-bit weights) ≈ 540 MB, plus the copy mappings ~16.6 M × 8 B ≈
  130 MB: ~0.7 GB in the database and resident, 2% of the run's memory. Only the (family, gene) rows of the copies a
  run touches need loading; the rest can stay on disk if that ever matters.

Had the conservation been estimated per gene across all phyla instead, every copy would need a profile alignment
across up to 70% divergence (an HMM-style `hmmalign`, ~10 ms per copy): ~46 CPU-hours, 1.5 h on 32 threads, for a
consensus that means little across phyla. The per-genus-averaged family estimate is 20× cheaper and answers the
right question.

## 4. The run-time side

1. **The weight view.** The ancestry cache's entry for a (species, gene) copy gains the copy's column weights
   expanded along the copy (one byte per base, from the mapping): ~1 kB per cached copy, 200 MB at the cache's
   limit of 200,000 copies, shared by the run as the sites are.
2. **Counting per record** (`NoteRecord`, the CIGAR walk that counts the ancestry sites already): for every
   aligned column its `w_within`; for every mismatch (`X`) its `w_within` and `w_among`, split by what the column
   is (an ancestry site with the species' base, with the congeners', neither; a polymorphic site; none).
3. **Features** (a `weights` group, in the default set after the truth test and the r226 build say so):

   | feature | what | separates |
   |---|---|---|
   | `conserved_mismatch_ratio` | the mean `w_within` of the taxon's mismatches over the mean `w_within` of its aligned columns: 1 for mismatches spread as errors are, below 1 for a relative whose mutations avoid conserved columns | errors and reads from other genera (47% of r226's false positives) from any close relative |
   | `conserved_mismatch_excess` | the mismatches at columns of `w_within` ≥ 4.6 nats (99%) per kb aligned, less the sample's error rate | the same, as a rate the depth features can use |
   | `ancestry_agreement_weighted` | the ancestry agreement with each site weighted by `w_among` | a derived state at a column the genera never change is worth more than one at a hypervariable column |
   | `ancestry_fixed_agreement_weighted` | the fixed-site agreement alike | the same for the strain-versus-congener contrast |
   | `nonsynonymous_conserved_share` (codon variant, section 5) | of the mismatches at AA-conserved codon columns, the share that change the amino acid | the strongest outsider signal, once the frame is trusted |

   The models learn the rest. The existing per-gene conservation features stay.
4. **The alignment.** `ShiftOf` / `SettleBySites` give each differing position its `w_among` (quantised to the
   half-units they use) instead of one; a read unsure between a species and its twin is then settled by the
   columns that reliably differ. Later, the same weights can enter the MAPQ of all reads; that is a separate
   decision, since it changes every output. `--no_site_weights` keeps the current scores for comparison.
5. **Long reads** go through the same walk and gain more: a 10 kb read has hundreds of mismatches, and their
   distribution over the columns is a statistic of its own.

## 5. The codon (amino-acid) variant

How complicated: about a third more than the DNA version, as the same machinery runs on a second alphabet.
- **Frame.** The copies are gene calls, so codon 1 starts at position 0, as `MismatchesByCodonPosition` already
  assumes. The build checks each copy's translation for internal stops and a length divisible by 3; a copy that
  fails (frameshifted MAGs, partial calls) gets DNA weights only. A 64-entry codon table is added to C++
  (`insilico_strains.py` has one in Python).
- **Columns.** The family reference's codon columns are its base columns in threes. Each aligned copy's codons are
  read off the same base alignment; a codon split by a gap of length not divisible by 3 is skipped. Per codon
  column: the share of translations that agree → `w_aa` (one 4-bit weight per codon column, ~180 MB more at r226).
  Aligning the proteins themselves instead of the DNA would need a substitution matrix, which WFA2 does not score;
  within a family at ≥50% protein identity the DNA alignment gives the same codon columns, so it is not needed.
- **Per mismatch.** The read's codon is the reference codon with the read's base at the mismatched position (the
  read's own three bases when the record covers the whole codon). Synonymous: weighted as a DNA mismatch at a
  variable column. Non-synonymous: weighted by `w_aa` of the codon column. Within a species the strains' mutations
  are mostly synonymous (dN/dS 0.1-0.3), a distant read's are not, and errors do not care; so a non-synonymous
  mismatch at an AA-conserved column is the strongest single piece of evidence the read is not a close relative.
- **Cost.** 2-3 days after the DNA version; the same build loop, the same per-record walk, one more table column.

## 6. Order of work and the checks before each step

| phase | days | what | the check that decides the next phase |
|---|---|---|---|
| 0 | 1-2 | **no protal change.** A Python prototype on the r226 v21 error reads (`model_logs/error_reads`, the pool's marker genes): family references by the 12-mer pairing the oracle uses (no WFA2 in Python), `w_within` and `w_among` per column, the mismatch-weight ratio and the weighted ancestry agreement per read, dN/dS of the mismatches | FP reads against TP reads: the AUC of `conserved_mismatch_ratio` (expect the cross-genus FP far above 0.5); the weighted against the plain ancestry agreement on the strain-versus-congener errors; the share of non-synonymous mismatches at AA-conserved columns. A signal that is not there on real reads is not built |
| 1 | 3-4 | the build step and the table (section 3), the loader, `--add_tables`, unit tests on hand-made families, the polarised two-member consensus | the mini world's columns against its true rates; the truth test's H3 fallback with outgroups |
| 2 | 2-3 | the per-record counting, the features, the oracle port, a gamma-rate option in `simulate_ancestry_world.py` (per-site rate classes, so the world has among-site variation at all) | protal = oracle on every row; the features' AUC in the truth test's hard strata; `conserved_mismatch_ratio` on reads of another genus planted into a sample |
| 3 | 2 | the weighted site shift, `--no_site_weights` | L1 of the truth test: moved reads by origin, the novel species' reads drawn onto the target (today 5-7.5%) |
| 4 | 2-3 | the codon variant | the same, with the nonsynonymous feature |
| 5 | one build | r226: training with the new groups, `--features auto` | the real-strain misses and the cross-genus FP against v21 |

About two working weeks to phase 4, with phase 0 as the gate. Nothing user-facing changes until phase 5; the
website then needs a line on the new table and the option.

## 7. Risks and limits

- **Reference quality.** A frameshifted or chimeric copy makes columns look variable (a lower weight, harmless)
  and, if it is the family reference, misaligns a whole family: choose the family reference among genus
  references whose translation is clean and whose alignments cover the others, not a random one when that fails.
- **Alignment reach.** At 0.3 divergence WFA2 with protal's scores still aligns marker genes; beyond that the
  genus is dropped from the family estimate and keeps its own columns. Families of one genus get the genus-level
  columns; genera of one species get the gene factor.
- **The alleles' sites.** A polymorphic site of the species (strain alleles) must not count as a conserved-column
  mismatch for the species' own strain: the per-record split by site class (section 4.2) keeps the two apart.
- **Determinism.** Genera and references chosen by hash, loops in sorted order, as everywhere in the build.
- **Memory.** ~0.7-0.9 GB at r226 for both variants, inside the budget.

## 8. The alternative asked about: deep alleles in the index

The misses that no feature reaches are the strains whose reads never seed on their species (one or two fragments,
reads lost at seeding to a congener). Indexing the stored alleles' k-mers would let those reads find the species.

- **How.** After the alleles are chosen, the edited copies are scanned for syncmers, and the k-mers not already
  in the index for that copy are inserted with the same (taxid, gene, position); both index passes see them. The
  alignment then runs against the representative's copy as now (the allele scores explain the differences), so
  nothing downstream changes but the candidates. The uniqueness tables (`unique_kmers.tsv`, the `lu`/`su` ranks)
  are computed against the full reference, which contains the allele genomes, so they stay consistent.
- **Memory.** r226 has ~187 entries per copy (one syncmer per ~5-6 bases). A substitution changes about `k`
  k-mers, of which about a fifth are syncmers: ~3-4 new entries per edit. With ~40 edits per copy with alleles
  (the test world's 15 per allele, 4 alleles) that is ~150 entries, +80% for those copies; if a third of r226's
  copies have alleles, +25-30% of entries: **+4-5 GB** on the packed layout (~16 GB of entries today). Restricting
  it to the alleles' edits ≥ 1% from the representative (the deep lineages, which are the misses) cuts that to
  about +1-2 GB.
- **Work.** 3-5 days (the index builder's two passes, the packed layout's tolerance of several entries per
  position, a `--index_alleles` option, tests), plus one cluster run for the memory and the accuracy.
- **Gain.** The FN anatomy's ceiling was +0.009 of F1 if every missed strain's reads found their species; a
  realistic third of that, since many of those strains still lose to a congener at the MAPQ filter. Worth doing
  after the site weights, whose phase 0 costs a day and targets the larger error class.

## 9. Implementation (2026-10-10)

The user chose the C++ implementation over the Python prototype of phase 0, since the r226 reference genes are not
on this machine (the share archive holds the error reads but no database). Phases 1, 2 and 4 of section 6 are in;
the weighted site shift (phase 3) is not. Uncommitted, on top of the changes of
[the ancestry test](2026-10-09-ancestry-true-positive-test/README.md) section 10; another session was editing the
same checkout at the time (the deep alleles in the index, `--index_alleles`), whose hunks these leave alone.

| part | where | what |
|---|---|---|
| the table | `src/SequenceUtils/ColumnWeights.h` | `column_weights.tsv`: per family and gene an `F` line (the family reference's species, the voting genera, per column the `within` and `among` codes as hex digits, per codon the `aa` code, the consensus bases), per copy a `C` line (its runs on the columns); `Table` (Find, Expand, Read/Write), `Cache` (a copy's columns expanded once per run), `Count` (a record's mismatches by column, the codon classification), `CodeOf` (round(2 w), w = −log(1 − p̂), p̂ = (k + 1)/(n + 2); 0 no estimate, 9 = 4.5 nats "conserved") |
| the build | `src/SequenceUtils/ColumnWeightsBuild.h`, `Build.h` `WriteColumnWeights` (after the strain alleles) | per gene the copies in memory, per family (one OpenMP task) `BuildFamily`: genus references and voters by hash, `MapTo` (WFA2, protal's scores, ends partly free, up to 0.4 among genera and 0.2 within), the votes per column and per codon, the mean of the genera's shares, every copy's runs; `--column_weights N` (genera that vote, default 10; 0: no table); the file packed into the database, accepted by `--add_tables` |
| the run | `GenomeLoader.h` (`ColumnWeightsOf`, `SetColumnWeights`), `RunProtal.h` (`LoadColumnWeights`, beside the other tables), `Profiler.h` (`RecordEvidence::cw_*`, `NoteRecord`, the `Taxon` accessors, `TaxonFeatures`) | five features in the `weights` group of the default set (91): `conserved_mismatch_ratio`, `conserved_mismatch_rate`, `ancestry_agreement_weighted`, `nonsynonymous_share`, `nonsynonymous_conserved_rate`; `column_weight_coverage` in no group; `--no_site_weights` leaves the table unread |
| the polarised consensus | `AncestrySites.h` `Consensus(..., outgroup)`, `Cache::Get(..., outgroup)` | with fewer than three congeners compared, the nearest congener's difference is a site only where the congener's base is the family's consensus |
| the scripts | `model_features.py` (`WEIGHT_FEATURES`, the default set, the named sets), `ancestry_oracle.py` (the table read and expanded, `count_columns`, the features), `ancestry_truth_test.py` (the `weights` lists, the model sets `weights`, `without weights`, `default`), `simulate_ancestry_world.py` (`--site_rates ALPHA`: a gamma rate multiplier per site, the same on every branch; `--genera_per_family N`: families of several test genera) | |
| tests | `tests/test_ColumnWeights.cpp` (7: the scale, the map and runs, a hand-made family, the round trip and expansion, a record's counts with codons, the features through the profiler, the polarised consensus), `test_model_pmml.py` | |
| docs | `features.md` (the group's section), `databases.md` (the table, the build step), `main.cpp` | |

The AA conservation needs no protein alignment: the copies are gene calls in frame, so a codon column is three base
columns of the family reference, and the amino acids are read off the DNA alignment (a codon split by a gap is
skipped). The memory at r226 is about two bytes per column of the family references (~1.1 GB) plus the copies'
runs; the plan's 0.7 GB would need the aa codes packed per codon, left for later.

**Two faults the first runs found, both fixed before the runs reported below:**
- A `std::string_view` of a temporary: the profiler took the copy's sequence as `GetGeneOMP(..).Sequence().View()`
  without keeping the sequence object, so the codon classification read freed memory; one mismatch in a few hundred
  came out wrong (protal's `nonsynonymous_share` 192/262 against the oracle's 193/263). Found by running protal's
  `Count` on the oracle's inputs in a harness, which gave the oracle's number. The ancestry cache keeps the object;
  the profiler now does too.
- The consensus base must be an outgroup. With one genus per family (the previous worlds) the "family consensus"
  was a single congener's copy, or the species' own, and polarising the two-congener genera by it removed their
  sites (the H3 fallback had no rows, H4 at b = 0 fell from 0.906 to 0.805). A column has a consensus base only
  from a majority of at least two among three or more genus references (`kMinConsensusGenera`, `kMinConsensusVotes`).

### The truth test with the weights (2026-10-10, WSL, 4 cores, ~15 min a run)

The user asked how effective the weights are on the previously used simulation. Two runs of
`ancestry_truth_test.py` with this binary, both worlds and reads byte-identical to the earlier runs where the
options are the same ([atp_cw_base/](atp_cw_base/summary.md), [atp_rates/](atp_rates/summary.md)):

| run | world | `--world-args` |
|---|---|---|
| `atp_cw_base` | the previous worlds: every site evolves at one rate, a family per test genus | `--long_samples 32` |
| `atp_rates` | the same design with a gamma(0.5) rate multiplier per site, the same on every branch, and three test genera per family | `--long_samples 32 --site_rates 0.5 --genera_per_family 3` |

**The implementation is exact.** protal's 91 features equal the oracle's recomputation on every target row of
both worlds, in the default, `--no_allele_scores` and shuffled runs, with error-free and paired-end reads (H1b/H2b
PASS), the five weights features included. The tables built for every copy, no alignment failed (2,400
family-gene rows in the rates world, 5,280 in the base world, 24,480 copies mapped).

**The models gain nothing from the weights in either world** (paired-end, AUC on the target rows; the hard stratum
is Q_deep against N_0.8/0.95; F1 over all test rows at the knob):

| model | base worlds | rates worlds |
|---|---|---|
| base | 0.858 / hard 0.671 / F1 0.853 | 0.856 / 0.669 / 0.848 |
| base + weights | 0.860 / 0.670 / 0.859 | 0.857 / 0.664 / 0.847 |
| all but weights (the previous default) | 0.920 / 0.797 / 0.909 | 0.925 / 0.802 / 0.905 |
| default (with weights) | 0.921 / 0.803 / 0.907 | 0.925 / 0.798 / 0.905 |
| HiFi base → default | 0.885 → 0.941 | 0.888 → 0.948 |

Over identity and depth alone the group adds +0.044 [0.023, 0.065] of AUC in the base worlds and +0.060 [0.044,
0.073] in the rates worlds, about what the ancestry group adds (+0.043, +0.056), but all four groups together add
+0.119 against +0.116 for three: the base set already carries what the weights see here.

**Why these worlds cannot show the weights' point.**
- Their genera have at most seven species and their families three genera, so no column reaches "conserved": the
  largest within weight is 2 nats (code 4 of 15) and the largest among weight 1.5 nats; the family references of
  the base worlds are single copies. `conserved_mismatch_rate` and `nonsynonymous_conserved_rate` are therefore 0
  on every row (AUC 0.500), and `ancestry_agreement_weighted` is the plain agreement (AUC 0.840 against 0.838).
  At r226 a genus has up to thousands of species and a family dozens of genera.
- The samples hold no reads from outside the database's genera. The weights target the cross-genus reads of
  r226's false positives (47% of them), whose fast columns are saturated; a strain and a novel congener at the
  same identity have the same young mutations at the same columns, which is what these worlds compare.
- `conserved_mismatch_ratio` still separates at 0.61-0.66 (higher for strains) and `nonsynonymous_share` at 0.50:
  with every site alike, synonymous and non-synonymous changes are as frequent in a strain as in a novel species.

**The polarised sites work where the family has outgroups.** In the rates worlds (three genera per family) the
two-congener genera's sites are the species' own derived states at 0.84-0.88 instead of 0.50 (the congener's own
derived states fell from ~3,200 to ~320-550 sites per genus; `l0_sites.tsv`), and `ancestry_agreement` separates
deep strains from novel species that left the stem early at 0.967 (b = 0) and 0.867 (b = 0.5) against 0.906 and
0.815 in the base worlds, whose single-genus families give no consensus. The price is recall: the sites the
consensus does not cover are dropped (stem recall 0.70 against 0.96), which the models did not miss. H3's
fallback prediction (the congener's own states enter) now fails in the right direction.

**What a test of the weights' point needs** (not done): worlds with 30-100 species per genus and 10 genera per
family, so that columns reach 99% conservation; samples with reads of held-out genera landing on the database's
species at 80-90% identity; and sequencing errors on conserved columns. The r226 build remains the real test,
with its error reads' `conserved_mismatch_ratio` the first thing to look at.

**Tests.** Unit 488 (3 skipped as always; `test_ColumnWeights.cpp` 7), e2e 149 on a mini database built by this
binary (its build line: 237 family-gene rows of 3 species), pipeline 6 (a GTDB-like build with the weights), scripts
41 + 8, ancestry world 19 (with the quick end-to-end run and the oracle's equality): all OK on WSL, 4 cores.

Committed as `f566cd0` (my hunks only, from a temporary index: another session's index-alleles work shares the
checkout; the committed tree builds and passes its unit tests on its own).

### Phase 3: the weighted site shift in the alignment (2026-10-10, after the commit)

`StrainAlleles.h ShiftOf` takes the copy's expanded columns (`column_weights::Columns`): of a candidate's
differences that no allele explains and no polymorphic site covers, a substitution at a column with an among code
below `kDiscountNone` (6, 3 nats, ~95% of the genus references agree) loses half a difference, one at `kDiscountFull`
(2, 1 nat) or below the whole difference (`DiscountHalves`, `SiteShift::discounted`). The alignment handler
(`AlignmentStrategy.h ScoreAlleles`) looks a copy's columns up once per handler and computes the shift for copies
without alleles too; `SettleBySites` is unchanged, so only the reads protal is unsure about (another species'
candidate within 3 mismatches) take it, and `--no_allele_scores` or `--no_site_weights` turns it off. Documented in
`running.md` and `features.md`; tests: the shift with columns in `StrainAlleles.ThePolymorphicSitesOfACopy`, and
`ColumnWeights.TheShiftDiscountsDifferencesAtVariableColumns` (two species four differences from a read each, the
one whose differences sit at hypervariable columns wins only with the table).

**The first version was wrong, and the rates worlds showed it at once**
([atp_shift_full_discount/](atp_shift_full_discount/summary.md)): it took a whole difference off at a hypervariable
column (among code 2 or less) and half at a middling one. In these worlds every column is below 3 nats, so every
candidate's unexplained differences were discounted alike and mostly in full, the candidates of an unsure read
became ties, and the strains' reads the shifts moved went to their own species at 0.505 instead of 0.999 (H7
settling FAIL); strain reads on the target fell from 0.946 to 0.934 without a twin and from 0.930 to 0.896 with one,
the novel species' reads from 0.866 to 0.811; the paired-end default model fell from 0.925 to 0.914 on the target
rows (hard stratum 0.798 → 0.778), HiFi unchanged (no shift for long reads). The lesson is in the code's comment:
the discount can never exceed half a difference, so that the count of differences still orders the candidates and
the columns' reliability decides only among those within the unsure margin.

