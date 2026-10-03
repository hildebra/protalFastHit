# What amplicon denoising suggests for protal's species calls

- **Date**: 2026-10-02.
- **Question** (the user's): false calls that grow with depth are a known problem. Amplicon denoisers
  dealt with it: AmpliconNoise (Quince et al. 2011, [PMC3045300](https://pmc.ncbi.nlm.nih.gov/articles/PMC3045300/)),
  then DADA2 (Callahan et al. 2016, [PMC4927377](https://pmc.ncbi.nlm.nih.gov/articles/PMC4927377/)) and
  UPARSE/UNOISE (Edgar 2013, 2016; [drive5.com/uparse](https://drive5.com/uparse/)). What did they do, and does any
  of it apply to protal's species detection? (The dada2.org link the user gave is the patient foundation for
  ADA2 deficiency, a disease, so the DADA2 paper was read instead. UNOISE's parameters are from the usearch
  manual's `unoise3` page.)
- **Code**: protal at `6d2e241` (`audit-fixes`), read to see how a call is made. The trainer is
  `scripts/random_forest_cmdline.py` of `6d2e241` with two more feature sets
  ([`model_features.patch`](model_features.patch)). No protal code was changed.
- **Data**: the training and test tables of the r226 v3 build (job 23887614,
  [report](../2026-10-02-r226-v3-training/README.md)), in `local/protal0.7.3_r226_v3/{training,test}/training_data*.tsv`
  (git-ignored).
  - Training: 198 samples, 1,000 to 10M read pairs, abundances Poisson-lognormal with σ 1.3.
  - Test: 78 samples, 500 to 5M read pairs, σ 2.0.
  - Species are drawn uniformly from GTDB r226 (`--congeners 0`). 20% of species and 32 clades are held out of the
    database.
- **Run**: [`run.sh`](run.sh) in WSL (protal-db-build env: Python 3.14, scikit-learn 1.9.1), about a minute per model.
  It writes the outputs here: [`compare_output.txt`](compare_output.txt), [`congeners_output.txt`](congeners_output.txt),
  [`ambiguity_output.txt`](ambiguity_output.txt), [`describe_species_held_out.txt`](describe_species_held_out.txt) and
  [`describe_test.txt`](describe_test.txt). The tables and models are in WSL `~/denoise_r226`.

## What the three methods do

They share protal's problem. Each read carries about the same chance of an error, so the number of distinct
erroneous sequences grows with the reads. Any rule of the form "a unit exists if it has reads" then inflates with
depth. protal shows the same shape: at r226 v1 the absent taxa per sample grew as depth^0.84.

- **AmpliconNoise** fits a mixture model by EM.
  - Reads are spread around L true sequences by an error-model distance (scale σ). Each sequence has a weight
    τ_j, its frequency.
  - A read is assigned by its posterior over the sequences, and τ is part of that posterior. An abundant
    sequence therefore explains its rare neighbours' reads.
  - Components that lose all their weight disappear.
  - The error model (base transition probabilities) was measured on mock communities.
  - Its chimera step, Perseus, uses abundance: a chimera's parents must be at least as abundant as the
    chimera, so only more abundant sequences are searched. The final call is a logistic regression trained on
    mock communities.
- **DADA2** works in the following steps.
  - It learns λ_ji, the rate at which sample sequence j is read as i, from base qualities. The rates are
    estimated from the data, alternating with the inference until both are consistent.
  - A sequence i next to a more abundant j becomes its own partition only if a Poisson test fails. The test
    asks: given n_j·λ_ji expected copies, how likely are a_i or more? The p-value must fall below Ω_A after a
    Bonferroni correction over the sequences tested.
  - Singletons are never called: the test is conditional on the sequence being seen, so their p-value is 1.
  - Two-parent chimeras (bimeras) are checked against more abundant parents.
- **UPARSE / UNOISE** go through the sequences in decreasing abundance.
  - A sequence becomes a new centroid only if no more abundant one explains it. In UPARSE that means it is
    more than 3% from every centroid.
  - In UNOISE, a sequence is an error of its neighbour if its abundance skew to that neighbour is below
    β(d) = 1/2^(αd+1), with α = 2 and d the number of differences. A one-difference neighbour must have under
    1/8 of the abundance, a two-difference neighbour under 1/32.
  - Sequences with fewer than 8 reads (minsize) are not considered.

**Common principle**: each candidate is judged against the more abundant candidates that could have produced
it, at an expected noise rate that falls with the distance between them.

- The bar scales with the parent's abundance, not with the sample's depth.
- The other half of the rule: a candidate with nothing nearby that could produce it is believed on little
  evidence (two reads in DADA2's mock communities).

## How protal decides now

Paths are in `src/`.

- **Read assignment.** Each read's primary record goes to the best-scoring taxon.
  - MAPQ = 1 + 40·(1 − s2/s1)·log10 s1 against the second best, whatever its taxon (`Alignment/AlignmentUtils.h:771`).
  - Exact ties get MAPQ 0 and fall to the MAPQ ≥ 4 filter (`Profiling/Profiler.h:3211`).
  - There is no EM and no prior from abundance.
- **Features.** The forest sees each taxon's own features. Only two concern other taxa: `congener_fit_share` and
  `other_genus_fit_share`. They say whether another taxon fits a read within one edit (the ZA tag), not whether
  that taxon is in the sample or how abundant it is.
- **Nothing compares a taxon with a more abundant relative.** The build's `gene_congeners.tsv` is a report that
  queries do not read.
- **Depth knobs.** They raise the threshold with the sample's total fragments. That stands in for "the relatives
  are abundant", but it applies equally to taxa with no relative in the sample.

## Where the analogy holds, and where it breaks

- **It holds** for false positives next to a species the database has, that is reads spilled over from a present
  relative. On species held out at knob 0.5 these are 208 of the 582 paired-end false positives (36%); on the test
  set 114 of 281 (41%).
- **It breaks** for false positives next to a species the database lacks (64% and 59%). Their "parent" is not a
  taxon protal can see, whereas amplicon denoisers are reference-free and always see the parent.
  - The nearest lesson is AmpliconNoise's: a read cloud centred away from the reference is another sequence.
    `excess_*` already measures that.
  - Such a taxon could be reported as an unknown member of its genus instead of as the species.
  - The limit at the species boundary stays: 79% of misses at r226 v1 were divergent strains.
- **The noise is biological, not sequencing error.** Two groups both sit beside a congener with ≥10× their
  fragments ([`ambiguity_output.txt`](ambiguity_output.txt), test set):

  | group | count | median identity | `excess_median` |
  |---|---:|---:|---:|
  | absent taxa (spill-over) | 7,383 | 0.94 | 0.045 |
  | present minor congeners | 56 | 0.985 | 0.000 |

  - So λ here is a cross-mapping rate per gene and per pair of genomes, not a quality-driven error rate. Spill-over
    lands on genes the parent's reference lacks, or where the parent's strain differs.
  - Per-read ambiguity does not separate the two groups: `congener_fit_share` has a median of 0.23 to 0.50 in both.
- **No dereplication.** Shotgun reads are not copies of one locus. The unit is fragments per taxon, or per gene.

## Experiment: abundance relative to the sample's relatives as features

All features are computed from the rows of the same sample, which protal has at run time (every taxon's fragments,
plus the taxonomy) ([`add_relative_features.py`](add_relative_features.py)):

| feature | definition |
|---|---|
| `genus_skew` | log10((n+1)/(m+1)), where m is the most fragments of another species of the genus (UNOISE's skew) |
| `family_skew` | the same against the largest species of another genus of the family |
| `genus_share` | the taxon's share of its genus's fragments |
| `genus_spill` | log10((n+0.5)/(10⁻³·Σ congeners + 10⁻⁴·Σ the family's other genera + 0.5)): DADA2's n_j·λ, with λ by rank only |

The models compared:

- **base**: normalized+adjacency, today's default. It reproduces the v3 retrain exactly: test F1 0.9542 at 0.5,
  0.9618 at the knob curve.
- **rel**: base plus the four features above.
- **spill**: base plus `spill_divergent` = max(0, −`genus_skew`) · max(0, `excess_median`). This counts the skew only
  as far as the reads differ from the reference beyond their errors ([`add_spill.py`](add_spill.py)).

### F1 at the knob of 0.5 and at the model's own knob curve

| read type | evaluated on | base 0.5 | rel 0.5 | spill 0.5 | base curve | rel curve | spill curve |
|---|---|---:|---:|---:|---:|---:|---:|
| pe | species held out | 0.9593 | **0.9693** | 0.9622 | 0.9691 | **0.9724** | 0.9694 |
| pe | test set | 0.9542 | **0.9632** | 0.9568 | 0.9618 | 0.9617 | 0.9609 |
| se | species held out | 0.9540 | **0.9659** | 0.9586 | 0.9644 | **0.9696** | 0.9656 |
| se | test set | 0.9471 | **0.9610** | 0.9543 | 0.9582 | 0.9608 | 0.9602 |
| pb | test set (22 samples) | 0.9736 | 0.9718 | 0.9732 | 0.9745 | 0.9758 | 0.9737 |
| ont | test set (22 samples) | 0.9701 | **0.9763** | 0.9722 | 0.9704 | 0.9698 | 0.9716 |

- **At a fixed knob of 0.5, rel matches or beats base with its depth knobs** (pe 0.9632 against 0.9618, se 0.9610
  against 0.9582, ont 0.9763 against 0.9704).
  - rel's knob curve is flatter: 0.29-0.83, against base's 0.05-0.93.
  - The best single test threshold moves from 0.68 to 0.51.
  - The parent's abundance does the depth knobs' job. It does not need to extrapolate past the deepest trained
    sample, which the r226 v1 evaluation named as the main risk on real samples.
- **False positives next to a species the database has fall by half**: pe test 114 → 51, se 125 → 48.
- **False positives next to a species the database lacks do not change** (pe test 167 → 161), as the analogy
  predicts.
- **By depth** (pe test, at 0.5):
  - At 1M pairs, false positives fall from 109 to 67.
  - At 2,000 pairs, misses fall from 19 to 8.

  The second effect is the other half of DADA2's rule. A low-abundance species with no larger congener is believed:
  misses in that group fall from 142 to 90.

### The catch: present species beside a more abundant present congener

Present species are counted by how many times more fragments the largest other present congener has, and how many
are missed at 0.5 ([`congeners_output.txt`](congeners_output.txt), test set):

| read type | congener has | present | missed base | missed rel | missed spill |
|---|---|---:|---:|---:|---:|
| pe | no larger congener | 4,671 | 142 | 90 | 118 |
| pe | 1-10× | 201 | 22 | 27 | 26 |
| pe | ≥10× | 52 | 11 | **36** | 19 |
| se | no larger congener | 4,571 | 172 | 94 | 137 |
| se | 1-10× | 190 | 17 | 32 | 26 |
| se | ≥10× | 49 | 14 | **36** | 19 |

- **rel misses 69-73% of minor congeners.**
- **The simulation rarely has them.** Species are drawn uniformly from GTDB, so congeners hardly ever share a
  sample. The paired-end training set has 35 present species beside a ≥10× congener, against about 20,800 absent
  taxa in that position. The forest learned "an abundant congener means absent", and F1 barely sees the cost.
- **Real samples have many.** In gut samples, species of *Bacteroides*, *Bifidobacterium* and other genera
  routinely share a sample at 10-1000× differences. There rel would lose real species that today's model finds.
- **Gating by divergence (spill) is no fix.** It cuts the cost by more than half (19 misses instead of 36) but keeps
  little of the gain. The forest does have the separating evidence (identity 0.985 against 0.94), but it cannot
  learn it from 35 cases.

**Singletons beside an abundant congener** ([`describe_species_held_out.txt`](describe_species_held_out.txt),
[`describe_test.txt`](describe_test.txt), base at 0.5):

| data | congener has | single-fragment calls | true | false |
|---|---|---:|---:|---:|
| species held out | ≥100 fragments | 40 | 0 | 40 |
| species held out | ≥10 fragments | 102 | 4 | 98 |
| test set | ≥100 fragments | 19 | 0 | 19 |

## What would apply, most useful first

1. **Train and test with congeners in the samples.**
   - This is DADA2's "Extreme" mock-community lesson: 27 strains over five orders of magnitude, some one nucleotide
     apart. UNOISE's figure has the same case: a low-abundance correct sequence wrongly discarded.
   - The option exists: `--congeners N` in `collect_training_data.py` and `build_gtdb_database.py` puts one genus of
     N species in every sample of a design point. It is 0 by default and was not used at r226.
   - Better still: several genera with 2-5 species each in every sample, at σ 2.
   - Until then no abundance-context feature can be judged. Even base's own 21% miss rate of minor congeners rests
     on 52 cases.
2. **Then the relatives features in protal.**
   - Implementation: one pass over a sample's taxa after their features are computed, before the model.
   - On these tables they add 0.009 (pe), 0.014 (se) and 0.006 (ONT) of test F1 at a fixed knob, and largely remove
     the need for depth knobs.
   - λ should fall with the distance between the two references, as in DADA2 and UNOISE, rather than be one value
     per rank. The build already compares congeners' gene copies (`GeneConservation.h`, for
     `gene_conservation.tsv` and `gene_congeners.tsv`). Keeping a distance per pair of congeners would give UNOISE's
     d.
   - Retest the minor congeners on the congener-rich design before making it the default.
3. **Abundance-weighted read assignment.**
   - This is AmpliconNoise's EM; for metagenomes, Pathoscope (Francis et al. 2013) and EMU (Curry et al. 2022) did
     the same.
   - The posterior of a read on taxon k is ∝ τ_k·L_k over its ZA alternatives, with τ from a first pass.
   - It helps with near-ties: reads within a few edits of two taxa.
   - It does not reach most of the spill-over seen here, at 0.94 identity. There the parent's reference lacks the
     read's gene or is further away, and only the taxon-level comparison (item 2) sees it.
   - It is more work: a second pass over the records, and ZA lists at most 4 alternatives within 5 edits.
4. **Singletons beside an abundant congener.**
   - DADA2 and UNOISE never call singletons. protal must, for shallow samples: 71% of present single-read taxa were
     called on the V2 test.
   - "One fragment and a congener with ≥100 fragments → absent" removed 40 false calls and no true ones on species
     held out (19 and 0 on the test set).
   - It is subsumed by item 2 if that is adopted, and is a cheap fallback otherwise. The minor-congener caveat
     applies here too, but with one fragment such a species is a coin flip anyway.
5. **Multiplicity instead of depth bins.**
   - This is DADA2's Bonferroni correction over the sequences it tests. In protal, the number of candidate taxa per
     sample is what grows with depth.
   - A per-sample threshold could keep the expected number of false calls under a target, for any depth: Σ(1 − p)
     over the calls, with calibrated p.
   - It is less needed if item 2 flattens the knob curve.
6. **A Perseus/bimera analogue.**
   - Flag a taxon whose hit genes are each covered more deeply by a different abundant relative.
   - This is the per-gene form of item 2. Later.
7. **Pooling across samples** (DADA2's pseudo-pooling), with prevalence across a run's samples as a prior.
   - Low priority: it makes one sample's calls depend on the others.

**Not applicable**: clustering at a fixed radius (OTUs), reference-free inference, dereplication, and error models
from quality scores (`excess_*` already uses the qualities).

No protal behaviour changed, so the website needs no update.
