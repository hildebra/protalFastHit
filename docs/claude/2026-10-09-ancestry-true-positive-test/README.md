# A true-positive test for the ancestry sites, strain alleles and polymorphic sites: a simulated world where the signal is known

**Status.** Designed, then implemented the same day (section 9, with the user's two additions: sensitivity by how
long ago a lineage diverged, and recombination with congeners), and run in full locally in 17 minutes (section 9,
[full_run/](full_run/summary.md)). The cluster job script is `scripts/mini_db/run_ancestry_truth_test.sh`, for larger
worlds.

**Data and commit.** The design reads the code at `6ae06fc` (protal 0.8.0, `a82df8f`):
- `AncestrySites.h`, for the sites and their consensus;
- `StrainAlleles.h` and `StrainAllelesBuild.h`, for the table, `BestAllele`, `Polymorphism`, `CountSites` and `ShiftOf`;
- `Profiler.h`, for `NoteRecord`, `NoteAlleles` and the feature accessors;
- `AlignmentStrategy.h`, for `ScoreAlleles` and `SettleBySites`;
- `scripts/mini_db/simulate_gtdb_release.py`, `simulate_reads.py`, `simulate_metagenomes` (`--from_manifest`,
  `--long_samples`) and `machine_learning_cmdline.py` (`--test-file`).

**The request (the user's, 2026-10-09).** Design a simulation that tests whether protal's ancestry, polymorphic-site and
strain-allele features actually work: a true-positive test.

**Short answer.**
- **A world grown on explicit trees.** The test world evolves the marker genes down explicit trees: a genus tree of the
  database's species, a genealogy inside each target species, and novel species attached to the target's stem at a
  chosen fraction `b`. Every site's history is then known.
- **What goes into the database.** The representative and a few "allele genomes" per species. The strains used for
  samples and the novel species are never in the release, so the honest split holds by construction.
- **Each test row is one taxon's row in one sample.** The target species S has either one of its own unseen strains
  in the sample (a true positive) or a novel species whose reads land on S (a false positive). The rows of the two
  classes share the database side (alleles, sites, cluster size), and their identities are matched.
- **Four levels, each against a known answer:**
  - the tables against the history;
  - each read's assignment against its origin;
  - each taxon's features against an oracle and against the analytic expectation;
  - the calls of models trained on a second world, with and without the feature groups.
- **Pre-registered predictions, including where the features should fail.** The theory below predicts a trap. A
  novel species carries the ancestral base at every site where the representative's own lineage is derived. The
  alleles carry that base there too, so two allele features and two polymorphic features count the novel species'
  reads as explained there. Only agreement at the fixed sites separates a late-branching novel twin from a deep
  strain at the same identity. The test measures whether the features as built recover that.
- **Cost.** About a day of scripting (Python, no protal change needed) and under an hour of compute at mini scale.

## 1. Why the existing worlds cannot answer it

**The mini world** (`simulate_gtdb_release.py`) evolves every species of a genus from the genus node, and every
genome of a species from the species node. Both are stars.

- **A held-out congener is trivially separated.** It shares none of S's derived states, so at S's ancestry sites it
  carries the congeners' base nearly everywhere: agreement ≈ 0. The real false positives at r226 are close relatives,
  and those agree at 0.68 (v17).
- **The strains share nothing but the representative's own mutations.** A strain's differences from the
  representative are the representative's private mutations (where every other genome has the ancestral base) plus
  its own. No other genome shares its own mutations. So no strain is "known" beyond what any outsider also shares.
- **The in-silico strains of the GTDB build** (`insilico_strains.py`) are mutated copies of one-genome species, and
  have no alleles at all.

**r226** has the real structure but no truth about where each site's history lies. The v18/v19 analyses inferred it
from the full reference, and the first inference was circular (v18). Even the honest number (fixed-site AUC
0.58-0.68, v19) does not say whether protal's code reaches what the data allow.

**What is needed:** a world where trees, not stars, set which genomes share which states, and where the evaluation
can compute what a perfect implementation would give.

## 2. What each feature should see: the site classes

**Notation.** For one marker gene, on S's representative copy, with the consensus over at least 3 DB congeners
(`kMinCongeners`):
- **L:** substitutions on S's stem, from its split with its nearest DB relative down to S's crown. Every S genome
  carries them.
- **b:** a novel species N_b leaves the stem at fraction `b`. It shares `bL` of these sites and lacks `(1-b)L`.
- **R:** substitutions on the representative's own lineage inside S, from the crown down to the representative.
- **R_Q ⊆ R:** those below MRCA(Q, rep). A strain Q carries the ancestral base there, which is also the congeners'
  base.
- **K_Q:** Q's own substitutions that an allele genome shares (above MRCA(Q, A)).
- **U_Q:** Q's private substitutions.
- **P_N:** N's own branch.

**How the code sees them:**
- protal's ancestry sites are `L + R`: the representative differs from every congener there.
- The polymorphic sites are the allele genomes' own lineages plus the part of `R` below MRCA(A, rep) for some allele
  genome A.

| | strain Q of S | novel species N_b |
|---|---|---|
| differences from S's rep | `R_Q + K_Q + U_Q` | `R + (1-b)L + P_N` |
| `ancestry_agreement` | `1 - R_Q/(L+R)` | `bL/(L+R)` |
| agreement at fixed sites (R made polymorphic by the alleles) | 1 | `b` |
| `ancestry_fixed_gain` (fixed minus all; even weights) | `R_Q/(L+R)` | `b·R/(L+R)` |
| `allele_explained_share` (Q's sister allele stored) | `(R_Q+K_Q)/(R_Q+K_Q+U_Q)` | `R/(R+(1-b)L+P_N)` |
| polymorphic sites with an allele's base (`polymorphic_known_share`'s numerator) | `R_Q + K_Q` | `R` |

Sequencing errors add to every denominator.

**What the table predicts. These are the test's hypotheses, not facts:**
- **The ancestry sites mix two histories.** Q disagrees at `R_Q`: a deep strain carries the congeners' base at the
  representative's own derived sites. This is the v19 observation that real strains side with the congeners at the
  congener sites (AUC 0.37-0.46). The world reproduces it by construction. That supports its realism.
- **Only the fixed sites remove `R_Q`.** There the margin is `1 - b`, whatever the strain's depth in S. Everything
  else shrinks as `b → 1` and `R_Q → R`. That pairing (a late novel twin against a deep strain) is the residual error
  class at r226.
- **`ancestry_fixed_gain` alone is not ordered by class.** `b·R` can exceed `R_Q`: a strain near the representative
  against a late twin. It separates only together with `ancestry_agreement`, whose sum approximates the fixed-site
  agreement.
  - The approximation is not exact. The gain weights sites by allele cover and leaves out uncovered sites;
    `ancestry_agreement` counts every site.
  - The test computes the true fixed-site agreement as a ceiling.
- **The allele and polymorphic features are confounded by `R`.**
  - Every genome outside the representative's subclade, a novel congener included, carries the alleles' base at `R`.
  - `allele_explained_share` favours Q only where `R_Q + K_Q > R`, that is, for strains with an allele relative.
  - A strain of a lineage without alleles (`K_Q = 0`) can score lower than a late twin.
  - `polymorphic_known_share` behaves the same way. Its numerator for N is all of `R`.
- **With 1-2 DB congeners the sites change.** With fewer than 3 compared, the sites are the nearest congener's
  differences, which include that congener's own derived states. There S's base is the ancestral one, so N carries
  it as well. N_0's agreement then rises from about 0 to about 0.5.
- **Allele scores move only unsure reads.** `SettleBySites` acts only on short reads with another species' candidate
  within 3 mismatches. Outside twin genera, at about 5% or more apart, almost no read qualifies.
- **In twin genera the build drops deep alleles.** It rejects an allele as far from the representative as the nearest
  congener's copy (`beyond_congener`), so a wide species with a DB twin keeps only its shallow alleles.
- **At 8 allele genomes, sister alleles may not be stored.** Farthest-first selection (K = 4) prefers deep alleles, so
  a near strain's sister allele may be missing from the table.
- **Fixed sites need an allele genome on the strain's side.** (Found by the generator's own test, seed 5.) A deep
  strain's `R_Q` sites are polymorphic only where some allele genome also branches off the representative's lineage
  above them. If every allele genome sits on the representative's side of the crown, the deep strain disagrees at
  "fixed" sites: in that world 0.81 instead of 0.99. The test records this coverage per strain
  (`at_rep_lineage_poly`) and bins the sensitivity by it.
- **`BestAllele` is a step, not a line.** A read is explained by an allele only when its shift is negative: it must
  carry more than half of that allele's edits in its span. A strain that shares a quarter of a stored allele's
  variants therefore gains almost nothing. So the planted P(f) line should rise steeply near f = 0.5, not as
  `explained = f` (H1 is corrected accordingly).

## 3. The world

**Generator.** A new `scripts/mini_db/simulate_ancestry_world.py`. It reuses `simulate_gtdb_release.py`'s marker table,
`mutate` (made to record each event), assembly and release writers, so the output is a release
`gtdb_to_protal_db.py` and `protal --build` take unchanged.

- **Markers:**
  - bac120 at the mini sizes;
  - every genome 200 kb with no marker loss;
  - `--gene_rates none`, so that gene-rate patterns cannot tell strains from species;
  - substitutions only (codon indels later, for the ancestry indel sites);
  - codon-aware, with no stop codons.
- **Per genus:** a random binary species tree.
  - The target S's nearest DB congener is at 0.05-0.08 marker distance, and the genus diameter is under 0.14, inside
    the 0.15 of `species_neighbours`.
  - Twin genera add a DB twin S' at 0.02-0.03.
- **Inside S:** a coalescent genealogy of 24 tips, with the widest pair 0.02-0.05 apart (drawn per genus), which makes
  wide species, as the FN anatomy found. The roles are drawn from the tips:
  - **rep:** the representative;
  - **allele genomes (`n_a`):** in the full reference only, drawn at random from the tips the strains do not reserve.
    This is what the hash split does; the build gets `--allele_genome_share 1`, since the release holds no other
    non-representative genome;
  - **Q_near ×3:** each the sister of an allele genome. With `n_a = 0` they become novel-lineage strains;
  - **Q_deep ×2:** in a clade without any allele genome, with MRCA(Q, rep) at the crown;
  - **Q_rep ×1:** inside the representative's subclade, the easy case;
  - **Q_leak ×1:** one of the allele genomes itself, the leaky positive control, kept apart in every analysis.
- **Novel species N_b on S's stem:**
  - `b` ∈ {0, 0.5, 0.8, 0.95, 1.0};
  - each own branch is drawn so that d(N_b, S rep) matches the Q_deep distances where `(1-b)L + R` allows;
  - **N_1.0 is a decoy:** it leaves at S's crown and evolves exactly as a deep strain, but is labelled a species of
    its own. No feature should separate it from Q_deep, so any separation is a leak in the world or the evaluation.
- **Factors over genera:**
  - DB congeners `n_c` ∈ {2, 6}, for the consensus against the nearest congener's differences;
  - twin ∈ {no, yes};
  - `n_a` ∈ {0, 2, 8};
  - 12 cells × 2 replicate genera = 24 genera.
- **Background:** 8 unrelated genera of 2-4 DB species each, for realistic samples.
- **Arms**, as extra query genomes of the same world:
  - **Q_rec** (built as Q_imp0.2 and Q_imp0.4, section 9): a Q_near whose 20% or 40% of genes carry a ~400 bp segment
    of the sister's copy, recombination with a congener; and N_imp, a novel species with segments of S's genomes;
  - **Q_ils:** a strain of a clade that keeps the ancestral base at 10% of S's stem sites (retained ancestral
    polymorphism). Run once with an allele genome of that clade and once without;
  - **planted genomes**, for four genera with `n_c = 6` and `n_a = 2`: S's representative with exactly 2% differences
    per gene placed by class:
    - P(f) has `f` of them at the alleles' polymorphic sites with an allele's base, the rest at non-sites;
    - A(a) gives the congeners' base at `1 - a` of the truth ancestry sites, the rest at non-sites;
    - f, a ∈ {0, 0.25, 0.5, 0.75, 1}.

  With `n_a = 2 ≤ K`, every allele genome is stored, so the planted sites are protal's.

**Outputs beside the release (`simulation/`):**
- `roles.tsv`: each genome's genus, target, role, `b`, its true distance to the target's representative, and
  `R_Q`, `K_Q`, `U_Q`, `L`, `R` per gene;
- `truth_sites.tsv.gz`: per target copy, every position where any genome of the genus differs, with the branch it
  arose on;
- `query_genomes.tsv`: the query and novel genomes, a genome table simulate_metagenomes reads;
- `tree.nwk`.

**A second world, `W_train`.** Another seed and 24 other genera, built as its own database. Only the call-level
models are trained on it.

**Honesty rules**, checked by the generator's tests:
- no query or novel genome is in the release or the full reference;
- every genome has the same length and marker set;
- query and novel focus genomes get the same depth distribution;
- the two classes are compared on the same taxon's rows.

## 4. The samples

**Paired design.** Each sample takes one focus genome per test genus, half from the Q roles and half from the N roles,
plus 10-20 background species.
- **Truth file (`--profile_truth`):** S's lineage for a Q, and N's own lineage, which the database lacks, for an N.
- **What a row then is:** S's row is a true positive or a false positive with everything on the database side equal.
- **Depth:** 3, 10, 30 or 100 marker fragments per focus genome, drawn independently of class.

| set | tool | samples | for |
|---|---|---|---|
| pe 2×150 | `simulate_metagenomes --from_manifest` (a written manifest: exact read pairs per genome; HS25, quality ~35) | 96 | everything |
| se | the same reads' R1 | 96 | features without mates |
| HiFi | `--long_samples` / `--long_genomes` | 48 | features on long reads (no allele scores there) |
| error-free pe | `simulate_reads.py --error_rate 0`, planted genomes and controls as focus genomes | 24 | the oracle, calibration lines |

- **pe profiling runs:** each pe set is profiled twice on the same database, by default and with
  `--no_allele_scores`.
- **Shuffled-allele run:** a third run uses a copy of the database whose `strain_alleles.tsv` was shuffled across
  species of the same gene (`--add_tables`). The positions stay valid, since substitutions only keep the lengths.
- **Per-read truth:** read names give each read's origin: `<contig>-<n>` for Illumina, `<accession>-<n>` for
  `simulate_reads.py`, and `g<i>x_<n>` through the long-genome table for HiFi.

## 5. What is measured

**L0, the tables against the history.**
- **`strain_alleles.tsv`** (from `protal --unpack_db`): each allele's edits against the true differences of the allele
  genome it matches, inside its range. Also checked: which genomes K = 4 chose, and the build's `Strain alleles:`
  line (`beyond_congener` in twin genera).
- **The ancestry sites:** computed by the oracle in protal's way (the consensus over the DB congeners' references,
  with `kConsensus` 0.9 and fallback to the nearest congener) and compared with the history, by class (stem, rep
  lineage, congener lineages).
- **No protal change needed.** A dump of protal's own site lists is only worth adding if L2 shows a disagreement
  that the site sets could explain.

**L1, reads.** For each focus read, the run's primary record:
- its taxon, MAPQ, kept or not, and whether S was a candidate at all (`ZA`/`ZC`);
- by default against `--no_allele_scores`.

What is reported:
- which reads changed;
- by origin, the share the scores moved to their true species;
- the per-sample `strain alleles:` counts against the reads that actually differ between the runs;
- in twin genera, the MAPQ gained by Q reads on S.

**L2, features.** S's rows from `.truth_annotated`. For each feature:
- **Separation:** the AUC of Q against N_b, overall, within bands of protal's `identity` (0.005 wide), and per
  stratum (`b`, `n_a`, `n_c`, twin, depth, Q role).
- **Increment over identity:** the gain over identity and depth in a logistic model, with CIs by resampling genera
  (rows of one genome are not independent).
- **The oracle:** a Python recomputation of every feature from the same SAM records, with the history's site sets
  and the stored alleles:
  - protal against the oracle tests the implementation;
  - the oracle with true sites against protal's sites tests the definition;
  - the oracle's true fixed-site agreement is the ceiling the features approximate.
- **Calibration lines:** on the planted genomes, `allele_explained_share` against `f` and `ancestry_agreement`
  against `a`.

**L3, calls.** GBM models (`machine_learning_cmdline.py`) trained on `W_train`'s tables and tested on `W_test`'s
(`--test-file`), with four feature sets:
- the default without `ancestry`, `alleles` and `polymorphic`;
- with `ancestry` added;
- with `alleles` added as well;
- the default set.

What is compared:
- row AUC, F1 at the knob, FP and FN, for all rows and for the hard strata;
- the default model on `W_test`'s `--no_allele_scores` and shuffled-allele runs.

The `lu_per_kb` and `uniqueness` features already see part of the ancestry signal through the index, so the increment
is over them.

**Design checks**, to read before any result:
- **D1:** identity alone has AUC ≤ 0.6 inside each matched band. Otherwise re-match.
- **D2:** N_1.0 against Q_deep is about 0.5 for every feature. Otherwise there is a leak.
- **D3:** S is a candidate equally often for both classes. Otherwise there is an alignment confound.
- **D4:** allele counts per copy are as designed.

## 6. Pre-registered predictions and pass criteria

| | prediction | pass |
|---|---|---|
| H1 mechanics, alleles | Q_leak on error-free reads: `allele_explained_share` ≈ 1, `polymorphic_novel_share` ≈ 0; P(f): explained rises with f, steeply near 0.5 (`BestAllele` is a step, section 2); protal equals the oracle | ≥ 0.95 and ≤ 0.02; monotone; oracle exact; else a bug in the table, `FromSamRecord`, `BestAllele` or `CountSites` |
| H2 mechanics, ancestry | A(a): `ancestry_agreement` = a (n_c = 6, error-free); sites per record = oracle | ± 0.03; ± 2% |
| H3 sites against history | consensus sites = stem + rep lineage at n_c = 6; with n_c = 2 the nearest congener's derived sites enter; N_0's agreement ≈ 0 against ≈ 0.5 | recall ≥ 0.95 of compared stem sites, ≤ 2% congener-lineage sites at n_c = 6 |
| H4 ancestry separation | AUC(`ancestry_agreement`) near 1 for b ≤ 0.5; falls for b ≥ 0.8 against Q_deep; follows the table | protal within 0.03 of the oracle in every stratum |
| H5 fixed sites | true fixed-site agreement ≈ 1 for Q and ≈ b for N_b where the alleles cover R; the gain alone not ordered by class; agreement + gain recovers most of it | increment of {agreement, gain} over agreement > 0 (CI) at b ≥ 0.8; gap to the oracle's ceiling reported; > 0.05 AUC → a direct fixed-agreement feature |
| H6 confound by R | `allele_explained_share` and `polymorphic_known_share`: Q_near > N_b, but N_b ≥ Q_deep possible, as the table says; at n_a = 8 about half the Q_near lose their sister allele to K = 4 | signs as predicted by stratum; the share of sister alleles stored from L0 |
| H7 settling | twin genera: unsure Q_near reads moved to S, N reads not newly drawn to S beyond the R-site shift, reads with a clear best never changed | precision of moves ≥ 0.9 by origin; no change outside unsure reads |
| H8 calls | the groups raise AUC and F1 in the hard strata (b ≥ 0.8 against Q_deep; twin genera); nothing for n_a = 0 (its values are 0 in both classes) | CI above 0 in the hard strata; none at n_a = 0 |
| H9 controls | shuffled alleles: the alleles and polymorphic increments vanish; N_1.0 against Q_deep: no separation | ≤ 0.01 AUC increment; AUC 0.5 ± 0.05 |
| H10 arms | Q_imp loses agreement in proportion to its genes with a congener's segment, N_imp gains it; Q_ils disagrees at "fixed" sites unless an allele genome of its clade is stored | sizes reported, no pass mark |

The decisions follow from which rows fail:
- **H1-H3 fail:** an implementation bug. Fix it before looking at r226.
- **H4-H5 hold but H8 shows nothing:** the models do not use a signal that is present. Check the feature forms (the
  gain against a direct fixed agreement), the weights, and the training data's share of hard rows.
- **All hold:** the remaining question is how much of this structure real data have and how much the stored alleles
  capture. Only the r226 builds answer that.

## 7. What it needs

**Scripts (Python, as the repo's other scripts):**
- `scripts/mini_db/simulate_ancestry_world.py`: trees, roles, the release, query, novel and planted genomes,
  manifests, long-read tables, truth files, `roles.tsv` and `truth_sites.tsv.gz`. Estimated 600 lines.
- `scripts/ancestry_truth_test.py`: L0-L3, the oracle (reusing `ancestry_sites.py`'s readers and
  `insilico_strains.substitutions`), the strata and bootstrap, and a summary in the report style. Estimated 700 lines.
- `scripts/mini_db/test_ancestry_world.py`:
  - the generator's tests: tree distances, roles, no query genome in the release, planted counts;
  - a smoke run of 2 genera, 6 error-free samples and H1/H2 directions, cheap enough for CI (~2 minutes).
- `scripts/mini_db/run_ancestry_truth_test.sh`: two worlds, two builds, the simulations, five profiling runs, training,
  evaluation.

**protal:** no change. A site dump only if L2 calls for it.

**Compute (estimates):**
- the worlds, seconds;
- two mini builds of about 300 species, minutes each at `-t 4`;
- about 2 M pe pairs per world, a few minutes to simulate and to profile per run;
- GBM fits on a few thousand rows, minutes.

Under an hour on 4 cores (`taskset -c 0-3`), or one cluster job.

## 8. What it cannot tell

- How much of this tree structure real species have, and how much of it 4 alleles from half a species' genomes
  capture. That is the r226 builds' question (the next one is the first with `504de80` and `ea6de54`).
- Recombination and indel spectra beyond the two arms.
- Allele scores for long reads, which do not exist yet.
- A genus as crowded as large GTDB genera: `--align_top` and the adaptive candidates are only checked by D3, not
  stressed.

## 9. Implementation (2026-10-09)

The user asked to implement the design, and to make it also test sensitivity: recently against long-diverged
lineages, and the potential for recombination with congeners.

### Files

| file | what |
|---|---|
| [simulate_ancestry_world.py](../../scripts/mini_db/simulate_ancestry_world.py) | the worlds: genus trees, coalescent genealogies, roles, recombination, the derived and planted genomes, the release, the truth tables, the samples (pe manifest, HiFi tables, error-free communities, truth files) |
| [ancestry_oracle.py](../../scripts/ancestry_oracle.py) | a line-by-line Python port of `AncestrySites.h`, `StrainAlleles.h` and the profiler's sums: recomputes the 13 features from a run's SAMs and the database's files |
| [ancestry_truth_test.py](../../scripts/ancestry_truth_test.py) | `run` (both worlds, builds, reads, five protal runs, tables, models; resumable) and `evaluate` (L0-L3, the verdicts, `OUT/report/summary.md` and a TSV per table) |
| [test_ancestry_world.py](../../scripts/mini_db/test_ancestry_world.py) | 16 tests: the generator's invariants, the oracle on hand-made cases, and an end-to-end `--quick` run with `$PROTAL` and `$SIMULATE` |
| [run_ancestry_truth_test.sh](../../scripts/mini_db/run_ancestry_truth_test.sh) | the SLURM job in the user's template: the work on the node's SSD, a tarball of the report back to `/hpc-home` |

No protal change.

### The two additions

**Recently against long-diverged lineages.**
- Each genus draws its species' width: the crown height `c_S` is 0.005-0.02, or 0.003-0.011 with a twin. It also
  draws its stem length `L`.
- Each strain records how long ago it shared an ancestor with the representative (`mrca_rep`) and with the nearest
  allele genome (`mrca_nearest_allele`, and which genome).
- The roles span the range:
  - Q_rep and Q_near diverged recently from the representative or from a known strain;
  - Q_lone diverged recently from a strain the database does not know;
  - Q_deep belongs to an old lineage without allele genomes.
- The novel species left the stem early (b = 0, 0.5) or late (0.8, 0.95, and the decoy at 1).
- The sensitivity table bins the Q rows by distance to the representative, by the MRCA with the nearest allele
  genome, and by the share of the representative's own derived sites they lack that the alleles cover.

**Recombination with congeners.**
- Three regimes, crossed with the other factors over the genera:
  - none;
  - low: per gene, an import from another S genome with probability 0.1, and from another species of the genus with
    0.03;
  - high: 0.3 and 0.15.
- Segments are 30 + exponential(400) bases, at a uniform place. A boundary codon that would be a stop keeps the
  recipient's bases.
- Donors can be any other species of the genus, the novel ones included. Every genome of the genus can receive:
  S's genomes (representative and allele genomes too), the congeners' representatives and the novel species.
- Directed arms in every genus:
  - Q_imp0.2 / Q_imp0.4: a strain with the sister's segments in 20% or 40% of its genes;
  - N_imp0.2 / N_imp0.4: N_0.5 with S's segments.
- `genes_congeneric` counts only segments from other species. The sensitivity table bins by it (the Q rows) and by
  the novel species' imports (the N rows).

### Changes against the design

- **A role added:** Q_lone, a strain whose sister tip is no allele genome. Q_rec became Q_imp and N_imp.
- **The truth uses protal's sites.** The truth's consensus sites, and the planted A_a genomes, are protal's (the
  oracle's port of its rule on the database's copies), not sites by position. The first smoke run planted by
  position, and missed `a` by up to 0.039. This rested on two differences, which L0 reports separately:
  - protal compares along shared 12-mers;
  - it needs 0.9 of the congeners compared at a position to agree.
- **H3's thresholds are kept as pre-registered.** Its failure is a property of the rule (below), not of the code.
- **The decoy check (H9) uses only genera without recombination.** There, N_1 takes segments from S as another
  species does, a real difference from a strain of S.
- **Models:** paired-end reads get all five feature sets, the other read types base and default
  (`--model-sets`, `--model-types`).

### Smoke run (`--quick`: 2 test genera per world, 8 pe samples; WSL, 4 cores, 74 s)

The protal build `~/bpar/build/protal`, whose `AncestrySites.h`, `StrainAlleles.h`, `StrainAllelesBuild.h`,
`Profiler.h`, `AlignmentStrategy.h`, `Build.h`, `RunProtal.h` and `Options.h` are byte-identical to HEAD (`6ae06fc`).

- **The implementation is exact.** protal's 13 features equal the oracle's on every target row, in every run:
  error-free, paired-end with errors, the shuffled table, and `--no_allele_scores`.
- **The table:** 99.5-100% of the stored alleles are exactly an allele genome's differences in their range.
- **The leaky control:** an allele genome's own reads are explained at 0.987, with no novel polymorphic sites.
- **The planted lines:**
  - P_f is a step, as predicted: f 0 / 0.25 / 0.5 / 0.75 give `allele_explained_share` 0 / 0.10 / 0.44 / 0.82.
  - A_a: `ancestry_agreement` within 0.013 of `a`.
- **The allele scores in the alignment:**
  - all 26 strain reads they moved went to the target;
  - they also drew 73 of 752 novel species' reads onto it, the confound section 2 predicted.
- **H3 failed: protal's sites are not S's history.** Over the genus without recombination, 0.956 of the sites lie on
  S's stem or the representative's lineage, and they recall 0.84 of those.
  - **Missed:** 423 sites because one of the six congeners carries another base (0.9 of 6 needs all six), 78 not
    compared, and 5 at positions fewer than 3 congeners reached.
  - **Taken that are not history:** 98 by consensus (homoplasy: S differs from every congener where its derived
    state arose above its split, e.g. where the sister reverted), and 30 that are the nearest congener's own states.
- **Model verdicts:** meaningless at this size (8 training samples).

### Full run

**Data.**
- Run locally at the user's request ("if under 30 min"): WSL, `taskset -c 0-3`, the `~/bpar` build, 17 min end to end.
  The second run, after the decoy fix, reused the train world, which `--no_decoy` leaves byte-identical.
- Command:

      ancestry_truth_test.py run --outdir ~/atp_full -t 4 --model-types pe,hifi --world-args "--long_samples 32"

- Each world: 36 test genera (2 or 6 congeners × twin or not × 0/2/8 allele genomes × recombination none/low/high)
  and 8 background genera; 96 pe samples (also profiled as se), 32 HiFi, 24 error-free.
- The report and every table are in [full_run/](full_run/summary.md) (the verdicts in
  [verdicts.tsv](full_run/verdicts.tsv)).

**Verdicts.**

| | outcome | verdict |
|---|---|---|
| H1b/H2b protal = oracle | every target row, all 13 features, in the default, `--no_allele_scores` and shuffled runs, error-free and with errors | PASS |
| H1a the table | stored alleles exactly an allele genome's differences: 100% | PASS |
| H1 leak, calibration | allele genome's own reads explained 0.965; P_f 0 / 0.052 / 0.285 / 0.551 / 0.691 (a step, as predicted); A_a within 0.009 of `a` | PASS |
| H3 sites = history | 6 congeners, no recombination: 0.965 on S's stem or the representative's lineage, recall 0.798 | FAIL |
| H3 fallback | 2 congeners: half the sites are the nearest congener's own states (0.502) | PASS |
| H4 ancestry, b ≤ 0.5 | AUC 0.904 / 0.807 overall; 1.000 / 0.901 with 6 congeners, alleles, no recombination | FAIL overall |
| H5 fixed sites, as registered | deep strains with their representative-lineage sites covered against b 0.8/0.95: `ancestry_agreement` 0.647 → with `ancestry_fixed_gain` 0.879 (+0.232 [0.149, 0.312]); the true fixed-site agreement alone 0.834 | PASS |
| H6 the confound | `allele_explained_share`, late twins against Q_near 0.904, against Q_lone 0.320 | PASS |
| H7 settling | 3,187 strain reads moved, all to their species; also 6,072 of 115,881 novel species' reads drawn onto the target (5.2%; 7.5% in twin genera) | PASS (and the confound) |
| H8 calls | default against base, AUC on the target rows: pe +0.073 [0.040, 0.112], hard strata +0.150; HiFi +0.067, hard +0.162 | PASS |
| H9 decoy | no feature separates the decoy from Q_deep (largest deviation 0.037) | PASS |
| H9 shuffled | the default model on a shuffled table: AUC 0.918 → 0.829, below base's 0.845 | PASS |

**What works.**
- **The code computes what it claims.** The oracle equality covers every feature and run, so a weak feature at r226
  is not an implementation fault.
- **The alleles carry most of the gain:**
  - pe models' AUC on the target rows: base 0.845, +ancestry 0.844, +alleles and polymorphic 0.922, default 0.918;
  - F1 over all test rows 0.839 → 0.903;
  - false calls of late twins (b 0.95) at the knob 0.335 → 0.159.
- **The polymorphic idea works where its premise holds.** H5 passes: `ancestry_fixed_gain` lifts deep strains against
  late twins by 0.23 AUC once the alleles cover the representative's lineage.
- **No leak.** The decoy check is clean.

**What does not.**
- **The ancestry group adds nothing to the models.**
  - Alone over identity and depth it adds +0.040 AUC, but the base set already carries it: `lu_per_kb` alone has AUC
    0.798, `ancestry_agreement` 0.810.
  - This matches the r226 v21 ablation (+0.001-0.002).
- **The consensus rule loses a fifth of S's derived sites with diverse congeners.** It needs 0.9 of those compared,
  i.e. all six. Of the missed history sites, ~85% have one congener with another base (L0 table).
  - With 2 congeners, half the sites are the nearest congener's own states. There a novel species that left early
    (b = 0) agrees at them: AUC 0.971 there against 0.993 with 6 congeners.
- **The alleles cost the strains they cannot explain.** The models learnt "explained means present", so strains of
  lineages the alleles do not reach lose sensitivity:

  | pe, share called at the knob | base | default |
  |---|---|---|
  | MRCA with the nearest allele genome < 0.002 (recently diverged from a known strain) | 0.716 | 0.984 |
  | 0.002-0.01 | 0.622 | 0.781 |
  | ≥ 0.01 (long diverged) | 0.312 | 0.281 |
  | species without allele genomes | 0.661 | 0.593 |
  | representative-lineage sites the alleles cover: all | 0.673 | 0.969 |
  | none | 0.669 | 0.601 |
  | Q_deep / Q_lone / Q_near / Q_rep | 0.452 / 0.740 / 0.722 / 0.995 | 0.594 / 0.740 / 0.991 / 0.985 |
  | distance to the representative ≥ 0.03 | 0.192 | 0.448 |
  | N_0.8 / N_0.95 / N_1 (false calls) | 0.158 / 0.335 / 0.395 | 0.074 / 0.159 / 0.285 |

  HiFi shows the same shape: recent 0.749 → 0.996, long-diverged 0.550 → 0.300, without allele genomes
  0.667 → 0.714.
- **K = 4 farthest first drops the near alleles.** With 8 allele genomes, a strain's nearest allele genome is stored
  in 0.54-0.56 of its genes, against 0.89-0.94 with 2 (`l0_nearest_allele_stored.tsv`).

**Recombination with congeners.**
- **Imported segments cost a strain sensitivity.** Q_imp, Q_near's genome with the sister's segments in 20-40% of its
  genes, is called at 0.873 by the default model against Q_near's 0.991 (base 0.694 against 0.722).
- **High recombination flattens the gain.** In genera with high recombination, base and default find the strains
  alike (0.800 / 0.807). The default model still halves the false calls there (0.177 → 0.085).
- **S's segments in a novel species barely raise false calls.** N_imp, N_0.5 with them, is called at 0.04-0.07 by
  either model, as N_0.5 is.
- **The allele scores of the alignment do little for the calls.** The default model loses only 0.005 AUC without
  them (`--no_allele_scores`). Yet they move 7.5% of a novel species' reads onto the target in twin genera, against
  4.0% for strains (`l1_reads.tsv`).

**What follows for protal** (not implemented):
1. **A fixed-site agreement feature.** The true fixed-site agreement adds +0.070 over identity and depth in genera
   with alleles. The models reach it only through `ancestry_agreement` + `ancestry_fixed_gain`, which H5 shows works;
   a direct feature would save the models that step.
2. **Alleles chosen for spread across the genealogy, not farthest first.** Then a strain's side of the tree, and its
   nearest relatives, are covered (the coverage bins above).
3. **A consensus that tolerates one dissenting congener at 6 or more.** For example, `kConsensus` as "all but one",
   to recover the ~20% of derived sites lost.
4. **Training worlds whose strains include long-diverged lineages without allele coverage.** Otherwise the models
   learn that unexplained means absent.

**Website.** No user-facing behaviour changed; the website needs no update.

## 10. The four changes, implemented (2026-10-09, evening)

The user asked for 1-4 above. All four are in, on top of `6ae06fc` (uncommitted with the test's files), and the
full test was run again on the same worlds and seeds, so that every number below compares with the full run of
section 9 on byte-identical worlds and reads.

| change | where | what |
|---|---|---|
| 1. `ancestry_fixed_agreement` | `Profiler.h` (`AncestryFixedAgreement`), `model_features.py` (the `polymorphic` group, 86 default features), the oracle, `features.md` | the share of the species' base at the ancestry sites fixed within the species, weighted as the gain is (n/(n+1) of the n alleles covering a site); -1 without the table, 0 without a site |
| 2. alleles chosen for coverage | `StrainAllelesBuild.h` (`Coverage`, `Select`), `--strain_alleles` help, `databases.md` | of the 16 sampled alleles, up to 4 chosen greedily to maximise the sum over the sample of the best coverage a chosen allele gives each: `Coverage(a, c)` = the edits of `a` that `c` shares, less the edits of `c` that `a` lacks, over `a`'s edits (what `BestAllele` would save a read of `a`'s genome, 0 when `c` shares no more than half of its own edits), in integer units so that ties are exact; every sampled allele weighs the same, whatever its distance from the representative. Farthest first before |
| 3. one dissenting congener | `AncestrySites.h` (`Agree`, `kTolerantFrom` 6), the two Python ports (`ancestry_oracle.py`, `ancestry_sites.py`), `features.md` | from six congeners compared at a position, all but one may carry the consensus base when the one carries a third base (neither the species' nor the others': its own change there); a congener with the species' base still blocks the site, as it may share the state by descent, and then a novel species below the site carries it too. Indels alike (another indel at the same position). The nine-in-ten rule otherwise |
| 4. training strains the alleles do not reach | `insilico_strains.py --multi-share` (`multi_genome_species`), `build_gtdb_database.py --insilico-multi` (default 0.5), `databases.md` | half of the species with real strains in the genome table get an in-silico strain of their representative as well: a mutated copy shares none of the alleles' variants, as a strain of a lineage GTDB has not sampled does not. The r226 build splits a species' other genomes by hash into alleles and simulated strains, so every simulated real strain had a close allele |

Tests: `AncestrySites.ConsensusToleratesOneCongenerWithABaseOfItsOwn`, the `Coverage` and the new `Select` order in
`StrainAlleles.TheBuildsSampleAndChoiceDoNotDependOnTheOrder`, the feature in the polymorphic features test, the
oracle's rule in `test_ancestry_world.py`, `--multi-share` in `test_insilico_strains.py`, the group lists in
`test_model_pmml.py`. Unit 479 (3 skipped as always), e2e 149 (on a mini database built by this binary; the
checkout's `data/mini_db` is a stale copy of the old layout), pipeline 6 (with `--insilico-multi 0.5` in effect),
scripts 41 + 13 + 8, ancestry world 17 (with the end-to-end quick run): all OK on WSL, 4 cores.

The design's item 4 for the test's own worlds needs nothing: the train world already has Q_deep and Q_lone, strains
without an allele relative; the change is for the GTDB build, so this run tests changes 1-3.

### The full run again ([full_run2/](full_run2/summary.md), 16 min, the same command and seeds as section 9)

The worlds and reads are byte-identical to section 9's (the generator and the simulators are seeded; the L0 rows of
the two-congener genera match to the site). What changed is the database's alleles (change 2), protal's sites
(change 3) and one feature (change 1). Every verdict of section 9 keeps its outcome; H1b/H2b (protal = the oracle)
holds with the new rule and the new feature on every row.

**Change 3, the consensus.** At six congeners without recombination the species' derived sites are recalled at
0.90-0.92 instead of 0.83-0.86 (`l0_sites.tsv`, genera 19/22/25/28/31/34). The price is homoplasy: sites where a
congener's own change left the five others and the species apart went from 59-74 to 161-199 per genus, so the share
of sites that are the species' history fell from 0.97-0.98 to 0.94-0.96. The sites still missed are the ones where the
dissenter carries the species' base (not tolerated, by design) and the positions no congener was compared at.
`ancestry_agreement` separates Q_deep from N_0 at 0.906 (0.904) and from N_0.5 at 0.816 (0.807): the rule was not
what limited H4, which stays FAIL overall (the two-congener fallback and recombination are).

**Change 2, the alleles.** With eight allele genomes a strain's nearest allele genome is now stored in 0.59 of its
genes against 0.55 (Q_near 0.52 → 0.56, Q_rep 0.37 → 0.53); with two, 0.92 as before. The stored edits fell from
116,863 to 106,298: nearer alleles. Four of eight cannot cover every lineage, so the gain is bounded; what moved is
which four.

**Change 1, the feature.** `ancestry_fixed_agreement` equals the oracle's ceiling in every comparison and stratum
(`l2_auc_pe.tsv`: 0.652 against 0.650 over all rows, 0.839 against 0.839 in genera with alleles, H5's stratum 0.837
against 0.837). The `polymorphic` group alone over identity and depth, in genera with alleles, now adds 0.107 [0.071,
0.144] of AUC against 0.031 [−0.008, 0.075] before; all three groups 0.167 against 0.162 (pe), 0.176 against 0.169
(se), 0.169 against 0.166 (HiFi).

**The models** (`l3_models.tsv`, the test world's rows, section 9's run in brackets):

| model | target rows AUC | hard: Q_deep vs N_0.8/0.95 | all rows F1 at the knob | sensitivity / FPR at the knob |
|---|---|---|---|---|
| pe base | 0.846 (0.845) | 0.647 (0.641) | 0.835 (0.839) | 0.785 / 0.068 (0.785 / 0.062) |
| pe default | 0.923 (0.918) | 0.803 (0.791) | 0.906 (0.903) | 0.877 / 0.043 (0.868 / 0.039) |
| pe alleles (no ancestry) | 0.926 (0.922) | 0.810 (0.801) | 0.906 (0.901) | 0.879 / 0.044 |
| HiFi base | 0.870 (0.872) | 0.682 (0.683) | 0.850 (0.851) | |
| HiFi default | 0.940 (0.939) | 0.854 (0.845) | 0.923 (0.921) | 0.901 / 0.049 (0.903 / 0.054) |
| pe default on shuffled alleles | 0.830 (0.829) | | 0.816 (0.822) | |

- The default model gains a little everywhere: +0.005 of AUC on the target rows, +0.012 in the hard stratum, +0.003 of
  F1 for pe; +0.009 in HiFi's hard stratum. The confidence intervals of a single run are ±0.03-0.07, so these are
  consistent, not established.
- The sensitivity by bin (`l3_sensitivity.tsv`, pe, called at the knob, default model; section 9's in brackets):
  strains whose MRCA with the nearest allele genome is 0.002-0.01 away 0.824 (0.781), no allele genome 0.621
  (0.593), none of the representative's lineage covered 0.625 (0.601), Q_deep 0.632 (0.594), Q_lone 0.762 (0.740),
  Q_rep 0.990 (0.975); recently diverged strains stay at 0.98. The false calls rose with them: N_0.95 0.198 (0.159),
  N_0.8 0.088 (0.074), the knob's FPR 0.043 (0.039). The model trades along the same curve; its AUC moved up.
- The ancestry group still adds nothing the base set lacks: `pe alleles` without it is as good as the default (0.926
  against 0.923), as in section 9 and the r226 v21 ablation. The direct fixed-site feature is in the `polymorphic`
  group, so its gain shows there, not in `ancestry`.

**What this says for r226.** Changes 1-3 are small, consistent gains here, where the world has no among-site rate
variation beyond codon position and four alleles out of eight. Change 4 acts on the training mix, which this test
cannot see; its effect, and whether the near alleles now stored help real strains, is the next GTDB build's question
(`--insilico-multi 0.5` is the default; `0` restores the old mix).

**Website.** No user-facing behaviour changed; the website needs no update. A database built from now on stores
other alleles than before (change 2), and a run on an old database computes the sites by the new rule (change 3)
and one more feature (change 1); models trained before this change lack the feature and keep working, as the model
reads features by name.

## 11. The deep alleles' k-mers in the index (2026-10-10)

The user asked for the alternative costed in the [site-weighted evidence plan](../2026-10-09-site-weighted-evidence-plan.md),
section 8, restricted to the alleles at least 1% from the representative. Implemented as `--index_alleles`
(default 0.01; 1: none):

- **Build** (`StrainAllelesBuild.h AlleleSeeds`, `Build.h`): in both index passes, after a copy's own k-mers, the
  k-mers of the copy with each deep allele's substitutions applied that the copy does not yield (windows over an
  indel, whose inserted bases the table does not keep, and over a base that is not A, C, G or T are left out, as
  are k-mers an earlier allele gave) go in under the species at the representative's positions, with the unique
  flag off (`SetFlagNonUnique`: a strain's k-mer is not the representative's). `unique_kmers.tsv` leaves them out
  of each copy's totals (`CountUniqueKmers` `left_out`), so a species' unique shares stay the representative's. The
  log: "Index alleles: N seeds of A strain alleles at least 0.01 from the representative's copy, of C gene copies".
- **Run**: nothing changes. A read of a deep strain now seeds on its species; its anchored alignment finds the
  link not exact and falls back to the whole-window alignment (`AnchoredAligner::NotApplicable`), scored with the
  alleles as before.
- **Tests**: `StrainAlleles.DeepAllelesGiveTheIndexSeedsTheRepresentativeLacks`; unit 479, e2e 149 (a mini
  database with 17,029 seeds of 288 deep alleles), scripts 41 OK.

**The full run again** ([full_run3/](full_run3/summary.md); the same worlds, reads and seeds as section 10's run,
[full_run2/](full_run2/summary.md), which is the comparison):

| | without the seeds (full_run2) | with the seeds (full_run3) |
|---|---|---|
| test world's index | 5,025,954 syncmers, 3.14 GB after pass 2 | +264,457 seeds of 3,605 alleles ≥ 0.01 of 1,918 copies (+5.3% of entries), 3.20 GB |
| strain reads with their best record on the target, no twin / twin (`--no_allele_scores`) | 0.938 / 0.887 | 0.944 / 0.888 |
| novel species' reads on the target, no twin / twin | 0.833 / 0.757 | 0.846 / 0.757 |
| pe base model: target-row AUC, hard stratum, all-row F1, sensitivity | 0.846, 0.647, 0.835, 0.785 | 0.858, 0.671, 0.853, 0.816 |
| pe default model: the same | 0.923, 0.803, 0.906, 0.877 | 0.920, 0.797, 0.909, 0.881 |
| HiFi base / default, all-row F1 | 0.850 / 0.923 | 0.857 / 0.921 |
| pe default, strains ≥ 0.03 from the representative called | 0.433 (base 0.212) | 0.429 (base 0.438) |
| pe default, Q_deep / long-diverged MRCA called | 0.632 / 0.250 (base 0.458 / 0.312) | 0.642 / 0.281 (base 0.535 / 0.391) |
| pe default, N_1 / N_0.95 false calls | 0.298 / 0.198 | 0.250 / 0.189 |

**What it says.**
- **The seeds find the deep strains' reads, and the base model gains most:** +0.018 of F1 over all rows, +0.03 of
  sensitivity, strains 0.03 or more from the representative called 0.21 → 0.44 and strains of old lineages 0.31 →
  0.39. Without the allele features, the reads that now seed are the only new evidence, and it is worth a lot.
- **With the full feature set the gain is small:** +0.003 of F1 (pe), −0.002 (HiFi, noise), the hard stratum
  unchanged. The allele features had already recovered these strains from the reads that did seed; the seeds add
  reads, not a new kind of evidence. The false calls of late twins fell a little (N_1 0.298 → 0.250).
- **The confound is there too:** the novel species' reads on the target rose from 0.833 to 0.846 (no twin), since a
  late twin carries the alleles' ancestral base at the representative's own sites and so matches their k-mers.
- **Memory:** +5% of entries in this world, where 70% of the copies with alleles have a deep one. The plan's r226
  estimate of +1-2 GB stands; the real strains' divergence decides it.
- **A side effect to know:** an allele seed whose core is the representative's k-mer's core (a substitution in the
  flex part) turns a single-entry core into a flex block, so the representative's entry moves from "short unique"
  to "long unique" (here 1,262k / 1,254k → 1,223k / 1,293k, the sum unchanged). The `ref` features see the split
  change for species with deep alleles; the models are trained on the same database, so they learn it.

**For r226.** The build's `--index_alleles 0.01` is on by default. What to look at: the "Index alleles" line
(seeds and memory), the real strains' miss rate by their distance to the representative against v21, and the
same-genus false positives beside held-out congeners (the confound). The site weights (the plan's phase 0)
remain the next step for the cross-genus false positives.

**Shared checkout.** Another session was implementing the column weights in the same checkout during this work
(`ColumnWeights.h`, hooks in `Build.h`, `Options.h`, `Profiler.h`, `AncestrySites.h`, `model_features.py`,
`ancestry_truth_test.py`). This run was built from a staged snapshot of HEAD plus this session's hunks only, in its
own WSL tree (`~/bidx`), as the shared tree held the other session's unfinished code. The checkout's files carry both
sessions' edits; the staged copy is in this session's scratchpad.
