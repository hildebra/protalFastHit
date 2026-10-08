# The r226 v17 build: the new features, the errors that remain, and the clouds around the held-out species

**Data.** The r226 v17 build (SLURM 24060119, q512n10, 2026-10-07 23:48 to 2026-10-08 03:32, 84 threads, protal 0.7.9
`ffbedb3` plus the congener gaps `45bc47b`): the share archive `local/protal0.7.9_r226_v17_share.tar.gz` (465 MB),
unpacked to `local/v17` (the logs, `model_logs/` with the error reads' tables and capped SAMs, the training and test
tables). `internal_taxonomy.dmp` and `heldout_species.txt` are v15's, copied in: the archive holds the held-out list
under `model_logs/`, and it is identical to v15's (same release, genome table and seed), so the taxonomy is too.

**Questions (the user's, 2026-10-08).** The build's F1 rose markedly; analyse the false positives and false negatives
and the ancestral-state patterns. And: does the training database hold out the species that form the "cloud" of near
relatives around each marker gene copy in a way that serves the training, or should protal record the clouds in a
file that steers the simulations?

Scripts here, run in WSL (one or two pinned cores each, niced; `~/soil13/venv`; the v15 report's error-record
pipeline reused on v17's SAMs with its parser made to read the per-sample FP and FN files):

    python3 holdout_clouds.py --build local/v17 > holdout_clouds.txt
    python3 ancestry_patterns.py --build local/v17 --v15 local/v15 > ancestry_patterns.txt
    python3 ablate_v17.py --build local/v17 --out ~/v17/ablate --read-types pe -t 2        # ablate_results_pe.tsv
    python3 ../2026-10-07-error-read-signatures/parse_error_records.py --dir local/v17/model_logs/error_reads --out ~/v17/err_analysis
    python3 ../2026-10-07-error-read-signatures/fragments.py --records ~/v17/err_analysis --build local/v17
    python3 ../2026-10-07-error-read-signatures/signatures.py --records ~/v17/err_analysis --build local/v17 > signatures_v17.txt
    python3 ../2026-10-07-error-read-signatures/sensitivity.py --records ~/v17/err_analysis --build local/v17 > sensitivity_v17.txt

## Summary

- **The build.** 3 h 44 min (v15 4 h 42 min): the training database 14.5 min, its congener gaps and foreign rates 3.8
  min, collection 2 h 39 min (21 of 116 simulations streamed into protal), training 3 min, reports 14 min. The
  in-silico strains took `--omega auto`: 0.828 of the real strains' substitutions lie on third codon positions, which
  omega 0.075 reproduces (0.15 before).
- **F1 at the knob** (pe 0.73 by the bootstrap rule, the others 0.5): species held out pe 0.976 / se 0.980 / pb 0.984 /
  ont 0.982; test set 0.983 / 0.984 / 0.989 / 0.989; soil hold-out pe 0.9745, shallow soil 0.9617, moderate 0.987,
  gut 0.991. v15: species held out 0.952, test 0.9645 (0.9698 at 0.5), soil hold-out 0.9655, shallow 0.9396. The
  rows are new (other reads, scenarios and strains), so the comparison is between builds, not models:
  [section 1](#1-what-moved-between-v15-and-v17) separates what the reads changed from what the features add.
- **The new features.** By the trainer's importance, the foreign-copy rates are 2nd and 3rd in every model (0.10 and
  0.06 against identity's 0.71 in pe); the ancestry sites 12th-13th (pe, se), 8th-9th (pb), 6th and 8th (ont, 0.015);
  the gaps and the untried candidates below 0.001. The ablation ([section 2](#2-what-the-feature-groups-add)) puts
  the pe model without the foreign group at design-test F1 0.962 (0.984 with it), soil 0.955 (0.972), shallow soil
  0.927 (0.963): the foreign group carries the whole jump over v15, and it carries it because it leaks (next bullet).
  The ancestry sites add +0.0005, the gaps +0.0004, all four new groups over the v15 set +0.002.
- **A leak in the foreign-copy features.** `foreign_scanned_share` is the simulation pool in disguise: the foreign
  scan read the genomes at hand, the 49,595 downloaded ones, so a database species has scanned copies when its genome
  was downloaded; every present taxon is such a species, most absent ones (118,000 of the 143,614 database species
  never simulated) are not. Among the *absent* rows, species that are present in some other sample have a median
  scanned share of 1.0, species never present 0.0; "scanned" equals "in the pool" in 87-93% of all rows
  ([foreign_leak.txt](foreign_leak.txt)). Inside the identity band its AUC for present against novel congener is
  0.89-0.91 at every identity, which no read signal has, and the trainer ranks it 2nd. On real data the shipped
  database's scan covers the same genomes, so the 83% of GTDB species without a downloaded genome would read as
  "unscanned", a push towards absent. The ablation ([section 2](#2-what-the-feature-groups-add)) puts F1 numbers on
  it; a scan of every genome's marker copies from the full reference would remove it.
- **The errors** ([section 3](#3-the-errors-and-the-ancestry-patterns)). The same two classes as in v15: the false
  positives are novel congeners (87% of pe's FP, identity median 0.977, 3 fragments), the false negatives strains
  (67% of pe's FN, with 19% in-silico and 14% representatives; identity 0.962, 2 fragments); 40% of both are
  singletons, and both sit in the soils. The ancestry sites separate the classes as designed (agreement 0.98-0.996
  for present taxa, 0.68 for absent congeners; the congener's base at 0.24 of the sites against 0.01) and add a little
  to the base features inside the band (pe AUC 0.959 → 0.965), but at the residual errors they point the wrong way:
  the missed strains carry the congener's base at 19% of the sites, the false positives that slip through at 13%.
- **The reads changed.** The Illumina reads of v17 come from `simulate_metagenomes`' own model, not ART: the
  representatives' reads align at identity 0.996 (v15 0.991), strains' at 0.993 (0.988), and the share of present
  strains inside the 0.95-0.985 band halved (37% → 18%); the novel congeners did not move (0.936 → 0.937). Part of
  pe's gain is cleaner reads, not a better model; PacBio and Nanopore reads are as in v15.
- **The clouds and the hold-out** ([section 4](#4-the-clouds-around-the-held-out-species)). The hold-out draws 30% of
  the simulated species at random, plus whole clades; it never looks at how near a species' congeners are. Of the
  7,324 species held out alone, 91.5% keep a congener in the training database (73.5% five or more): their reads have
  a congener to land on. Of the kept species 53% have a held-out congener (91% in genera of 21 or more species). By
  distance, the held-out congeners mostly sit far from their kept neighbours: of the absent rows beside a held-out
  congener with 3 or more fragments, 2.1% are at identity 0.985 or above and 6.5% at 0.97-0.985; present strains are
  at 0.985 or above in 82% of the rows, at 0.97-0.985 in 15%. The torn clouds (a near-identical sister held out
  alone, its twin kept) are a few percent of the novel-congener rows, and they are where the remaining FP sit.
  protal records the clouds (`species_neighbours.tsv`, `congener_gaps.tsv`) only once a database is built; to steer
  the simulations they would be computed from the converted references before the hold-out (`species_clouds.tsv`),
  used to keep species complexes together, written into the collector's tables as the novel congener's distance
  (`meta_novel_distance`) so the trainer can report and weight by it, and used to draw the samples' congener groups
  by distance band. Not implemented.
- **The ancestry report of the build failed** for every read type on one line: a congener copy carrying an `N` at a
  position where it differs from the species' copy gave base code 4, out of a 4-wide vote table (`ancestry_sites.py`,
  `Sites.__init__`). It had read 30,208 counted records of 2,157 taxa and 32,222 alleles of 13,009 taxon-genes by
  then. Fixed in `401c4f5` (`insilico_strains.substitutions` and the sites skip a position with another letter on
  either copy; a test covers it). The training database's full reference was removed after the failed step, so the
  alleles' half cannot be re-run for v17; the next build makes the report.

## 1. What moved between v15 and v17

v17 is `ffbedb3` + `45bc47b`; v15 was `ce85bd7`. Between them: the 15 features against false positives in complex
communities (`6a4b16c`), the ancestry sites, the congener gaps, the foreign rates and the untried candidates (81
features against 55); scenarios of varied richness and depth with a moderate scenario (`fb6c42d`); the knob by
bootstrap (`c922d0f`: pe 0.73) and the training design's strains as the test set's; in-silico strains with
`--omega auto` (0.075); and new reads for every type (Illumina without ART, PacBio and Nanopore in
`simulate_metagenomes`). The rows do not pair with v15's, so the builds are compared on their own held-out estimates.

| | v15 | v17 |
|---|---|---|
| pe species held out / design test / soil hold-out / shallow soil hold-out | 0.952 / 0.9645 / 0.9655 / 0.9396 | 0.976 / 0.9833 / 0.9745 / 0.9617 |
| se species held out / design test / soil hold-out | 0.963 / 0.9622 / 0.9635 | 0.980 / 0.9836 / 0.9758 |
| pb species held out / design test / soil hold-out | 0.960 / 0.9785 / 0.9588 | 0.984 / 0.9894 / 0.9835 |
| ont species held out / design test / soil hold-out | 0.954 / 0.9723 / 0.9483 | 0.982 / 0.9886 / 0.9795 |

The reads' identity by class of row, 3 or more fragments ([ancestry_patterns.txt](ancestry_patterns.txt), last
section of each read type):

| pe | v15 median (10%, 90%) | in the 0.95-0.985 band | v17 median (10%, 90%) | in the band |
|---|---|---|---|---|
| present, the representative | 0.9913 (0.983, 0.994) | 0.11 | 0.9962 (0.994, 0.998) | 0.02 |
| present, another genome | 0.9879 (0.973, 0.993) | 0.37 | 0.9934 (0.979, 0.997) | 0.18 |
| present, in-silico strain | 0.9888 (0.975, 0.993) | 0.32 | 0.9939 (0.981, 0.997) | 0.16 |
| absent, beside a held-out congener | 0.9361 (0.916, 0.964) | 0.25 | 0.9371 (0.917, 0.968) | 0.27 |

The sequencing error of the new Illumina model is about half ART's (a representative's own reads 0.4% from it
instead of 0.9%), so the strains moved out of the band the novel congeners occupy; se (the pairs' first mates) moved
less (strains 0.9898 → 0.9914, band 0.31 → 0.22); PacBio and Nanopore are unchanged (strains 0.997 and 0.959 in
both). Whether 0.4% is the right error for the user's instruments is a question for the simulator, not the model;
either way v17's Illumina gain is in part the reads.

[holdout_clouds.txt](holdout_clouds.txt) adds the kept congeners' distances: the nearest kept congener of a
database species is at a median k-mer distance of 0.09 on the marker genes (10% 0.029, 25% 0.046); 3.5% of the
rows' taxa have one within 0.02 and 27% within 0.05.

## 2. What the feature groups add

[ablate_v17.py](ablate_v17.py) refits v17's pe tables as the build trained them (gradient boosting, 250 rounds at
0.1, 63 leaves, balanced classes, scenario rows weighted 0.25; 2 threads, 8-20 min a variant), with species held out
(5 folds by taxon) and on the test table by the final model; the knob as the v15 report's ablation chose it (the
highest weighted held-out F1 over 0.30-0.85 if it gains 0.002 over 0.5). [ablate_results_pe.tsv](ablate_results_pe.tsv);
[ablate_split.py](ablate_split.py) ([output](ablate_split_pe.txt)) splits the test rows by design and scenario.

| variant | features | knob | species held out, F1 at 0.5 | test, pooled, at knob | design test | soil hold-out | shallow soil hold-out | moderate | gut |
|---|---|---|---|---|---|---|---|---|---|
| v17, all | 81 | 0.5 | 0.9759 | 0.9761 | 0.9838 | 0.9719 | 0.9625 | 0.9867 | 0.9914 |
| without ancestry | 78 | 0.71 | 0.9755 | 0.9773 | 0.9830 | 0.9750 | 0.9631 | 0.9879 | 0.9910 |
| without gaps and untried | 76 | 0.5 | 0.9758 | 0.9757 | 0.9838 | 0.9715 | 0.9619 | 0.9861 | 0.9907 |
| **without foreign** | 78 | 0.74 | **0.9476** | **0.9545** | **0.9621** | **0.9546** | **0.9274** | 0.9700 | 0.9871 |
| without all four new groups | 70 | 0.72 | 0.9458 | 0.9526 | 0.9605 | 0.9517 | 0.9257 | 0.9687 | 0.9868 |
| the v15 set | 55 | 0.77 | 0.9453 | 0.9512 | 0.9584 | 0.9508 | 0.9229 | 0.9684 | 0.9858 |

(v17 as built: design test 0.9833 at 0.73, soil hold-out 0.9745, shallow 0.9617; the refit reproduces it. v15 as
built, on its own rows: 0.9645, 0.9655, 0.9396.)

- **The foreign group carries the jump.** Taking it out costs 0.022 of design-test F1 (0.9838 → 0.9621), 0.017 on the
  soil hold-out and 0.035 on shallow soil, with species held out 0.028. What it carries is the simulation pool
  ([section 3](#3-the-errors-and-the-ancestry-patterns), [foreign_leak.txt](foreign_leak.txt)): without it v17's
  model stands at 0.962 on its design test and 0.955 / 0.927 on the soils, the level of v15 (0.9645 / 0.9655 /
  0.9396 on v15's rows, which were easier in the scenarios: 6 samples of one richness and depth each, no moderate
  scenario) and of v17's own 55-feature refit (0.9584 / 0.9508 / 0.9229).
- **The ancestry sites are worth +0.0005** (species held out 0.9755 → 0.9759, design test 0.9830 → 0.9838) and the
  gaps +0.0004: real in direction, within noise in size on these rows. All four new groups over the v15 set without
  the foreign group: +0.002 (species held out 0.9453 → 0.9476, design test 0.9584 → 0.9621, shallow soil 0.9229 →
  0.9274).
- So the honest v17 is a model at v15's level with cleaner Illumina reads, a harder set of scenarios, and two
  small new features. The foreign group is out of the default set since `13188c2` (the other session's fix); its
  redesign (foreign sources from every genome's marker copies, every copy defined) is what would make it a feature
  rather than a label.

## 3. The errors and the ancestry patterns

[ancestry_patterns.py](ancestry_patterns.py) ([output](ancestry_patterns.txt)) joins each build table with the
model's calls (`model_logs/trained_model*.calls.tsv.gz`: training rows with species held out, test rows by the final
model, at the knob) and classes the rows by the collector's meta columns, as the v15 analysis did.

**The budget at the knob** (pe 0.73; the test rows here include the scenarios' hold-out samples, so their F1 is
below the console's design-only 0.9833):

| pe | rows | present | TP | FP | FN | F1 | FP: novel congener / near novel / other | FN: strain / in-silico / representative |
|---|---|---|---|---|---|---|---|---|
| training, species held out | 394,309 | 75,512 | 73,625 | 1,456 | 1,887 | 0.9778 | 1,265 / 115 / 76 | 1,270 / 263 / 354 |
| test, final model | 165,137 | 28,912 | 28,078 | 508 | 834 | 0.9767 | 451 / 40 / 17 | 572 / 131 / 131 |

The FP have a median 3 fragments (35% one, 28% ten or more) at identity 0.977; the FN 2 fragments (41% one, 21% ten
or more) at 0.962; 74% of the FP and 73% of the FN are in the two soil scenarios, 13% and 15% in the moderate one.
se, pb and ont are alike (FP 86-91% novel congeners; FN 67-74% real strains).

**The ancestry sites by class** (rows with 5 or more fragments, pe; the other read types are alike, with more sites
per record on long reads):

| | present rep | present strain | present in-silico | absent, novel congener | absent, near novel | absent |
|---|---|---|---|---|---|---|
| sites per record | 12.1 | 11.5 | 13.0 | 8.7 | 5.7 | 5.3 |
| agreement (the species' base) | 0.996 | 0.978 | 0.989 | 0.676 | 0.650 | 0.680 |
| the congener's base | 0.001 | 0.010 | 0.002 | 0.244 | 0.138 | 0.221 |
| identity | 0.996 | 0.993 | 0.994 | 0.939 | 0.925 | 0.945 |

A present taxon's reads agree with its reference at 98-100% of the sites where it differs from its nearest congener;
a novel congener's at 68%, and they carry the congener's base at 24% (the rest is neither: the novel species' own
states, and errors). 93% of the rows with 3 or more fragments have sites (97% of the novel-congener rows, 72% of
the "near novel" ones, whose closest relative is farther than a genus). The agreement of a novel congener rises as
the kept congener nearest the reference gets nearer (0.75 within 0.02, 0.66 beyond 0.1): the farther the reference's
own nearest congener, the more of its derived states the novel species shares.

**Inside the identity band 0.95-0.985**, which holds 64% of pe's FP and 52% of its FN, AUC for present (strains and
representatives) against absent novel congener, real strains only:

| feature | whole band | 0.95-0.96 | 0.96-0.97 | 0.97-0.98 | 0.98-0.985 |
|---|---|---|---|---|---|
| ancestry_agreement | 0.818 | 0.665 | 0.685 | 0.732 | 0.755 |
| ancestry_sites_per_record | 0.677 | 0.600 | 0.667 | 0.748 | 0.800 |
| ancestry_congener_share (novel high) | 0.816 | 0.694 | 0.712 | 0.743 | 0.778 |
| lu_per_kb | 0.853 | 0.757 | 0.770 | 0.803 | 0.849 |
| em_own_share | 0.769 | 0.666 | 0.723 | 0.800 | 0.852 |
| identity | 0.865 | 0.624 | 0.641 | 0.638 | 0.587 |
| foreign_scanned_share (the leak) | 0.898 | 0.907 | 0.900 | 0.885 | 0.887 |
| gap_within_min_share | 0.692 | 0.572 | 0.620 | 0.675 | 0.688 |
| gap_position (novel high) | 0.822 | 0.665 | 0.725 | 0.741 | 0.741 |

A boosted classifier inside the band, samples grouped: base features 0.9593 AUC; with the ancestry sites 0.9650 (AP
0.917 → 0.929); with the foreign features 0.9917; with both 0.9923 (se 0.967 → 0.972 → 0.993; pb 0.934 → 0.947 →
0.983; ont 0.986 → 0.990 → 0.997). The ancestry sites add what unique k-mers do not; the foreign features add the
pool.

**At the residual errors** (rows with 3 or more fragments, pe): FN taxa have ancestry agreement 0.754 and the
congener's base at 0.190 of their sites; FP taxa 0.832 and 0.128; right calls 0.991 / 0.002 (TP) and 0.668 / 0.240
(TN). So the strains that are still missed are the ones whose reads side with the congener at a fifth of the
sites (identity 0.960, `excess_scaled_median` 0.029), and the congeners that still pass are the closest ones
(identity 0.972, agreement 0.83). As in v15, no feature separates the two residues (AUC FN vs FP: identity 0.25,
`ancestry_agreement` 0.30, `ancestry_congener_share` 0.68, `excess_scaled_median` 0.72, `gap_position` 0.70, all
pointing the "wrong" way: the missed strains are the more divergent organisms).

**The error reads' sources** ([signatures_v17.txt](signatures_v17.txt), [sensitivity_v17.txt](sensitivity_v17.txt);
the build's SAMs keep 20 fragments per taxon and reason, so the per-taxon views are the ones to read):

| | pe | se | pb | ont |
|---|---|---|---|---|
| FP taxa whose main source is a congener (v15) | 0.855 (0.803) | 0.856 (0.921) | 0.954 (0.813) | 0.941 (0.813) |
| ... from beyond the genus (v15) | 0.145 (0.197) | 0.144 (0.079) | 0.045 (0.187) | 0.058 (0.187) |
| FP taxa fed only by held-out species (v15) | 0.884 (0.842) | 0.875 (0.817) | 0.942 (0.890) | 0.956 (0.887) |
| FP taxa fed by an FN species of the sample (v15) | 17 (42) | 13 (62) | 6 (15) | 3 (20) |
| counted FP fragments from species not in the database | 0.975 | 0.973 | 0.951 | 0.987 |
| FN taxa: own reads that seeded on them and aligned (v15) | 0.998 (0.958) | 0.997 (0.986) | 0.975 (0.973) | 0.910 (0.835) |
| FN taxa losing half or more of their aligned own reads to a congener (v15) | 0.224 (0.225) | 0.263 (0.297) | 0.155 (0.150) | 0.147 (0.140) |

The beyond-genus class shrank for the long reads (19% → 5% of the FP taxa) and for pe (20% → 15%): the suspect
copies, the foreign rates and the gaps took the reads of other families. What remains is the novel congener, more
purely than in v15. The aligner keeps the missed strains' reads that seed on them (ONT 91%, up from 84%: the
adaptive candidates); a fifth to a quarter of the FN taxa still lose half of their aligned reads to a congener that
fits as well, as in v15.

## 4. The clouds around the held-out species

[holdout_clouds.py](holdout_clouds.py) ([output](holdout_clouds.txt)).

**How the hold-out is chosen** (`build_gtdb_database.choose_holdout`): whole clades first (2 phyla, 4 classes, 6
orders, 8 families, 12 genera, drawn among those with two or more species to simulate and at most 2% of the species),
then 30% of the remaining species the genome table can simulate, at random, the same share per domain. Nothing in
it looks at how near a species' congeners are. The genome table simulates 24,980 of the database's 143,614 species;
the other 118,634 are kept in every training database and never present in a sample.

**What the taxonomy says** (v17's held-out list is v15's):
- Of the 7,324 species held out alone, 91.5% have a congener kept in the training database (73.5% five or more):
  their reads have a congener to land on, the FP class. 8.5% have none, and of those 14% have no kept species in
  their family either: their reads land farther away or nowhere.
- Of the 133,014 kept species, 53% have a held-out congener (4.5% in genera of 2 species, 13% in genera of 3-5, 41% in
  6-20, 91% in 21 or more); 442 kept species are the last of their genus.
- The clade hold-outs (3,276 species) have no kept species at their clade's rank by construction: deeper novelty.

**What the distances say** (v17's tables; the reads' identity stands in for the held-out species' distance to the
kept reference it lands on):
- Absent rows beside a held-out congener, 3 or more fragments: identity median 0.937; 2.1% at 0.985 or above, 6.5%
  at 0.97-0.985. Present strains: 82% at 0.985 or above, 15% at 0.97-0.985.
- So the hold-out tears a near-identical cloud (a sister within a strain's distance of its kept twin) in a few percent
  of the cases, and those are exactly where the FP remain (identity 0.977): at 0.98 the model is taught "absent" by
  a held-out twin and "present" by a divergent strain, and no feature can satisfy both.
- The kept clouds are dense: 27% of the database species have a kept congener within 0.05, 3.5% within 0.02
  (`db_nearest_congener`).

**Answer.** The hold-out does not break the clouds systematically: a held-out species keeps its congeners with the
same probability as any species, so the training sees novel congeners at every distance, and the near-identical
twins are a few percent. It also does not *use* the clouds: it cannot tell the trainer that a given absent row is a
twin's, nor place held-out species by distance band, nor keep species complexes together. protal already records
the clouds, but after the hold-out is chosen: `species_neighbours.tsv` (each species' 16 nearest congeners within
0.15, from the references) and `congener_gaps.tsv` (per gene copy its nearest congener by alignment, the min and
median gaps) are written by `--build` into each database. What would let them steer the simulations:

1. **A `species_clouds.tsv` from the converted references before the hold-out** (the gaps pass on
   `reference.fna`, which took the training database 4 minutes at r226 with the foreign rates): per species its
   nearest congener and the distance, the congeners within 0.01, 0.02 and 0.05, and the within-species spread where
   the full reference has other genomes.
2. **`choose_holdout` reads it**: species complexes (a congener within 0.01 on most genes) are held out together or
   kept together, never torn; the held-out species are reported by their nearest kept congener's distance, and
   `--holdout-bands` could demand a share in each band (0-0.02, 0.02-0.05, 0.05-0.15, farther).
3. **`collect_training_data.py` writes `meta_novel_distance`** for an absent taxon beside a held-out congener: the
   gap between the two (from the finished database's table, which knows the held-out species), so that the trainer
   reports the FP rate by distance band, weights the twins' rows down (a model should not be punished for calling a
   species' twin), and the genus-level fallback has the number it needs.
4. **The simulator's congener groups by distance**: `--congeners 0.25:2-5` draws genus members at random; drawing
   them by band (near, mid, far) makes every sample hold the three situations the model must tell apart.

The last two are where the gain is; the first two make them possible before the database exists.

## 5. Implemented afterwards (2026-10-08, uncommitted at the time of writing)

Items 1-3 above and the scan without the leak, for the next r226 build to measure:

- **Species clouds before the hold-out.** `protal --write_species_neighbours FILE --db FOLDER` writes the species
  neighbours table (what `--build` stores) from a converted folder without building; `build_gtdb_database.py` runs it
  on the converted release (`species_clouds.tsv`, also in `model_logs/`), or takes `--species-clouds FILE`.
- **Complexes held out whole.** `choose_holdout` joins congeners within `--holdout-complex-distance` (0.01) into
  complexes (union-find on the clouds) and draws units, a complex or a species alone, until the domain's share is
  reached; `heldout_species.txt` gains the distance to the nearest kept congener and the complex, `holdout.txt` the
  bands. Without complexes the draw is the old one (same seed, same species).
- **`meta_novel_distance`.** `collect_training_data.py --species_clouds` writes, for a taxon whose genus has a
  held-out species in the sample, the distance to the nearest one; the trainer reports FP and FN by distance band and
  takes `--twin-weight` / `--twin-distance` (1 and 0.01 by default: reported, not weighted).
- **The foreign scan from the full reference.** `simulate_metagenomes --tiles L:S --tile_fasta full_reference.fna.zst
  --tile_per_header 10` tiles every species' marker copies alike and lists every header's records and tiles;
  `scripts/foreign_rates.py --full-reference` counts the reads and lists every copy, reached or not; the profiler
  counts a listed copy with no read as unreached. A copy's reads are its foreign reads plus one genome's worth of its
  own (own reads over the records tiled), not every own read: those scale with the species' genome count, and the
  share would have carried the cluster-size rule back in (the other session's point). The build scans the training database (its full reference lacks the
  held-out species) and, before the models go in, the finished one; on by default, `--no-foreign-rates` skips it.
  `foreign` stays out of the default set; a named set with it is among the `--features auto` candidates.
- Not done: item 4, the simulator's congener groups by distance band (the designer draws genus members at random;
  the clouds analysis shows every band is represented, so the trainer's report by band comes first).
- **Consensus ancestry sites** (the user's suggestion, after section 3): protal compares each copy with its congeners
  (the gaps' nearest and the 16 species neighbours, one vote each); where three or more were compared at a position
  it is a site when nine in ten carry one base other than the species' (all of three to nine), with fewer the
  nearest congener's difference as before. The build's ancestry report gains the same definition
  (`species_base_at_consensus_sites`) beside the nearest-congener and majority-of-three sites, so the next build
  shows all three on its error reads.
- **Chained comparison and indel sites** (the user's follow-up): the paired 12-mers are chained across changes of
  diagonal of up to 60 bases, so the stretches past an indel are compared (before, the comparison ended at the first
  indel) and the indel is located by extending the exact matches inwards from both anchors. An indel of three bases
  or more that the congeners share against the species (the same consensus rule) is an indel site; a record aligned
  five bases beyond it counts the congeners' state when it has a gap of that kind and length within four bases of it,
  the species' when it has no gap near. Two features, `ancestry_indel_sites_per_record` and
  `ancestry_indel_congener_share` (80 in the default set). The report script (`ancestry_sites.py`) does not chain
  yet: its indel mirror is left for later. On the mini database (its copies carry codon indels) the base sites of
  Mockella alpha went from 10.8 per record at 0.957 agreement to 9.2 at 0.929: the old comparison filled the stretch
  between two main-diagonal anchors that straddled an indel and its compensating indel with misaligned bases, so
  it counted spurious sites there, which a strain's read "agreed" with; the chain compares that stretch on its own
  diagonals. Unit 459 (two test expectations corrected), Python 53, e2e 142, pipeline test_a pass.
- **What the next build said** ([2026-10-08-r226-v18](../2026-10-08-r226-v18/README.md), the other session's):
  the consensus, chained and indel sites separate the classes no better than v17's nearest-congener sites (band
  agreement AUC pe 0.815 against 0.837), the indel features are empty for short reads, the leak-free foreign scan has
  no leak and adds 0.001 AUC. And the ancestry report's fixed-site AUC (0.81-0.86) was circular: a missed real
  strain's own GTDB genome was usually among the six alleles. Fixed afterwards: the converter names each
  full-reference record's genome (`>taxid_geneid accession`), and `ancestry_sites.py` leaves a read's own source
  genome (its `xg` tag) out of the alleles per record, so the next build's report gives the honest number.

Tests: `test_foreign_rates.py` rewritten, `test_gtdb_build.py` gains the complex hold-out, the pipeline's `test_a`
asserts both scans and the steered hold-out, `test_f` runs `--no-foreign-rates`; `tests/test_CongenerGaps.cpp` the
unreached listed copy.
