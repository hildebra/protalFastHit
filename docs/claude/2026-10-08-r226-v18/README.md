# The r226 v18 build: F1 across the builds, and what the ancestry sites do

**Data.** The r226 v18 build (SLURM 24089809, q1536n4, 2026-10-08 14:52 to 19:16, 84 threads). Its protal,
simulator and scripts were at `4b13640` at the start and `e24dcff` at the end: the HPC checkout moved during the run,
and the build's warning says so. Its share archive `local/protal0.7.9_r226_v18_share.tar.gz` (517 MB) is unpacked to
`local/v18`. The earlier builds come from `local/v10` to `local/v17`; v16 was never run to the end. The scripts here
read `audit-fixes` at `e24dcff`. They ran in WSL on one pinned core, niced, in `~/soil13/venv`.

    python3 f1_history.py --builds local/v10,local/v11,local/v12,local/v13,local/v14,local/v15,local/v17,local/v18 > f1_history.txt
    python3 ancestry_v18.py --build local/v18 --v17 local/v17 --types pe > ancestry_v18_pe.txt   # and se, pb, ont
    python3 fixed_sites_check.py --build local/v18 > fixed_sites_check.txt
    python3 ../2026-10-08-r226-v17/foreign_leak.py --build local/v18 > foreign_leak_v18.txt

**Questions (the user's, 2026-10-08).** The build was "partially successful": analyse the results and how F1 changed
over the last builds. Then: did the ancestry patterns help, and should anything change in how protal uses them?

## Summary

- **v18 is at the honest level of v14 and v15, not below it.** Paired-end test F1 is 0.962 (v15 0.9645, v14 0.9649).
  v17's 0.983 was the foreign-copy leak ([2026-10-07-congener-gaps](../2026-10-07-congener-gaps/README.md#the-foreign-features-leak-2026-10-08)).
  The v17 report's refit without the leaking group gave 0.9621 on its design test, 0.9546 on the soil hold-out and
  0.9274 on shallow soil. v18 gives 0.9622, 0.9541 and 0.9269. The step down from v17 is the leak's removal, not a
  regression.
- **Long reads gained over v15, without the leak.** PacBio: test 0.9785 → 0.9816, soil hold-out 0.958 → 0.968, shallow
  soil 0.956 → 0.963. Nanopore: test 0.972 → 0.979, soil 0.948 → 0.959. This is where the ancestry features rank
  highest (5th-8th in the long-read models). Short reads are flat on the test set; their scenarios have been harder
  since v17, so they read slightly lower.
- **The ancestry sites help a little, as in v17.** Inside the identity band 0.95-0.985, a boosted classifier gains
  AUC +0.0068 (pe), +0.0043 (se), +0.0118 (pb) and +0.0033 (ont) over the base features. v17's pe ablation put that
  at +0.0005 of F1. The trainer ranks `ancestry_agreement` 12th (pe, se), 8th (pb) and 5th (ont).
- **v18's site changes did not add to that.** v18 changed the sites three ways: consensus over the congeners, the
  comparison chained across indels, and two indel-site features.
  - They cut the sites per record by a third to two thirds, and the novel congeners' agreement fell (pe 0.676 → 0.605).
  - The per-feature AUC inside the band fell slightly for pe, se and pb (agreement 0.837 → 0.815, 0.847 → 0.826,
    0.805 → 0.786) and stayed for ont. The classifier's gain is the same as v17's.
  - The indel sites are 0 in the median row of every class for short reads. For long reads they are rare
    (0.03-0.07 per record on present taxa), and they add 0.0004-0.001 of AUC.
- **The ancestry report's headline is circular.** The fixed sites' AUC of 0.79-0.86 (FN strains against FP
  congeners) rests on the missed real strains' own genomes being among the "alleles" that define the polymorphic
  sites.
  - Of their mismatches, 80-92% fall at polymorphic sites that cover 3.5-3.9% of their bases. In 52-64% of their
    records with two or more mismatches, every mismatch does (short reads, PacBio).
  - The in-silico strains, whose genomes the full reference does not hold, gain nothing from the fixed sites.
  - The honest value is unknown: it needs a leave-one-genome-out measurement ([section 3](#3-the-fixed-sites-are-circular-in-the-report)).
- **The residual errors are where they were.** At pe's knob, 76% of the FP are novel congeners and 69% of the FN real
  strains, two thirds of both in the soils. On the ancestry features the FN look more like congeners than the FP do
  (agreement AUC FN vs FP 0.29-0.33 across read types), as in v17.
- **The new foreign scan has no leak, and little signal.** In-pool and out-of-pool species now have the same scanned
  share (0.97 vs 0.98; v17: 1.0 vs 0.0), and the leak check's agreement with pool membership fell from 0.87-0.93 to
  0.30-0.58. In the band it adds AUC +0.001. It is in no model (not in the default set), yet the build
  ran the scan twice, for 21 minutes.

## 1. F1 over the builds

[f1_history.txt](f1_history.txt): F1 with species held out at 0.5, on the independent test set and on the
scenarios' hold-out samples at the model's knob. The rows of each build are its own: simulations, hold-outs and
scenarios differ, so this compares builds, not models.

| pe | v10 | v11 | v12 | v13 | v14 | v15 | v17 | v18 |
|---|---|---|---|---|---|---|---|---|
| species held out | 0.9643 | 0.9598 | 0.9597 | 0.9496 | 0.9519 | 0.9519 | *0.9760* | 0.9462 |
| test set | 0.9630 | 0.9614 | 0.9701 | 0.9652 | 0.9649 | 0.9645 | *0.9833* | 0.9622 |
| soil hold-out | | | 0.9594 | 0.9535 | 0.9646 | 0.9655 | *0.9745* | 0.9541 |
| shallow soil hold-out | | | 0.9371 | 0.9253 | 0.9408 | 0.9396 | *0.9617* | 0.9269 |
| gut hold-out | | | 0.9954 | 0.9950 | 0.9958 | 0.9962 | *0.9910* | 0.9795 |
| knob | 0.5 | 0.5 | 0.5 | 0.7 | 0.82 | 0.8 | 0.73 | 0.74 |

| test set / soil hold-out | v12 | v13 | v14 | v15 | v17 | v18 |
|---|---|---|---|---|---|---|
| se | 0.9715 / 0.9671 | 0.9592 / 0.9511 | 0.9612 / 0.9645 | 0.9622 / 0.9635 | *0.9836 / 0.9758* | 0.9599 / 0.9517 |
| pb | 0.9807 / 0.9608 | 0.9767 / 0.9469 | 0.9791 / 0.9577 | 0.9785 / 0.9578 | *0.9894 / 0.9835* | 0.9816 / 0.9684 |
| ont | 0.9738 / 0.9533 | 0.9720 / 0.9412 | 0.9759 / 0.9439 | 0.9723 / 0.9483 | *0.9886 / 0.9795* | 0.9794 / 0.9586 |

(v17 in italics: inflated by the foreign leak. Its refit without the foreign group, from the v17 report: pe species
held out 0.9476, design test 0.9621, soil 0.9546, shallow 0.9274, moderate 0.9700, gut 0.9871.)

What moved between the builds, from their metadata and reports:
- v10-v11 (0.7.5): six groups.
- v12 (0.7.6): scenarios, `--features auto`.
- v13 (0.7.7): gradient boosting.
- v14: the `ref` group, scenarios at varied depths.
- v15: the sample's complexity.
- v17 (0.7.9 + congener gaps): the fp-feature groups, ancestry, gaps, foreign (leaking) and untried. Also new Illumina
  reads (no ART), scenarios varied in richness and depth plus a moderate scenario, the bootstrap knob, strains 0.5,0.2.
- v18: foreign out of the default set (its scan without the leak, unused), the consensus, chained and indel ancestry
  sites, the hold-out steered by species clouds (102 complexes of 247 species held out or kept whole), and
  `meta_novel_distance`.

Reading it:
- **pe.** The test set has held at 0.962-0.970 since v12. v18 is level with v14 and v15. Its species-held-out F1
  (0.946) is the lowest yet, but v18's hold-out is a new draw (complexes kept whole), so it scores other species.
- **The soils** (0.954, 0.927) sit below v15 (0.9655, 0.9396) and level with v17 without foreign (0.9546, 0.9274).
  Since v17 the soil samples vary in richness and depth, which makes them harder. v18 matches v17's honest refit
  on every scenario but gut (0.9795 vs 0.9871, in 4 hold-out samples; 50 FN, where v17 as built had 19).
- **More FP at the same FN.** pe test FP rose from 64 (v15) to 106 at about the same FN (253 → 243), with the knob at
  0.74 against v15's 0.8. FP per sample is still 1.4.
- **The long reads.** PacBio and Nanopore are the honest gain: +0.003 and +0.007 on the test set, +0.010 on soil
  against v15. v17's report found their reads unchanged since v15, so the gain is from the models' features. The
  ancestry sites are the candidate: they rank 5th-8th there, against 12th for short reads.

## 2. What the ancestry features do in v18

[ancestry_v18_pe.txt](ancestry_v18_pe.txt), [_se](ancestry_v18_se.txt), [_pb](ancestry_v18_pb.txt),
[_ont](ancestry_v18_ont.txt). They use the v17 report's `ancestry_patterns.py` for the calls and the classes of rows.

**By class**, medians of rows with 5 or more fragments, v18 (v17):

| | present rep | present strain | present in-silico | absent novel congener |
|---|---|---|---|---|
| pe sites per record | 8.4 (12.1) | 7.5 (11.5) | 10.8 (13.0) | 2.9 (8.7) |
| pe agreement | 0.996 (0.996) | 0.976 (0.978) | 0.990 (0.989) | 0.605 (0.676) |
| pe congener's base | 0.000 (0.001) | 0.010 (0.010) | 0.002 (0.002) | 0.311 (0.244) |
| pb sites per record | 61.6 (87.9) | 55.7 (85.8) | 77.2 (95.5) | 28.9 (66.4) |
| pb agreement | 1.000 (1.000) | 0.982 (0.984) | 0.995 (0.995) | 0.590 (0.657) |
| pe / pb indel sites per record | 0 / 0.051 | 0 / 0.044 | 0 / 0.069 | 0 / 0 |

The consensus sites keep the positions where the congeners agree, and the chain drops the spurious sites between
anchors that straddled an indel. Both separate the classes further: novel congeners agree less and carry the
congener's base more. They also cost evidence, with a third to two thirds fewer sites per record and nearly none on
a novel congener's reads.

**Inside the identity band 0.95-0.985**, rows with 3 or more fragments, present (high) against absent novel congener:

| AUC | pe v18 (v17) | se v18 (v17) | pb v18 (v17) | ont v18 (v17) |
|---|---|---|---|---|
| ancestry_agreement | 0.815 (0.837) | 0.826 (0.847) | 0.786 (0.805) | 0.887 (0.884) |
| ancestry_sites_per_record | 0.695 (0.689) | 0.712 (0.700) | 0.642 (0.639) | 0.747 (0.743) |
| ancestry_indel_sites_per_record | 0.618 | 0.631 | 0.611 | 0.702 |
| identity | 0.871 | 0.876 | 0.847 | 0.649 |
| lu_per_kb | 0.878 | 0.890 | 0.835 | 0.911 |

| boosted classifier in the band, samples grouped, AUC | pe | se | pb | ont |
|---|---|---|---|---|
| base features | 0.9576 | 0.9679 | 0.9436 | 0.9854 |
| + v17's three ancestry features | 0.9640 | 0.9726 | 0.9543 | 0.9883 |
| + all five (with the indel sites) | 0.9644 | 0.9722 | 0.9554 | 0.9887 |
| + foreign (the scan without the leak) | 0.9587 | 0.9685 | 0.9449 | 0.9858 |
| v17: base → + ancestry | 0.9593 → 0.9650 | 0.967 → 0.972 | 0.934 → 0.947 | 0.986 → 0.990 |

- **The ancestry group earns its place.** It gains most for PacBio (+0.012 AUC in the band), the model where it
  ranks 7th-8th, and least for ont, whose band holds few novel congeners (2,057 of 42,078 rows).
- **v18's definition is not better than v17's.** Agreement is a little weaker alone in three read types, and the
  classifier's gain is the same.
- **The two indel features add 0.000-0.001.** For short reads they are almost always empty.

**At the errors** (rows with 3 or more fragments), the AUC of FN against FP, FN high:
- `ancestry_agreement`: 0.31 / 0.33 / 0.29 / 0.30 (pe / se / pb / ont);
- `identity`: 0.21 / 0.27 / 0.22 / 0.34;
- `excess_scaled_median`: 0.77 / 0.72 / 0.80 / 0.79.

The missed strains are the more divergent organisms: their reads side with the congener at about a fifth of the
sites (pe FN agreement 0.734, FP 0.854). No site definition changes that, because the reads behind the two errors
look alike. That is v15's and v17's finding again.

**The errors by scenario** (pe):
- In the test set, 67% of the FP and 69% of the FN are in the two soils; the design test has 106 FP and 243 FN.
- Of the FN (training and test), 67% are real strains, 17% in-silico strains and 16% representatives.
- The species complexes are no source: 0.5-1.5% of the FN have a database congener within 0.01, against 0.2-0.3% of
  the TP.

## 3. The fixed sites are circular in the report

**How the report defines them.** The build's ancestry report (`scripts/ancestry_sites.py`) takes up to 6 other
genomes' copies of each species' gene from the full reference as "alleles". A position is polymorphic where any of
them differs from the representative, and a congener site is fixed where it is not polymorphic. Its console line
ranks the fixed sites first: AUC 0.809 (pe), 0.839 (se), 0.860 (pb) and 0.814 (ont) for the FN taxa's own reads
against the FP reads.

**Why that is circular.** A missed real strain was simulated from a GTDB genome of the species, and that genome is
usually one of the 6. Its differences from the representative are then polymorphic by construction.
[fixed_sites_check.py](fixed_sites_check.py) ([output](fixed_sites_check.txt)) splits the FN taxa's own records by
what was simulated (the tables' `meta_rep_genome` and `meta_insilico_strain`):

| pe | FP genus | real strain | in-silico strain | representative |
|---|---|---|---|---|
| records with polymorphic sites | 0.333 | 0.830 | 0.035 | 0.385 |
| bases at polymorphic sites | 0.012 | 0.035 | 0.001 | 0.013 |
| mismatches at polymorphic sites | 0.146 | **0.798** | 0.003 | 0.019 |
| records whose every mismatch (≥ 2) is polymorphic | 0.119 | **0.520** | 0.036 | 0.008 |
| species' base at congener sites | 0.914 | 0.899 | 0.964 | 0.995 |
| species' base at fixed sites | 0.928 | 0.993 | 0.965 | 0.995 |

| AUC against FP genus, per taxon (≥ 10 sites) | congener sites | fixed sites |
|---|---|---|
| real strain, pe / se / pb / ont | 0.46 / 0.48 / 0.36 / 0.39 | 0.83 / 0.86 / 0.88 / 0.83 |
| in-silico strain, pe / se / pb / ont | 0.75 / 0.77 / 0.72 / 0.76 | 0.68 / 0.71 / 0.63 / 0.68 |

- **The real strains match their own genome.** 80% of their mismatches (88% se, 92% pb) fall at polymorphic sites
  that cover 3.5% of their bases. Half of their records with two or more mismatches have every one of them there,
  against 12% for the FP. That is the strain's own genome in the allele set. ONT's errors dilute it (52%).
- **The in-silico strains get nothing.** Their genomes are not in the full reference, and their species have one
  genome, so no alleles. The fixed sites are then the congener sites, and they gain nothing.
- **So 0.83-0.88 is an upper bound, and the in-silico strains are no lower bound either**: they have no alleles at
  all. A novel strain of a species with known strains lies in between.
  - The FP genus records show that between-species differences are enriched 12-fold at the polymorphic sites (15% of
    mismatches on 1.2% of the bases).
  - A strain's few differences should fall on the fast sites even more often. How much of the lift survives is
    unknown.

## 4. Recommendations for the ancestry sites in protal

1. **Keep the ancestry group in the default set.**
   - It is the third read-level signal inside the band after `lu_per_kb` and identity, and it gains most for the long
     reads, which are v18's honest gain.
   - Its F1 worth at v18 has not been measured. v17's pe ablation gave +0.0005 with the leaking group masking it. A
     v18 refit with and without it, per read type, would put a number on each; see the end.
2. **Leave the consensus sites as they are, or go back to the nearest congener's: it does not matter.** The report's
   three definitions (nearest, majority of 3, consensus) give AUC 0.549-0.564 at the residual errors, and the run-time
   features give the same classifier gain. If anything changes, weight a site by the share of congeners carrying the
   other base rather than dropping it below 0.9. That keeps the evidence the threshold throws away on reads with few
   sites, which are 40% of the errors.
3. **The two indel features can stay, but expect nothing from them.** They are empty for short reads and add at most
   0.001 AUC for long reads. The chained comparison is worth keeping for the spurious sites it removes.
4. **Do not build the polymorphic-site mask into protal until it is measured without the circularity.** The
   measurement:
   - The converter writes each full-reference record's genome accession (it has it when it writes the record, in
     `gtdb_to_protal_db.py`: a sidecar or a header suffix).
   - `ancestry_sites.py` drops the simulated genome's copies from a read's alleles. The read's source genome comes
     from its contig (the build's `genome_contigs.tsv.gz`).
   - The next build's report then gives the honest AUC for real strains.

   Only if that clearly beats identity's AUC inside the bands should protal store a per-copy mask. And then the
   training database's mask must leave out every genome the simulations draw strains from, as the strain-alleles plan
   noted, or the models learn the circularity the report shows here.
5. **Make the foreign scan opt-in again, or put the group in the default set if an ablation shows a gain.** It has no
   leak and little signal: +0.001 AUC in the band. As built it costs 21 minutes (two scans) for features no model
   uses.

## 5. Strain alleles in the index: not in v18

The user asked whether v18 adds alleles to each reference as differences to the reference gene, and whether that helped.
It does not: protal's index holds the representative's copy of each marker gene only, at `4b13640` and on every branch.
The build reads the other genomes' copies (`full_reference.fna.zst`) only for the k-mer uniqueness check, the
conservation factors, the suspect copies, the foreign scan and the ancestry report. Strain alleles in the index were
proposed in the v7/v8, v9, v13 and v14 reports and in the error-read report (items 2 and 4, "for re-scoring, not for
seeding"), but never built. v18's F1 owes them nothing.

What they would reach in v18 (the FN by what was simulated, training and test rows,
[ancestry_v18_*.txt](ancestry_v18_pe.txt), [fixed_sites_check.txt](fixed_sites_check.txt)):

| FN | pe | se | pb | ont |
|---|---|---|---|---|
| real strains (a GTDB genome other than the representative) | 3,834 (67%) | 1,915 (70%) | 2,081 (82%) | 2,249 (77%) |
| in-silico strains (one-genome species: no allele exists) | 951 (17%) | 442 (16%) | 252 (10%) | 319 (11%) |
| representatives (alleles change nothing) | 929 (16%) | 393 (14%) | 217 (9%) | 351 (12%) |
| real-strain FN taxa whose species has other genomes (>= 10 sites) | 90% | 91% | 98% | 96% |

- **What they could reach:** about 60-80% of the FN, the real strains of species with other genomes. The ceilings the
  earlier reports gave (+0.009 overall at v9, +0.011-0.013 in deep soil at v14) are of that order.
- **The same circularity as the fixed sites.** Every simulated real strain is a GTDB genome, so an index with every
  genome's alleles would hold each simulated strain exactly. Its reads would align at 0.99+, and the models would be
  trained and tested on strains the database already knows. In use, a sample's strain is only as close as GTDB's
  nearest genome of the species.
- **So the training database's alleles must leave out every genome the simulations draw strains from**, as for the
  fixed-site mask. Only the shipped database would hold them all.
- **How much is honestly gained is unknown.** Next build's leave-one-genome-out ancestry report (`1c184fb`) measures
  the nearest piece of it: how many of a missed strain's differences the species' *other* genomes share. That is the
  number to read before building alleles into the index.

## What is not done

- **No refit ablation of v18.** The F1 worth of the ancestry group, its indel features and the leak-free foreign
  group, per read type, would come from refitting v18's tables with and without them, as `ablate_v17.py` does (pe ~1
  hour on 2 cores, the long reads ~20 minutes). The user's rule is to ask before local runs.
- **The gut hold-out drop** (pe 0.9795 with 50 FN; v17 as built had 19) is unexamined: 4 samples, and a new hold-out draw.
- **The leave-one-genome-out fixed-site measurement** (recommendation 4) needs the next build.
