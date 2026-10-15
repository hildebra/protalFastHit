# protal since 0.6.0a: what the changes can do to F1

2026-10-15. Branch `audit-fixes` at `dd01eb1`, against the tag `0.6.0a` (`014f4a9`, 2026-08-10): 463 commits,
`src/` +51,684 / −11,328 lines in 159 files.

**The question** (the user's): audit protal, specifically the changes introduced since 0.6.0a and how they could
affect F1.

**How it was done.**
- **Code:** the path a read takes to a call, read in five parts: read input and sample definition; index and
  seeding; alignment, candidate choice and MAPQ; profiling features and the model; training data and evaluation.
  Each part was read against `git diff 0.6.0a..dd01eb1` by its own reviewer. The findings rated medium or higher
  were then checked against the code again for this report; the others are marked "read once".
- **Measured effects:** taken from the earlier reports in this folder, linked where used. Nothing was built or run
  for this report, so no number below is new.
- **Line numbers:** `src/` at `dd01eb1`. Another session was editing `scripts/*.py` and `docs/databases.md` in the
  working tree during the audit (the fixes of the [GTDB build audit](2026-10-15-gtdb-build-audit/README.md)), so the
  script line numbers are those of `git show dd01eb1:<file>`, and some script findings may already be fixed there.
  This report does not repeat that audit's findings G1-G22.

## Summary

- **Almost all of the F1 gain since 0.6.0a came from one change, 0.7.0's model trained for its own database.** On the
  benchmark world 0.6.0a scores paired-end F1 0.634 and 0.7.0 0.903. Run with 0.6.0a's model, 0.7.0 detects as
  0.6.0a does. 0.7.1 to 0.7.5 added +0.012 together. At GTDB r226 the simulated test F1 has been flat at about
  0.95-0.97 since 0.7.4. 0.8.x's feature groups gained +0.001 to +0.004 each, much of it within noise
  ([section 1](#1-what-has-been-measured)).
- **The speed work did not change what is called.** The exact shortcuts give the same results: the k-mer screen, the
  indel bound, the anchored flanks, the WFA2 patches, the compact key map, the SIMD kernels and the parallel input
  decoders ([section 4](#4-checked-and-found-safe)). F1 is not at risk there.
- **The risks are somewhere else, in three places:**
  1. **What is evaluated is not what users run** (E1-E4):
     - the models are scored on the training database, never on the database that ships;
     - the test set uses the same held-out species as training;
     - the headline F1 counts only taxa that have a row;
     - "by species" cross-validation does not hold out the sample.

     So the reported F1 is optimistic by an amount nobody has measured.
  2. **Train/run agreement is checked by feature names only** (C1-C3). Most run-time defaults that change alignments
     are not recorded with the database. Read-type detection and placeholder base qualities can send a sample
     through the wrong model, or put its strongest feature far outside its training range, without a warning.
  3. **Two built-in biases** (B1, B2):
     - the strain-allele score shift only ever helps, so it favours species with many genomes;
     - the default sparse-first gates drop short and divergent reads with all their evidence.

     Both pay on the simulations, and both have an untested cost on real data.
- **Real data already shows a gap.** On twelve real mock samples (reduced r226 database, 761 species), per-sample F1
  is 0.85-0.97 and pooled F1 0.816 (0.927 counting only species with a read). The simulated test sets of the same
  database score 0.960-0.967. A few misses had plenty of reads; in simulation such species are never missed:
  - *Akkermansia muciniphila*: 1,112 reads on 38 genes at 2.3×, P 0.484;
  - *E. coli* in the same sample: 981 reads, P 0.375;
  - *Bacteroides ovatus*: 71 reads on 44 genes, P 0.469.

  The unknown share was 0.13-0.20 where 0.006 was expected. The long reads' 9 false positives in 4 samples were
  *Escherichia* and *Salmonella* congeners at low abundance ([section 1](#1-what-has-been-measured)).
- **The first checks to run, cheapest first** ([section 5](#5-what-to-check-next)):
  1. profile the test samples on the finished database;
  2. run the test set once with `--no_allele_scores`;
  3. record the alignment-side defaults in `build_settings.tsv` and stamp the feature code in the models;
  4. treat all-`!` qualities as missing, and warn when the long-read type was guessed from qualities alone.

## 1. What has been measured

**Between versions, on one simulated world** (900 species, 135 in no database; each version with its own models;
[version_changes.md](../version_changes.md#benchmarks)). Paired-end F1 at knob 0.5:

| version | full database | 290-349 species missing | source |
|---|---|---|---|
| 0.6.0a | 0.634 | 0.629 | [v0.7.5 benchmark](2026-10-03-v075-benchmark/README.md) |
| 0.7.0 | 0.903 | 0.895 | [v0.7.3 benchmark](2026-10-02-v073-benchmark/README.md) |
| 0.7.1 | 0.908 | 0.900 | same |
| 0.7.2 | 0.911 | 0.906 | same |
| 0.7.3 | 0.913 | 0.902 | [v0.7.5 benchmark](2026-10-03-v075-benchmark/README.md) |
| 0.7.5 | 0.915 | 0.902 | same |

No version after 0.7.5 was run on this world. 0.6.0a cannot read today's databases, and no real sample has been
profiled with both 0.6.0a and a later version.

**At GTDB r226, the build's own test sets.** Each build drew its own sample mix, so the F1 values are not comparable
across builds. The reports say so each time, and several "steps" turned out to be the mix:
- v24b's drop is a harder mix;
- most of v21's rise is an easier one;
- v17's jump was the `foreign` leak.

| build (version) | pe F1 | what it shows | source |
|---|---|---|---|
| v1 (0.7.1 era) | 0.948 test | first r226 models; FP-limited | [r226 build evaluation](2026-10-02-r226-build-evaluation/README.md) |
| v8 (0.7.4) | 0.9671 test | + sample depth feature | [v7/v8](2026-10-03-r226-v7-v8-evaluation/README.md) |
| v10 (0.7.5) | 0.9630 test | species priors off | [v10](2026-10-04-r226-v10-evaluation/README.md) |
| v17 (0.7.9) | 0.983 test | the `foreign` leak; 0.962 without it | [v17](2026-10-08-r226-v17/README.md) |
| v18 (0.7.9) | 0.962 test | the honest level | [v18](2026-10-08-r226-v18/README.md) |
| v21 (0.8.0) | 0.9548 held out | alleles, ancestry +0.001-0.002 each | [v21](2026-10-09-r226-v21/README.md) |
| v23 (0.8.1) | 0.9595 held out | site groups +0.001-0.003 together | [v23](2026-10-11-r226-v23/README.md) |
| v24b (0.8.2) | 0.9513 held out | harder, more realistic mix | [v24b](2026-10-12-r226-v24b/README.md) |
| v27 (0.8.3) | 0.9431 test, C 1 | gates up to pe C 4 cost no F1 | [v27 gates](2026-10-13-r226-v27-c-gates.md) |

**On real samples** (the reduced r226 database of 761 species, 0.8.4 and earlier):
- **[mock local test](2026-10-13-mock-local-test/README.md):** per-sample F1 0.85-0.97, three false positives, *E. coli*
  missed in both HiFi samples. That miss was fixed by `allele_top_identity`
  ([report](2026-10-13-consensus-allele-identity/README.md)).
- **[sparse gates on the mocks](2026-10-14-mock-sparse-gates/README.md):** pooled F1 0.816-0.831 over 10 samples.
  - 58 of the 80 misses had no read on their own genes: hCom2 strains lost in culture, and D6331's rarest species.
  - Counting only species with a read, F1 is 0.927.
  - The near-misses with real evidence are the ones in the summary, P 0.375-0.484.
- **The unknown share** was 0.13-0.20 in the D6331 pe samples against 0.006 expected, while on simulated samples it is
  within 0.02. The local test calls this "a real-read effect" and leaves it open. Its first candidate cause is that
  real strains' reads align less often than simulated ones.

The real-data numbers come from a reduced database, a handful of samples and one or two runs each. They show the
direction of the gap, not its size.

## 2. The changes since 0.6.0a, by stage, and what is known of their F1

"Default" means on in a default run on a database built by `build_gtdb_database.py` at `dd01eb1`. The "F1 evidence"
column says what was measured and where.

| stage | change (version) | default | F1 evidence |
|---|---|---|---|
| input | se, pb, ont read types, a model each (0.7.0) | yes | new capability; per type at r226 |
| input | read type guessed from the first 200 reads (0.7.3) | without `--read_type` | **never measured**: training always passes the type (C2) |
| input | ISA-L / BGZF / zstd / parallel gzip decoders, lanes as one sample (0.7.8-0.8.4) | yes | same reads, checked (section 4) |
| index | index format 2, correct unique flags; no k-mer with N (0.7.0); short-unique fix `78f7a92` | yes | in every r226 build |
| index | packed values, 2-bit genes, compact key map (0.7.1-0.8.4) | yes | same records |
| seeding | sparse-first gates; stride pick; pe:4, se:4, pb:16, ont:8 (0.8.3, 0.8.4) | GTDB builds | r226 v27: no F1 loss up to these gates; mocks: same calls (B2) |
| seeding | unique k-mers counted from the alignment, `unique_cores.bin` (after 0.8.4) | yes | strain world only ([report](2026-10-14-length-free-features/README.md)) |
| alignment | anchored flanks, k-mer screen, indel bound, WFA2 patches (0.7.3-0.8.2) | yes | exact; same alignments |
| alignment | Singletrack traceback (0.8.2) | `--build` | 0.01% pe, 47-49% ONT records differ among equal alignments; never evaluated alone |
| alignment | mates across neighbour genes, mate rescue (0.7.1) | yes | inside 0.7.1's +0.0055 |
| alignment | adaptive candidates (0.8.0) | yes | not measured alone |
| alignment | strain-allele scores, polymorphic half-differences (0.8.0) | yes, with alleles | strain world: +9/−2 strains; r226 shift **never measured** (B1) |
| alignment | weighted site shift (0.8.2 on, 0.8.3 off) | no | neutral to +0.003 (small real world); not recorded with the database (C1) |
| profiling | the presence model's features: 98 in the default set (0.7.0 → after 0.8.4) | yes | per group +0.001-0.004 at r226 |
| profiling | depth margin 0.08 (0.7.1) | yes | abundance only |
| profiling | suspect gene copies left out (0.7.5) | yes | 0.04% of copies at r226 |
| profiling | read EM, relatives, distance (0.7.4) | yes | in 0.7.5's +0.0019 (within noise) |
| model | knob curves over depth (0.7.3); model per gate (0.8.4); forest default again (after 0.8.4) | yes | forest vs GBM equal on the mocks (0.815 / 0.816) |
| training | held-out species and clades, scenarios, in-silico and real strains, depth points | yes | defines the reported F1 (E1-E4) |
| training | own Illumina model instead of ART (0.7.9); HiFi Q by length (0.7.3); Ultima (0.8.4) | yes | realism gaps (C4) |

## 3. Findings

Ranked by the damage they can do on real samples, times how likely that is. **Verified** means checked against the
code for this report; **read once** means one reviewer's reading.

### E. What is evaluated is not what users run

**E1. The models are never scored on the database that ships.** High; verified.
- **Evidence:**
  - `build_gtdb_database.py` (`dd01eb1`) ~2406-2411: the finished database "is never profiled for training or
    evaluation".
  - The training and test collections both run on `training_db` (`collect_command`, ~4154-4169).
- **How the databases differ:**
  - **Held-out species:** the training database lacks the 30% held-out species and the held-out clades
    (`--holdout 0.3`). The shipped database has them.
  - **Strain alleles:** the training database takes them from half of each species' genomes
    (`--allele-genome-share 0.5`), the shipped one from all of them (`--final-allele-genome-share 1`). At v21 that
    meant 29,019 species with alleles against 33,279.
  - **Unique k-mers:** a k-mer a kept species shares only with a held-out congener is unique in the training database
    and not in the shipped one. `gtdb_to_protal_db.py --exclude_species` drops held-out species from the full
    reference too.
- **Mechanism:** about 30 of the 98 default features describe the database in use:
  - `ref`, `neighbourhood`, `gaps`, `alleles`, `polymorphic`, `weights`, `strain_identity` and `ancestry`;
  - the unique k-mer rates;
  - `unexpected_congener_fit_share`.

  `model_features.py` says so for `ref` and `neighbourhood`. Where the shipped database differs, a sample's features
  shift relative to everything the model was trained and tested on:
  - species whose congeners were held out get fewer unique k-mers and lower MAPQ;
  - with twice the alleles, more reads become unsure, and the allele shift (B1) settles them;
  - the raw unique k-mer rates fall. Only the enrichments normalise for that.
- **Direction and size:** unknown. Nothing reads the finished database before it ships: no parity check, no test.
- **Check:** profile the test collection once on the finished database. Score only rows of species that were not held
  out, and compare with the same rows on the training database. That is one protal run per read type.

**E2. The "independent" test set shares the training's held-out species and genome pool.** Medium-high; verified.
- **Evidence:** both collections get the same `training_db`, genome table, held-out list (`--novel_species heldout`),
  clouds and in-silico strains (`collect_command`, ~4154-4169). Only the seed differs (`args.seed + 1000`, ~4269).
- **Mechanism:**
  - A real sample's novel species are other species. The test set's are the same ones training saw.
  - The `ref` rates are constant per species, and the trainer itself warns that they "can name a species". A model can
    therefore learn which database species turn into false positives beside these particular held-out congeners, and
    the test set rewards it.
  - This is the class of leak that hid the `foreign` features' leak in v17.
- **Partly guarded:** the "by species" CV holds a taxon's rows out of its fold, so the CV numbers are less exposed
  than the test set's. The reports quote both.
- **Also:**
  - The test seed is fixed (1001 with the default `--seed 1`), so builds v12-v28 reuse the same test communities
    whenever the genome table and design are unchanged.
  - Defaults were chosen on test-set differences within ±0.001-0.003: the scenario weight, gradient boosting, the
    distance features, `AUTO_MIN_GAIN`.
  - The global knob is vetoed by test F1 (`machine_learning_cmdline.py` ~1232-1237: a knob that scores lower on the
    test set than the default is dropped), and the same test F1 is then reported.
- **Check:** score a test set drawn from a different held-out list (a second split of the same release), or report
  test F1 separately for absent species whose held-out congener never appeared in training.

**E3. The headline F1 counts only taxa that have a row.** Medium; read once, consistent with the
[v25 report](2026-10-12-r226-v25/README.md).
- **Evidence:**
  - The trainer scores the rows of the test table, micro-pooled (`metrics()`).
  - A present species without any record has no row, so it is neither TP nor FN.
  - The marker-truth F1 is the only number that counts such species, and only if their marker genes got reads. It
    broke twice: v25 printed 0.633, and the build audit's G1 and G2 lose it.
- **Mechanism:** a change that strips a present taxon of all its records improves the headline. At v24b, 37 of the
  180 pe test misses had no row.

  Users count every truth species. The mocks give 0.816 that way and 0.927 over species with a record.
- **Also:** micro-pooling weights the deepest samples most, and the design's deepest test point is 5M pairs
  (`--test-read-pairs`). Users' samples are 5-50M pairs, and only the scenario hold-outs (4 samples each) reach them.
- **Check:** make the F1 over all truth species with marker coverage the headline, and keep it on every path a rerun
  takes (G1, G2).

**E4. "By species" cross-validation keeps the sample in training.** Medium; verified.
- **Evidence:** `folds()` groups the `species` scheme by `taxon` alone (`GroupKFold` on `df["taxon"]`). The trainer's
  docstring says each row is scored by forests that saw neither its sample nor its species; that holds only for the
  two schemes taken separately.
- **Mechanism:** six sample-level features take one value per sample and are not rounded:
  - `sample_log_fragments`, `sample_log_taxa`, `sample_low_identity`, `sample_identity`;
  - `sample_read_length` and `sample_read_length_sd` (`Profiler.h:4033-4051`).

  A tree can learn a per-sample offset from the sample's other species. The forest's `min_samples_leaf` is 1
  (`MODEL_DEFAULTS`). The knob, the depth knobs and `--features auto` all take this estimate.
- **Measured** ([v13](2026-10-06-r226-v13-soil/README.md)):
  - The v13 report found the mechanism: "cross-validation by species does not see it (every sample is in every
    fold)".
  - Rounding `sample_log_fragments` to 0.25 raised the boosted pe model's shallow-soil hold-out F1 from 0.881 to
    0.931.
  - The forest exploits it less: its by-sample estimate for shallow soil was 0.925, against 0.930 by species.
  - Only the depth was rounded. The other five sample-level features, and `sample_read_length` (new since then),
    are exact values.
- **Check:** group the folds by sample and species together, or report the by-sample estimate beside it when choosing
  features and the knob.

### C. Train/run agreement

**C1. A model is tied to the code that computed its features by their names only.** Medium; verified.
- **What is checked:**
  - `ModelContractProblemInXml` (`Profiler.h:3291-3330`) checks that the names exist, that a TRUE label exists and
    that a taxon can be scored.
  - `build_settings.tsv` records only the traceback, the fragment count and the gates' C, score, hits and pick
    (`BuildSettings.h:154-170`). Runs follow these and warn on a mismatch (`Run.cpp` `ApplyBuildSettings`).
  - The parity check re-profiles saved SAMs at build time, so it cannot see a change in the alignment-time tags (ZA,
    ZF, ZU/ZT/ZS, ZC, MAPQ).
- **What changes alignments without being recorded:** the weighted site shift, allele scores, `--align_top`,
  `--adaptive_candidates`, mate guidance, `max_score_ani` (0.9, ONT 0.85), `--keep_foreign_genes`,
  `--keep_suspect_copies` and `--no_unique_cores`. Each is silently whatever the binary's default is.
- **Already happened once:**
  - A database built with 0.8.2 (r226 v24b) was trained with the weighted site shift on.
  - Since 0.8.3 the default is off, so a default run on that database aligns differently from its training, without
    a warning.
  - The shift lowered short-read species depths by 1-2% ([v24b](2026-10-12-r226-v24b/README.md)).
- **Also:**
  - `--profile_only` takes the model's gate from the run's options, not from the SAM (`Options.h:1815`), and its
    gate warning is suppressed (`Run.cpp:202`).
  - A `--no_unique_cores` run on a database with the table only logs a line (`Load.cpp:339-361`). Its ZS/ZU/ZT then
    mean the anchor's seeds, not the covered cores the models were trained on.
- **Fix:**
  - Write every alignment-side setting the features depend on into `build_settings.tsv`, follow it with `auto`, and
    warn on a mismatch as for the traceback.
  - Stamp a feature-code version into each model's PMML (or `build_settings.tsv`) and warn when the binary's differs.
  - Write the gate, the pick and the unique-k-mer scheme into the SAM header for `--profile_only`.

**C2. Read-type detection can send ONT reads through the PacBio model, silently.** High for the samples it hits;
verified.
- **Evidence** (`ReadTypeDetection.h:131-151`):
  - A file is long-read if more than a quarter of its first 200 reads exceed 1 kb.
  - ONT is recognised only by a UUID read id or ` runid=` in the header (`IsOntName`). Otherwise a median read quality
    of Q25 or more means PacBio, and FASTA or Q0 throughout means PacBio.
  - The decision is printed as a `Note:`, not a warning (`Options.h:2146`).
- **What it catches:**
  - **ONT reads from SRA/ENA:** their ids become `SRR….N`. Unless the original header with `runid=` is kept, they are
    judged by quality.
  - **R10.4.1 sup and duplex reads:** median Q20-30. Duplex ids are `uuid;uuid`, which is no UUID.
  - **SRA Lite:** every base Q30, so PacBio.
  - **ONT FASTA:** always PacBio.
- **What a wrong type costs:**
  - the identity floor (0.9 instead of 0.85; the ONT proxy identity counts indels twice);
  - the k-mer screen instead of the indel bound;
  - the PacBio gate (16 instead of 8);
  - the HiFi model;
  - FASTA reads filled with Q30 instead of Q18 (`ReadType.h:33-35`).

  Expect many misses: ONT reads fail the floor or look like relatives' reads.
- **Also:**
  - A nanopore file whose reads are nearly all ≤1 kb (cfDNA, aDNA) is aligned and scored as single-end Illumina,
    with no note.
  - The training always passes `READ_TYPE` (`collect_training_data.py`), so no benchmark ever tests the detection.
- **Fix:**
  - Look for a UUID anywhere in the header (SRA keeps the original name after the id in most dumps), and accept
    `uuid;uuid`.
  - Warn, not note, when a long-read type was decided by quality or by the absence of qualities.
  - Say in the docs and on the website that `--read_type` matters for ONT data downloaded from SRA/ENA.

**C3. Placeholder base qualities break the long-read models' main feature.** Medium-high for those inputs; verified.
- **Evidence:**
  - `MeanErrorProbability` (`Profiler.h:609-619`) gives `!` (Q0) an error probability of 1.0.
  - `FillMissingQuality` fills only empty quality strings (`SeqReader.h:91-93`). Detection, by contrast, treats
    all-`!` as "no qualities" (`ReadTypeDetection.h:111`).
- **Mechanism:**
  - A FASTQ with `!` throughout (pbsim3, some BAM-to-FASTQ exports without QVs) gets `ReadExcess` ≈ −0.99 for every
    read.
  - In the r226 v23 PacBio model (GBM), `excess_scaled_median` carries 70.7% of the split gain
    ([v23](2026-10-11-r226-v23/varimp_groups.txt), which lists pe and pb only). The ONT model's share was not checked
    here.
  - The sample's every taxon sits far outside the training range, and a tree extrapolates as a constant: the calls are
    arbitrary. No warning.
- **Same mechanism, smaller:**
  - SRA Lite's constant Q30 on ONT reads makes their own differences look like divergence (excess +1-3%, beyond
    `kHighExcess` 0.02).
  - Phred+64 files are read as Phred+33; no encoding check exists.
- **Fix:** treat a quality string of all `!` like a missing one (fill with the read type's quality, as for FASTA) in
  the reader. Warn when a sample's median quality lies outside its read type's training range.

**C4. Inputs outside the training ranges.** Medium; the ranges were verified, the effects are plausible.
- **Single-end short reads:**
  - The se model trains on 150 and 250 bp reads and Ultima only (`--se-read-setups`, ~3090). pe trains on 75-250 bp.
  - Any single short-read file is se. A 150 bp model gave 75 bp reads 2-6 times the false positives
    ([report](2026-10-14-length-free-features/README.md)), and se's model has no 75 or 100 bp rows to learn from.
- **Depth:**
  - The pe design ends at 30M pairs with one sample per setup, the scenarios near 40M. Long reads end at 6 Gb (design)
    and about 12 Gb (scenarios).
  - Deeper real samples fall into the deepest leaf, while the evidence of absent taxa keeps growing (about
    depth^0.84).
- **Base qualities:**
  - The simulator draws each Illumina error at exactly its quality's probability (calibration 1.0; MiSeq 1.5), and
    no default setup is NovaSeq (4-bin qualities, poly-G).
  - The `excess_*` features measure divergence beyond the qualities, so real miscalibration shifts them.
- **ONT:** trained on one error model at 0.97 accuracy (`--ont-setup …:0.97`). R10.4.1 sup reads at Q20+ are more
  accurate, so a relative's reads at ~97% look more like a present species' own.
- **What the simulator does not make:**
  - adapters (a fragment is at least a read long);
  - PCR duplicates, GC or coverage bias (fragment starts are uniform), chimeras and contamination;
  - in-silico strains with indels or differences in gene content.

  These shift the coverage features (`depth_cv`, `breadth_ratio`, `gene_presence_ratio`) and the identities.
- **Read length:** real reads trimmed to many lengths have a `sample_read_length_sd` the fixed-length training reads
  never show.
- **Fix:**
  - Add 75 and 100 bp to `--se-read-setups`, and a NovaSeq setup with miscalibrated qualities to the design.
  - Log per sample which features fall outside the model's training range (the trainer knows the ranges; store them
    in the PMML or beside it).

### B. Built-in biases

**B1. The strain-allele shift only ever helps, so it favours species with many genomes.** Medium; verified.
- **Evidence:**
  - `SiteShift::Half()` = `2*best - 2*known - variable - discounted` is zero or negative (`StrainAlleles.h:786`), and a
    copy without alleles gets 0 (`ShiftOf`).
  - The shift applies to every candidate within the unsure margin (`SettleBySites`, `AlignmentStrategy.h:962`) and
    is folded into the score that ranks the records and sets MAPQ (`AlignmentOutputHandler.h:443, 500`).
- **Mechanism:**
  - An unsure read of a one-genome species B, within 3 mismatches of a congener A with alleles, can move to A. That is
    an FN for B and an FP for A.
  - Even when A was already best, the wider gap raises MAPQ past the profiler's `m_min_mapq` 4, so the read now
    counts in A's `fragments`.
  - The more genomes a species has (up to 16 alleles), the stronger the pull.
- **Measured at v21** (scores and features together):
  - real strains of species without alleles were missed 10-15% more (pe 4.7% → 5.2%);
  - false positives of species with alleles beside a held-out congener rose (pe 2.58% → 2.79%, ont 3.56% → 3.96%).
- **Never measured:** the alignment shift on its own at r226. The v21 report's recommendation 2, "profile the test set
  once with `--no_allele_scores`", is still open.
- **With E1:** the shipped database's alleles come from twice the genomes the models were trained with, so this bias
  is stronger in user runs than in any evaluation.
- **Check:** the `--no_allele_scores` run, split by species with and without alleles, on the training and on the
  finished database.

**B2. The sparse-first gates drop short and divergent reads, with all their evidence.** Medium-low; verified.
- **Evidence:**
  - **The keep rule:** a read (a pair counts both mates) is kept if one of its sampled lookups hits: every C-th
    syncmer's core in the index, with a flex cell sharing ≥14 of 16 flank bases (`ChainAnchorFinder.h:1103-1112`,
    `:337-340`).
  - **What a dropped read leaves:** no record and an empty failed-candidate list (`Classify.h:751-775`).
  - **Ungated seeding** has no absolute flank threshold.
- **Estimate** (binomial, independent errors; overlapping k-mers make it worse; not measured). Dropped at C 4:

  | reads | 3% from the reference | 5% | 10% |
  |---|---|---|---|
  | se 75 bp (~3 picks) | 5% | 17% | 59% |
  | se 100 bp | 2% | 10% | 49% |
  | pe 2×75 | 0.3% | 3% | 35% |
  | pe 2×150 (~14 picks) | ~0 | 0.03% | 8% |

  ONT reads under ~500 bp also lose picks to flank indels without `--flex_indels`.
- **Mechanism:**
  - **FN:** thin taxa lose reads, and 60-76% of r226's misses had ≤2 reads of their own.
  - **FP:** the most divergent reads, those that would fail on a congener and count against it in
    `failed_candidate_rate`, are dropped most.
  - Training and runs gate alike, so the models are calibrated for the simulated mix of divergence and read length,
    not for a real one.
- **Measured:** no F1 loss up to pe/se C 4, pb C 16, ont C 8 on the r226 v27 test set (150 bp reads mostly), and the
  same calls on the mocks.
- **Untested:** se 75-100 bp reads, and samples dominated by species 5% or more from their references.

### Lower-severity findings

Read once unless marked.

| # | finding | evidence | effect |
|---|---|---|---|
| L1 | Mate guidance completes a pair on the guide's taxon and ignores the other mate's own best taxon; split-gene pairs get the consensus MAPQ (0.6.0a gave ~0) | `MateGuidance.h:144-146`, `ReadConsensus.h:88-100` | FP for the guide's taxon when it is wrong; no double count |
| L2 | Ranking and MAPQ use Bitscore (+2/3/1/2), WFA2 optimises 4/6/2: equal-penalty alignments differ in MAPQ; a read overhanging its gene's end loses ~2 per clipped base against a congener whose gene covers it | `AlignmentUtils.h:167, 713` | gene-end FN/FP; the anchored alignment changed MAPQ for 0.14% of pairs ([round 1](2026-09-29-performance-profiling/README.md)) |
| L3 | Adaptive candidates (7) are taken in anchor order (taxid-first) and only in the best alignment's genus; a split genus (`Genus_A`) is never tried | `ChainAnchorFinder.h:1291`, `KmerLookup.h:131` | strains of split genera missed |
| L4 | Exact ties go to the lower internal taxid (primary record, ZA, `fragments_all`); MAPQ 0, so `fragments` unaffected | `Classify.h:858` | small, deterministic |
| L5 | Long reads: `fragments` counts gene segments, not reads (docs say "a long read once"); `--fragment_count soft` is near a no-op for pb/ont | `Profiler.h:1363`, `features.md:101` | docs; training and runs agree |
| L6 | `unique_cores.bin` fit test accepts up to 3 missing cores (`positions.size() + 4 > capacity`); a mismatch falls back to the anchor's counts per gene, so one sample can mix schemes; `Mismatched()` is never reported | `UniqueCores.h:262-266` (verified) | only with a table of another database |
| L7 | On an index from before `78f7a92`, single-value cores count as unique without the check | `KmerLookup.h:385-392`, `Seedmap.h:1011` | FP on old databases; no warning |
| L8 | Reduced (`--genes`) databases keep k-mers unique that other taxa's left-out genes share | `Build.h:2411-2413` | FP; opt-in |
| L9 | The AVX-512 `--flex_indels` kernel has never run on AVX-512 hardware; `--test-c-gates` flex rows run on Zen 4 | `FlexScan.h:246-314` | opt-in; `PROTAL_SIMD=avx2` until verified |
| L10 | N in a read is looked up as A/T (the index has no N k-mers) | `KmerIterator.h:217-226` | pre-existing; candidates only |
| L11 | Mate files in different orders but equal counts only warn on names | `Align.cpp:573-577` | user error, mispaired |
| L12 | A build with `-march=native` turns on WFA2's AVX2 extend kernels, which change end-free ties | `lib/wfa2-lib.cmake` | MAPQ and allele shifts differ from the database's build |
| L13 | `--depth_identity_margin`'s help says such reads "still count for detection", but `depth` and `low_identity_share` are features | `Options.cpp:60` | opt-in; misleading help |

## 4. Checked and found safe

- **Input:**
  - Gzip read in parallel pieces: a piece counts only where the one before ends, and CRC and length are checked.
  - Many-member gzip: a wrong cut falls back to one thread.
  - BGZF and zstd frames: block checks; the file must end on a frame or block.
  - Lanes as one sample: a batch ends at every file end.
  - Mates are read in step, and different read counts fail the sample.
  - Malformed FASTQ fails, CRLF is handled, lowercase is uppercased (0.6.0a seeded it as A).
  - No read is filtered before seeding.
  - The sample-level features count aligned fragments, not reads read.
- **Seeding:**
  - The scalar, AVX2 and AVX-512 syncmer and flex kernels give the same codes, ties and tails; tested at every length.
  - The compact and direct key maps give identical starts and slots; failure is loud.
  - The grid sparse load keeps the flags and flex cells.
  - The uniqueness pass, orientation and N handling agree between index and lookup.
  - No floating point in seeding; prefetches change nothing.
  - At C 1 gated seeding equals ungated.
  - A database without `build_settings.tsv` runs ungated; one with gates but no pick uses hash.
- **Alignment:**
  - The k-mer screen, the indel bound and the anchored flanks are exact, with the right per-operation bounds.
  - The halved penalties and X-drop change nothing; 78,777 replayed WFA2 calls were identical (`aa36fc3`).
  - Rounding never makes the budget stricter.
  - The Singletrack fallback is per alignment.
  - Caches and copies make no thread dependence.
  - Pairing across genes needs the same taxon and an expected neighbour.
- **Profiling:**
  - Features reach the model by name, so their order does not matter; a missing one stops the run (exit 2).
  - Training dumps and scoring share one feature writer.
  - No NaN or infinite path in the default getters.
  - Overlapping mates are counted once.
  - ExactSum and sorted aggregation give thread and order independence.
  - Variants are computed with `--no_strains` too, so the collections see the same features.
  - The single-seed failed-candidate change is inert under the stride pick.
  - The read EM and the depth margin feed features only. They reassign no `fragments` and decide no call.
  - The 100-sweep cap is the same in training.
  - The suspect-copy scan flags 0.04% of r226's copies ([v7/v8](2026-10-03-r226-v7-v8-evaluation/README.md)).
- **Training:**
  - The truth only labels; no read name or quality reveals a read's source.
  - `meta_*` columns are never features.
  - The allele/strain split uses the same FNV-1a hash in Python and C++.
  - Held-out species stay out of the unseen tables, the column weights and the training database's neighbours.
  - Training and test communities use different seeds.
  - The earlier leaks (priors, foreign, depth identity) are out of the default set.

## 5. What to check next

In order of what they tell per unit of work. None needs new code except 3, 4, 6 and 7.

1. **Profile the test collection on the finished database** (E1, B1): one `--profile_only` is not enough, because the
   alignment differs, so one protal run per read type. Score the rows of species that were not held out, against
   the same rows on the training database. If F1 moves by more than the noise (±0.002), the build should evaluate
   the finished database by default.
2. **Run the test set once with `--no_allele_scores`** (B1), on both databases, split by species with and without
   alleles. This is v21's open recommendation.
3. **Record and check the alignment-side settings** (C1): the weighted site shift, allele scores, `align_top`,
   adaptive candidates, mate guidance, the identity floors and the unique-k-mer scheme, in `build_settings.tsv`;
   a feature-code version in each model; the gate and scheme in the SAM header.
4. **Fix the read-type and quality traps** (C2, C3): a UUID anywhere in the header, `uuid;uuid`, all-`!` as missing,
   and a warning when the type rests on qualities. These are small, and they are what a user downloading public ONT
   data meets first.
5. **Make the headline honest** (E2-E4):
   - F1 over all truth species with marker coverage as the headline (E3);
   - folds grouped by sample and species for the knob and the feature choice (E4);
   - the knob veto off the reported test set, and the test seed rotated between builds (E2).
6. **Widen the training ranges** (C4, B2):
   - se 75 and 100 bp;
   - a NovaSeq setup with miscalibrated qualities;
   - one deeper point per read type;
   - in-silico strains with indels.
7. **Log when a sample's features fall outside the training range** (C3, C4): the cheapest guard against everything
   in C, and it would have flagged the mocks' near-misses if that is what they are.
8. **Use the mock benchmark as the regression gate between versions** at r226 (`scripts/mock_benchmark`). It is the
   only real-data F1 there is. The near-misses with many reads (*Akkermansia*, the Dohlman *E. coli*, *B. ovatus*) are
   worth a feature-by-feature look against the training rows of present species at the same depth: that is where E1,
   C3 and C4 would show.

**Website.** This report changes no code. The website's options page should say that `--read_type` matters for ONT
reads downloaded from SRA/ENA (C2); `docs/running.md` should say it too.

## 6. Follow-up, 2026-10-15: what was done

The user asked for every check of section 5 to be implemented (1 inside `build_gtdb_database.py` as a final check, 2
with `--share-logs`), for L1 to become a vote of both mates, and for the *Akkermansia* near-miss to be explained. The
work went into side branches, each tested on its own, then merged into `f1-integration` (from `audit-fixes` at
`4d89d86`) and tested together. Nothing has run at r226 scale: the next r226 build is the first measurement of what
the build-side changes do to F1.

| Finding, check | What changed | Commit (branch) |
|---|---|---|
| E1, check 1 | the build profiles the test samples again on the finished database and scores the species not held out beside the training database's rows (`model_logs/final_db_check.tsv`, a warning beyond 0.002 of F1); on by default, `--no-final-check` | `dcd8ae1` (`f1-build-pipeline`) |
| B1, check 2 | with `--share-logs`, the short-read test samples profiled with `--no_allele_scores` too, split by species with and without strain alleles (`model_logs/no_allele_scores.tsv`, in the share archive) | `dcd8ae1` |
| C1, check 3 | `build_settings.tsv` records the run settings (allele scores, site shift, `align_top`, adaptive candidates, mate guidance, identity floors, unique k-mers, the three profiling ones); runs follow them (`auto`) and warn when told otherwise; every SAM header carries the settings it was aligned with, which `--profile_only` checks; a feature-version table (`protal --feature_versions`) that the trainer writes into each model and protal checks at load | `b39ce7f` (guards), `f3eb3f8` (trainer), `dcd8ae1` (build passes it) |
| C2, C3, check 4 | ONT reads from SRA/ENA (UUIDs anywhere in the header, `uuid;uuid`) and placeholder qualities (all `!`) | `9b7357d` (`claude/exciting-swartz-159a28`) |
| E3, check 5 | the headline F1 counts every truth species the marker genes can see (`metrics.json` `headline`, first in the report and in `summary.txt`) | `f3eb3f8`, `dcd8ae1` |
| E4, check 5 | the knob, the knob curve, the calibrated calls and `--features auto` chosen on folds that hold out samples and species together (`--cv-scheme samples+species`) | `f3eb3f8` |
| E2, check 5 | no knob veto by the test set; the test set's seed drawn per build (kept for reruns, `--test-seed`) | `f3eb3f8`, `dcd8ae1` |
| C4, check 6 | single-end 75 and 100 bp, a NovaSeq setup whose qualities overstate accuracy by 3 Phred, one deeper sample per read type, codon indels in the in-silico strains (+~11% of the collections' alignment, estimated) | `dcd8ae1` |
| C3, C4, check 7 | every model holds its features' training ranges; protal warns per sample when taxa near the knob fall outside them and names the features in `.profile.log` (`outside_training_range`) | `f3eb3f8`, `b39ce7f` |
| check 8 | the Dohlman D6331 near-miss traced to PCR piles on single genes; spike mode (`--spike_genes auto`) leaves such genes out when the sample's coverage shows the pattern ([report](2026-10-15-mock-near-misses/README.md)) | `21c5c38` (`f1-spike-genes`) |
| L1 | `--mate_guidance consensus` (both mates vote; `guide` the old behaviour), recorded in `build_settings.tsv`; a database that records none keeps `guide`, which its models were trained with ([report](2026-10-15-mate-consensus/README.md)) | `651d995`, merge `53e25ab` |
| E3 (releases) | `build_gtdb_releases.py` took a scenario's or the marker truth's last summary row as "held out"; now each from its own row, with the marker-truth F1 beside them | `1094c6d` |
| (found on the way) | a map naming a profile folder not yet there lost the first sample's `.truth_annotated` (and would fail spike mode's file) | `b61e16e` |

**Measured so far** (reduced r226 database, WSL, models of `5369c99`):
- Spike mode: only the contaminated D6331 library goes into it (10/0/2 → 12/0/0); no other mock or simulated
  sample's probabilities move.
- Mate consensus with models trained on `guide`: calls within 1-2 species on the simulated sets, mocks identical, test
  set F1 0.9597 → 0.9618; +2.65% of the paired-end alignment loop's instructions.
- Trainer, c4test pe: knob and F1 within ±0.002. The GBM loses 0.003 of held-out F1 when the sample is held out
  too, the forest 0.0003, as E4 expected ([report](2026-10-15-trainer-honest-f1.md)).

**Still open:** the mock benchmark as a regression gate (check 8's first half); every r226-scale number above. B2 and
L2-L13 were taken up afterwards ([report](2026-10-15-f1-audit-lower-findings/README.md)): B2 measured and fixed
(`--sparse_first_min_picks`), L2 an option (`--gene_end_clips`), L6, L7, L10-L13 fixed, L3 counted, L4, L5, L8, L9
documented.
