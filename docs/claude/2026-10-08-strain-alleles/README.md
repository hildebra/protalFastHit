# Strain alleles: the species' other genomes as diffs of the reference genes, in the scores and as features

**Question (the user's, 2026-10-08).** After the v18 report (section 5: no build had alleles; they would reach 60-80% of
the false negatives), the user asked how to store extra strain alleles. Then: "implement the diff code in C++ (steps
2-4), I will run this with the next DB build; agreed on circularity and leaks; can this be coupled to the same genomes
used in the ancestry checks?"

**Short answer.** Implemented on branch `strain-alleles` (from `audit-fixes` `7b44001`), not yet built at GTDB scale:
- `--build` writes `strain_alleles.tsv`: per species' copy of each marker gene, up to 4 alleles of the species' other
  genomes, as edits of the representative's copy.
- A run scores each short read's candidate alignments with their species' best allele, so a read of a known strain
  wins on its species. It also counts two `alleles` features, which are in the default set (82 features).
- The GTDB build takes the alleles from half of the genomes, chosen by a hash of the accession, and simulates strains
  only from the other half. No simulated strain is its own species' allele, so the circularity the v18 report found in
  the fixed sites cannot arise here.
- The coupling to the ancestry checks: the ancestry report takes its alleles from the same hash half
  (`ancestry_sites.py --allele-genome-share`). Its polymorphic sites are those of genomes protal can know, never a
  simulated strain's own, with the read's own genome left out on top (`1c184fb`).

## Design

**What is stored** (`src/SequenceUtils/StrainAlleles.h`).
- An `Edit` of the representative's copy is a substitution (position and base), an insertion of n bases before a
  position, or a deletion of n bases from one, in 4 bytes.
- An `Allele` is the stretch `[begin, end)` of the representative it covers, plus its edits there.
- `Table` holds per (taxid, gene) up to 4 alleles. It answers `Of(taxid, gene)` with a span and `Get` with a view of the
  edits.
- `strain_alleles.tsv` has one line per copy: `taxid<TAB>gene<TAB>begin-end:edit,edit;begin-end:...`. An edit reads
  `123G`, `300i3` or `450d6`.
- The table is packed into `database.protal` like `congener_gaps.tsv` (`BundleSources`). `--add_tables` accepts it, and
  the `--no_bundle`, `--unpack_db` and `--add_tables` help texts name it.

**How `--build` makes it** (`src/SequenceUtils/StrainAllelesBuild.h`, `Build.h WriteStrainAlleles`). The pass runs
right after the congener gaps, as one pass over the full reference on all threads, framed as for the conservation
factors.
- **Records it reads:** every record of a database species and gene from a genome that may give alleles
  (`AlleleGenome`: the FNV-1a 64 hash of the accession, the record's second header word since `1c184fb`, top 53 bits as
  a fraction below `--allele_genome_share`).
- **Records it skips:** copies identical to the representative's, and sequences seen before for that copy (by hash, up
  to 64).
- **The diff:** the rest are aligned against the representative's copy with WFA2 (4/6/2, both ends 15% free). Each end
  is trimmed to its first run of 10 matches, because ends-free alignment can skip a prefix of only one copy. A copy that
  starts elsewhere and has a flank of its own would otherwise give its flank as edits; the first unit-test run caught
  this.
- **What is rejected:** an allele covering less than half of the representative, with more than 0.1 of its columns
  differing, or as far from the representative as the copy's nearest congener's copy (`congener_gaps.tsv`). The last
  would be a misassigned genome or a congener's gene, whose reads the species would then take. A congener with an
  identical copy admits no allele.
- **The sample:** per copy, the 16 distinct alleles with the least hashes. The sample does not depend on the order the
  threads read the frames in.
- **The choice:** up to `--strain_alleles` (4) of the sample, farthest first from the representative and from the
  alleles chosen before (k-center).
- The log line counts every filter: `Strain alleles: N alleles (E edits) of C gene copies of S species, up to 4 each (of
  R full-reference copies...)`.

**How a run uses it.** The table is loaded whenever reads are aligned or profiled (`RunProtal LoadStrainAlleles`), not
only for profiling, because the scores need it.
- **The scores (short reads; `--no_allele_scores` off).**
  - After a candidate alignment is made and checked, `SimpleAlignmentHandler::ScoreAlleles` takes the read's differences
    from the alignment's columns (`FromColumns`) and finds the copy's best allele (`BestAllele`).
  - The best allele is the one with the least shift: the read's differences it does not explain, plus its own
    differences inside the read's span that the read lacks, minus the read's current differences. A substitution
    matches by position and base; an indel by kind and length within 4 bases.
  - `AlignmentInfo::allele_shift` (≤ 0) enters `Score()`, which counts the explained differences as matches. So the
    candidates' ranking, the mates' pairing, MAPQ and `AS` all see it, and the CIGAR and the identity stay the
    reference's.
  - The log's `strain alleles:` line counts the candidates on copies with alleles and those an allele scored higher.
  - Long reads are not changed (their chains score elsewhere); they get the features.
- **The features (every read type).**
  - `RecordEvidenceCollector::NoteAlleles` takes each kept record's differences from its CIGAR and SEQ (`FromSamRecord`)
    and finds the best allele the same way.
  - `RecordEvidence` sums records, differences, explained differences, the edits gained and the aligned columns, all
    integers and so independent of threads.
  - `allele_explained_share`: of every kept record's differences, the share their copies' best alleles explain.
  - `allele_identity_gain`: the edits gained over every kept record's aligned columns.
  - Both are -1 without the table, and 0 where nothing is explained, a species without alleles too.
  - Group `alleles` in `scripts/model_features.py`, in the default set.
  - `allele_copy_share` (the share of kept records on copies with alleles) is a column of the training table in no
    feature group: see the next paragraph.

**Availability is not a feature.** Whether a species has alleles at all is whether GTDB has other genomes of it: the
cluster size that the `priors` group carries, opt-in because the simulation cannot test it (the r226 v9 evaluation:
its whole gain was the rule that a divergent read cloud on a one-genome species is a relative the database lacks). The
first version of the features had -1 for a species without alleles and a share of records on copies with alleles in
the group: both would have told the models that rule again. Now the two features give such a species 0, the value of
reads its alleles would not explain, and `allele_copy_share` stays out of every group. The scores still move reads to
species with known strains, which is the point in use as in training.

**How training stays honest** (`scripts/build_gtdb_database.py`).
- `--allele-genome-share` (0.5) and `--strain-alleles` (4) go to both `protal --build` runs (`ALLELE_ARGS`).
- The genome table loses every non-representative genome of the allele half (`split_allele_genomes`, right after the
  representatives are read, before the in-silico strains): those genomes give alleles and are never simulated.
- A species left with one genome gets an in-silico strain, which has no alleles, as a species known from one genome has
  none in use.
- The finished database takes its alleles from the same half, so the features mean in use what they meant in training.
  An index of every genome's alleles in the finished database alone would shift them between training and use.
- The stage keys of both databases include the allele options.
- The ancestry report gets `--allele-genome-share` when the build has alleles. Its alleles are then the genomes' that
  give protal's alleles. `1c184fb` already leaves the read's own source genome out.
- **What it costs:** half of the downloaded strain genomes are no longer simulated, so more species are simulated
  through in-silico strains. The ratio is in the console's genome-table line and `genome_table.txt`.

**What the user asked about coupling.** The coupling is the hash rule, not a list file:
- `strain_alleles::AlleleGenome` (C++) and `gtdb_to_protal_db.allele_genome` (Python) are the same FNV-1a rule.
- Both tests pin the same values: `GCF_000005845.2` 0.278, `GCA_000001405.29` 0.090, `GCA_900000000.1` 0.696,
  `GCF_000195955.2` 0.508.
- protal's alleles, the genome table's split and the ancestry report's alleles therefore agree without a file passed
  around.
- The report keeps up to 6 alleles per copy where protal keeps 4 chosen farthest first. Its polymorphic sites are a
  superset of protal's alleles' edits, from the same genomes.

## Files

- New: `src/SequenceUtils/StrainAlleles.h` (table, edits, read differences, best allele, the hash rule),
  `src/SequenceUtils/StrainAllelesBuild.h` (diff, sample, choice, the per-record filters),
  `tests/test_StrainAlleles.cpp` (8 tests).
- Changed:
  - `src/Options.h`: `--strain_alleles`, `--allele_genome_share`, `--no_allele_scores`, the table's file names, validation
    and help texts.
  - `src/Build.h`: `WriteStrainAlleles`; `WriteCongenerGaps` returns its table; bundling and `--add_tables`.
  - `src/SequenceUtils/GenomeLoader.h`: the table.
  - `src/RunProtal.h`: loading, the handler's setting, the sample's log line.
  - `src/Alignment/AlignmentUtils.h`: `allele_shift` in `Score()`.
  - `src/Core/AlignmentStrategy.h`: `ScoreAlleles`, its counts.
  - `src/Profiling/Profiler.h`: `NoteAlleles`, the evidence, the three accessors and features.
  - `scripts/model_features.py`: the `alleles` group, default set and named sets.
  - `scripts/build_gtdb_database.py`: the options, `split_allele_genomes`, `ALLELE_ARGS`, the keys, the ancestry
    report's share.
  - `scripts/mini_db/gtdb_to_protal_db.py`: `allele_genome`.
  - `scripts/ancestry_sites.py`: `--allele-genome-share`.
  - Tests: `scripts/test_model_pmml.py`, `scripts/test_ancestry_sites.py`, `scripts/mini_db/test_gtdb_build.py`,
    `scripts/mini_db/test_gtdb_pipeline.py`.
  - Docs: `docs/databases.md`, `docs/features.md`, `docs/running.md`, `docs/version_changes.md`.

## Tests (WSL, 4 cores, niced; `strain-alleles` at the commit)

- **C++ unit tests** (`tests/test_StrainAlleles.cpp`, 8 new):
  - the hash rule's values, the same as Python's;
  - the table's write, read and errors;
  - an allele's edits from its alignment: a substitution, an insertion, a deletion, both ends' overhangs, and an
    insertion before a substituted base (the table keeps the substitution first);
  - the sample and choice, independent of order and threads, and the build's filters;
  - what an allele explains of a read, from a SAM record and from an alignment's columns alike;
  - the score's shift;
  - the handler: a read of a known strain loses to a congener by the references alone and wins on its species with the
    alleles, its CIGAR unchanged;
  - the profiler's features with and without the table, on 1 and 3 threads.

  All 477 unit tests: 474 passed, 3 skipped.
- **Python:**
  - `test_model_pmml.py` pins the default set with `alleles`, `allele_copy_share` in no group, and a v18 table failing
    clearly;
  - `test_gtdb_build.py` pins the hash values and `split_allele_genomes`;
  - `test_ancestry_sites.py` checks the report's alleles restricted to the share.

  68 tests OK.
- **The CI suites** (`PROTAL_TESTS_REQUIRED=1`, before the last change below, otherwise the same code): the mini
  database, the end-to-end tests (146 OK), the mock community, the mini database and GTDB build tests (83 OK), and the
  script tests (81 OK).
  - The end-to-end database had alleles ("665 alleles (7277 edits) of 347 gene copies of 3 species"), so the end-to-end
    tests ran with the allele scores.
  - The mini GTDB build: "Strain alleles: 6953 alleles (131661 edits) of 3558 gene copies of 32 species, up to 4 each
    (of 19728 full-reference copies of the database's species and genes: 9301 of genomes outside
    --allele_genome_share 0.5, 3464 identical to the representative's, 7 repeated, 0 covering less than 0.5 of it or
    more than 0.1 apart, 2 as far as the nearest congener's copy or farther)".
- **Found by the tests:**
  - The first unit run: ends-free alignment aligned a copy's own flank as edits. Hence the trimming to runs of 10
    matches.
  - The first pipeline run: an insertion and a substitution at one position came out of the alignment in the order the
    table's check refuses. protal then refused the training database's table and the foreign scan's run failed. Hence
    the sort in `Diff` and its test.
- **The last change**, availability out of the features (above): afterwards all unit tests, the three script files and
  the pipeline's `test_a` passed again.

## Cost at GTDB scale (estimates, not measured)

- **The pass:**
  - It reads the full reference as the conservation pass does (40 s at r226 v18, 4.9 GB compressed for the training
    database).
  - Half the records stop at the hash. Of the rest, most are identical to the representative or repeated (big
    species).
  - Perhaps 10-20 million WFA2 alignments of near-identical copies of ~1 kb at 20-50 µs: 3-15 s at 84 threads.
  - The sample holds at most 16 alleles per copy, likely ~1 GB at r226.
- **The table:**
  - About 26,000 species have other genomes' copies that differ (the conservation pass: 25,929 at v18).
  - Some 3 million copies, at an assumed ~2 alleles of ~10 edits each.
  - About 400 MB of text, ~300 MB in memory at run time, a few seconds to parse beside the other tables.
- **The run:** per candidate on a copy with alleles, a pass over the read's columns and a lookup of the alleles' edits
  in its span: small beside the alignment.

## What to look at in the next build

- **Both index logs** (`index_and_package.log`, `training_db_index.log`): the `Strain alleles:` line, with how many
  copies each filter took (genomes outside the share, identical, repeated, too divergent, beyond the congener) and the
  allele count.
- **The console's genome-table line:** how many genomes went to the alleles. **`genome_table.txt`**: the share of
  real against in-silico strains after the split.
- **The collection's protal logs:** `Sample ... strain alleles: N candidate alignments on gene copies with known alleles,
  M of them scored higher`.
- **`allele_copy_share` in the training tables:** how strongly it alone separates present from absent rows. That is
  the cluster-size signal the features were kept from.
- **The trainer's varimp** of `allele_explained_share` and `allele_identity_gain`, and an ablation without the `alleles`
  group. The scores themselves cannot be ablated after the fact: they move reads. v18 against the next build is the
  comparison, with its other changes.
- **The ancestry report:** the fixed sites against the allele half's genomes. With the split and the leave-one-out,
  that is the honest number the v18 report asked for.

## What is not done

- No GTDB-scale build: the gain is unmeasured.
- Long reads get the features but not the scores.
- Indels in alleles are matched by kind, length and position within 4 bases, not by their bases.
- The table is text. At r226 a binary form would load faster, if its parse shows in a run's start.
