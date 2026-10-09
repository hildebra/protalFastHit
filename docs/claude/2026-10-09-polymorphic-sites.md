# The species' polymorphic sites: in the alignment of unsure reads, as features, and in the in-silico strains

**Data and commit.** Implemented on the branch `polymorphic-sites` from `audit-fixes` at `156ccda` (the strain alleles,
`504de80`). It was built and tested in WSL on 4 cores; the commands are below. It has not been built at r226: v19
(`a81cee3`) predates even the strain alleles ([2026-10-09-r226-v19](2026-10-09-r226-v19/README.md)).

**The request (the user's, 2026-10-09), after the v19 report's section 7:**
1. Use the ancestral and polymorphic sites to bias the alignments protal is unsure about.
2. Give the classifier the agreeing and disagreeing polymorphic sites seen per species, and the fixed-site agreement.
3. Make the in-silico strains realistic at the congener sites.

## What was built

### 1. Unsure reads take the species' polymorphic sites into their scores

[StrainAlleles.h](../../src/SequenceUtils/StrainAlleles.h) gains the polymorphic sites of a copy (`Polymorphism`).
These are the positions where any of its alleles differs from the representative, with the bases the alleles carry
there; an insertion or a deletion makes an indel site.

A candidate's site shift (`ShiftOf`, in half differences) has three parts:
- its best allele's shift, as in `504de80`: the read's differences that allele explains are matches, and the allele's
  own edits that the read lacks are differences;
- of the differences that allele leaves, those another allele of the copy has: matches too (strains recombine);
- those at a polymorphic site where no allele has the read's base: half a difference (the species varies there).

A difference where the species is fixed weighs in full.

The ancestry sites need no extra weight in the alignment. Where the species and its congener are both candidates, their
scores already differ by exactly those sites. The question the alignment lacked was how much a difference costs on
the species, which is what the polymorphic sites answer.

**Only unsure reads** (`SimpleAlignmentHandler::SettleBySites`):
- `AlignAnchor` computes the shift of each candidate on a copy with alleles, but keeps it pending
  (`AlignmentInfo::site_shift_pending`).
- Once all of a read's candidates are aligned, a read is unsure when a candidate of another species comes within
  `kUnsureMismatches` (3) mismatches of the best by the references alone. Then every candidate within that margin takes
  its shift (`site_shift`, which `Score` counts as `(match + mismatch) × site_shift / 2`).
- A read one species fits clearly better keeps the references' scores. So do its `AS`, its MAPQ and the identity
  filter.

In `504de80` every candidate took its best allele's shift, unsure or not. The log's `strain alleles:` line now also
counts the unsure reads and those whose best species the shifts changed. Short reads only, as before; long reads are
aligned through their chain elsewhere.

### 2. The `polymorphic` features (default set, 85 features)

In the profiler (`NoteAlleles`, kept records with SEQ on a copy with alleles):

| feature | what |
|---|---|
| `polymorphic_known_share` | of the polymorphic sites the records cover, the share where the read has a base (or an indel within 4 bases) one of the alleles has |
| `polymorphic_novel_share` | the share where it has another base or indel |
| `ancestry_fixed_gain` | the share of the species' base at the ancestry sites fixed within the species (no allele differs there), less that at all its ancestry sites, both weighted by the alleles covering each site (n/(n+1); 0 without one: unknown) |

- **Unknown and empty values:** all three are -1 without the table and 0 without a site, a species without alleles too
  (the `alleles` convention).
- **Diagnostics in no group:** `allele_sites_per_kb` and `ancestry_fixed_share` are columns of the training table. They
  say whether a species has alleles at all, which is the cluster size.
- **Where the gain comes from:** `ancestry_fixed_gain` is the v19 report's fixed-site signal in a form that is 0 without
  alleles. A real strain's congener-like bases fall on polymorphic sites, so its agreement rises without them. A novel
  congener's fall on fixed sites.

### 3. In-silico strains with congener sites

[insilico_strains.py](../../scripts/insilico_strains.py) `--clouds` gives each one-genome species its nearest congener
with a representative in the table. Each marker gene's congener sites are where that congener's copy differs (the
copies compared along their shared 12-mers, `aligned_bases`).

`--congener-share` of the gene's substitutions are then moved to those sites with the congener's base:
- the number is drawn per gene, binomially;
- each base placed sets back another of the gene's substitutions, so the gene keeps its divergence;
- a base that would make a stop codon is not placed.

`auto` (the default) measures the share on the same sample of real strains as the omega: their marker substitutions
that sit at a congener site and carry the congener's base, with at least 2,000 substitutions. Without clouds, or
without enough real strains to measure on, the share is 0 and the strains are made as before.

`insilico_strains.tsv` gains `congener` and `congener_substitutions`. `build_gtdb_database.py` now makes
`species_clouds.tsv` before the in-silico strains, whenever there are any or the hold-out is steered, and passes it to
them. The in-silico stage's key includes the clouds.

## Tests

- **Unit tests:** 476 passed, 3 skipped (was 474). The new ones in
  [test_StrainAlleles.cpp](../../tests/test_StrainAlleles.cpp):
  - the polymorphic sites of a copy (bases, indels, cover);
  - how a read stands at them;
  - the site shift (best allele, a known variant of another allele, a variable site, an indel);
  - the half-difference score;
  - the handler: a strain's read within 2 mismatches is settled on its species (site shift -12). A read of the other
    species' own gene, 10 mismatches clear, keeps its scores, with the shift pending.
  - the profiler's features on a reference with species neighbours, against hand counts (1 and 3 threads).
- **Script tests** ([test_insilico_strains.py](../../scripts/test_insilico_strains.py), `CongenerSites`):
  - a species whose congener differs at every 20th marker base puts 0.5 ± 0.06 of its marker substitutions there with
    `--congener-share 0.5`, and under 0.05 with 0;
  - the genes keep their divergence;
  - no stop codon is made;
  - the summary names the congener.
- **Feature sets** ([test_model_pmml.py](../../scripts/test_model_pmml.py)): the default set and named sets with
  `polymorphic`, the diagnostics in no group, and a clear error for a table of 2026-10-08.
- **Pipeline** ([test_gtdb_pipeline.py](../../scripts/mini_db/test_gtdb_pipeline.py)):
  - the species clouds come before the in-silico strains on the console;
  - `insilico_strains.log` measures a congener share, or says there were too few real-strain substitutions to measure
    on. It never says the clouds were missing.
- **CI suites** (WSL, 4 cores): mini database OK, end to end 146 OK (its database has strain alleles, so the settling
  ran), mock community OK, scripts 84 OK, pipeline and GTDB builds 83 OK. The one failure on the first run was the new
  console assertion reading `build_metadata.tsv`, which truncates the note; it now reads the step's log, and `test_a`
  passes.

Commands (scratchpad scripts; every job `taskset -c 0-3 nice -n 5`):

    bash gaps_test.sh                     # sync to ~/gaps, ninja -j4, protal_tests
    python3 -m unittest scripts/test_insilico_strains.py scripts/test_model_pmml.py scripts/test_ancestry_sites.py \
        scripts/mini_db/test_gtdb_build.py
    python3 -m unittest -k test_a_build_rerun scripts/mini_db/test_gtdb_pipeline.py

## What the next r226 build should show

- **The HPC checkout:** it must hold this branch and `504de80`. `build_metadata.tsv` names the commit.
- **`insilico_strains.log`:**
  - the congener share measured on the real strains;
  - the substitutions placed.
- **The per-sample `strain alleles:` lines:**
  - how many reads were unsure;
  - how many changed species (the strains rescued, or congeners' reads moved).
- **The trainer's importances:**
  - of `polymorphic_known_share`, `polymorphic_novel_share` and `ancestry_fixed_gain`;
  - whether `ancestry_agreement` loses rank to `ancestry_fixed_gain`.
- **The ancestry report's congener-site AUC for the in-silico strains:** v19 0.71-0.77, expected nearer the real
  strains' 0.36-0.47 now that their sites are realistic.
- **F1 against v18/v19:**
  - the FN of real strains (69% of pe's FN at the knob in v18);
  - the gut hold-out's 50 FN.
- **What `allele_sites_per_kb` alone separates:** the cluster-size check, as for `allele_copy_share`.
