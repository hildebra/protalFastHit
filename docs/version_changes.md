# What changed from protal 0.6.0a to 0.8.4

The shipped 0.6.0a (tag `0.6.0a`, commit `014f4a9`, August 2026), the 0.7 series (0.7.0 on
2026-09-30 to 0.7.9 on 2026-10-07), 0.8.0 (2026-10-09), 0.8.1 and 0.8.2 (both 2026-10-11), 0.8.3 (2026-10-13) and 0.8.4 (2026-10-14, branch `audit-fixes`) are compared here:
what changed, and how every version scores on the same simulated data in detection, abundance,
strains, speed and memory.
The numbers come from the benchmark and audit reports in [`docs/claude/`](claude/README.md), which
hold the commands, scripts and raw tables; this page only collects them. Nothing here was run again
for this page.

| Version | Commit | Date | Commits since the one before | In short |
|---|---|---|---|---|
| 0.6.0a | `014f4a9` | 2026-08-10 | | paired-end reads only, one model, raw database files |
| 0.7.0 | `4b21427` | 2026-09-30 | 98 | the audits' fixes; single-end and long reads; the single-file database; twice as fast; a database build-and-train pipeline |
| 0.7.1 | `1c11a00` | 2026-10-01 | 30 | mates paired across neighbouring genes; the reads' other candidates as model features; depth margin 0.08; genes at 2 bits per base |
| 0.7.2 | `fd1e719` | 2026-10-02 | 21 | divergence and conservation features; one sample profiled on all threads; `gene_congeners.tsv` |
| 0.7.3 | `b146951` | 2026-10-02 | 17 | one binary for every CPU; gene neighbours from every genome; strain rows from long reads; HiFi training reads |
| 0.7.4 | `8b2e376` | 2026-10-03 | 23 | species calls judged against the sample's congeners (relatives, distance, read EM); long reads through their chain; prefetched lookups |
| 0.7.5 | `b892133` | 2026-10-03 | 17 | against false positives: the sample's depth and the reads' divergence as default features, suspect gene copies, species priors; the index packed in memory |
| 0.7.6 | `13aac36` | 2026-10-05 | 17 | GTDB-scale speed (paired-end run 131 → 58 s, HiFi 104 → 25 s on r226); outputs independent of record order; species priors opt-in; reduced marker sets; in-silico strains and denser depths in training |
| 0.7.7 | `f2bafd9` | 2026-10-05 | 33 | GTDB-scale speed and memory (paired-end run 58 → 42 s, HiFi 25 → 15.6 s, peak 38.1 → 34.2 GB on r226): AVX2 flex scan, unaligned reads counted in the SAM header (`--write_unmapped_reads`), the index packed key by key as it loads; read EM capped at 100 sweeps; `.fq.zst` input |
| 0.7.8 | `5de3224` | 2026-10-06 | 16 | GTDB-scale speed (paired-end run 42 → 33 s, HiFi 15.6 → 12.3 s on r226, the same outputs): AVX2 tie masks in the seed lookup, flat gene tables, the k-mer screen from packed genes and reads, only shared seeds sorted; gzip with ISA-L (paired-end runs no longer wait for their input); `--add_model` in place; gradient-boosted models by default |
| 0.7.9 | `ffbedb3` | 2026-10-07 | 37 | against false positives in complex communities: the reads' consistency, the genes' shape, the database's neighbourhood, the sample's complexity and the ancestry sites (which side a read takes where the species differs from its congeners) as default features, untested at r226 until the next build; strain MSAs packed (159 → 37 MB per dense sample), spilled to disk and merged over runs; builds keep the reads behind every model error (`--error-reads`, `--share-logs` with the ancestry report), calibrate the in-silico strains' dN/dS on real strains, vary the scenarios and choose the knob by bootstrap; `simulate_metagenomes` makes Illumina, Ultima and long reads itself and streams large samples; `--build` resumes, packs the index straight into `database.protal`; every test suite in CI |
| 0.8.0 | `a82df8f` | 2026-10-09 | 39 | strains told from congeners: the species' other genomes as alleles of its marker genes (in the scores of the reads protal is unsure about, and as features), where the species varies, the congeners' gene copies (`congener_gaps.tsv`), a divergent read's untried candidates, ancestry sites by the congeners' consensus and across indels (85 default features, untested at r226); profiles end with the share of the genomes the called species do not explain (`?`); Nanopore candidates refused by an exact indel bound (−45% of the alignment's instructions), AVX-512 seeding, BGZF input on several threads; builds without the foreign leak, species complexes held out whole, ~118 output files instead of ~11,900, the build on many cores |
| 0.8.1 | `b3ac271` | 2026-10-11 | 22 | 0.8.0's features measured at r226 (v21-v23: strain alleles and ancestry sites +0.001-0.002 of F1 each, the site groups together +0.001-0.003) and acted on: the ancestry sites prebuilt by `--build` (profiling −47% of its instructions, the same profiles) and counted where the species' known genomes share them; the columns' conservation in the family from GTDB's own alignments (`weights`, 91 default features); up to 16 strain alleles per gene, the finished database's from every genome; the deep alleles' index entries off; more real strains and fewer in-silico ones in training; NCBI species taxids in the database and the profiles (`--ncbi_taxids`) |
| 0.8.2 | `53aa78c` | 2026-10-11 | 15 | start-up without the 12-s wait at r226 (the profiling's tables beside the alignment, column weights and strain alleles as binary tables, every table read once and only where used; the same outputs); WFA2 faster with every alignment the same (−26% of its instructions for paired-end, −39 to −48% for long reads); WFA2's Singletrack traceback, the default of `--build`; the weighted site shift on by default |
| 0.8.3 | `23bacd9` | 2026-10-13 | 31 | sparse-first gates per read type, recorded by the database and followed by runs (at r226 −35% of the paired-end alignment's time at C = 1, −62% of HiFi's at C = 8); GTDB builds gated by default, `--test-c-gates`; gzip input in pieces on several threads (system ISA-L); a sample's lanes as one sample; F1 over the species the marker genes see; the weighted site shift off again |
| 0.8.4 | `e990c69` | 2026-10-14 | 35 | the stride pick by default with gates pe:4, se:4, pb:16, ont:8 and models at other gates (r226 v27); the identity against each read's best strain allele; features that do not change with the read length or a sparse index (`length_free`, `read_length`, `ZS`; at the knob a sparse index at C 4 calls as ungated, 75 bp reads still call more false positives); GTDB builds on 75-250 bp paired-end, single-end Illumina and Ultima reads, random forests by default, both kinds of models in the database; opt-in a grid sparse index and a compact key map; zstd and many-member gzip input on several threads |

## Since 0.8.4

- **Unique k-mers counted from the alignment** (2026-10-14). `--build` stores which of every gene copy's index cores
  are unique (`unique_cores.bin`, two bits a core, ~0.7 GB at r226; `--no_unique_cores`: none), and a run counts a
  record's `ZS`, `ZU` and `ZT` from the gene's cores its alignment covers with exact k-mers instead of from its
  anchor's seeds ([databases.md](databases.md), [features.md](features.md#independent-of-the-read-length-and-the-seeding-length_free-read_length-2026-10-14)).
  The counts are those the whole index's seeds give, also on a sparse index or with `--sparse_seeding`, so every
  unique-seed feature is the same there and no seed expectation is written. A database without the table, or
  `--no_unique_cores`, counts the anchors' seeds as 0.8.4 did. The anchor moves to a new tag, `ZL` (1 a single seed, 0 a
  rescued mate). Models trained on 0.8.4's counts should be retrained.
- **The unique k-mer rates back in the default features** (2026-10-14). With the counts from the alignment the eight
  rates (`lu_per_kb`, `lsu_per_kb`, the six gene rates) no longer change with the seeding. `length_free` now puts the
  two unique seed enrichments beside them instead of in their place: 98 default features (90 before).
  - **Why:** on the strain world, models with both were best at equal false positives and called the fewest false
    positives at the knob. Their sparse index calls as the whole index does: 0.903 against 0.847 for the previous set
    with the anchors' counts ([report](claude/2026-10-14-length-free-features/README.md)).
  - **Cost:** the rates still fall with 75 bp reads (×0.79: fewer whole k-mers fit a short read).
- **Ancestry sites in ~0.65 of the memory** (2026-10-14). `ancestry_sites.bin` was 1.13 GB at r226 v27, held whole
  beside the alignment. Version 2 stores the positions as Rice-coded gaps, the index in 4-byte keys and offsets and the
  identity in 16 bits; on the reduced r226 database's table that is 0.65 of the bytes (0.69 on the small real world; ~0.73 GB at r226), with the same sites. A run reads
  either version; the run log gives the table's version and size.
- **Failed candidates over anchors of two seeds on both sides** (2026-10-14). `failed_candidate_rate` and
  `failed_gene_share` leave out the reads anchored on a taxon by single seeds only (`ZL:i:1`), as their failures already
  were ([features.md](features.md#independent-of-the-read-length-and-the-seeding-length_free-read_length-2026-10-14)).
  - **Whole index:** no single-seed anchors, so its rates do not change.
  - **Grid's sparse index at C 4:** an absent species' mean rate on the strain world rises only from 0.155 to 0.163
    (0.340 on the whole index), and its calls do not change. The rest of the fall is in the failures themselves: a
    relative's read meets the taxon there on a single grid seed.
- **protal compiles in parallel** (2026-10-14; building from source only: protal runs the same). `main.cpp` was all of
  protal in one translation unit, 97-125 s on one core and compiled again whole after any edit. A run's stages are now
  units of their own in `src/Run/`, and the options' parsing is in `src/Options.cpp`
  ([development.md](development.md#source-layout)); cPMML is built without `-flto`, which cost every link 4.4-5.2 s.
  - **Rebuilds on 4 cores:** an edit within one stage 24-36 s (was 103-187 s), an edit every unit sees (`Options.h`)
    53-111 s (96-149 s), from scratch 86-90 s (105-157 s) ([report](claude/2026-10-14-compile-time/README.md#implementation)).
  - **Cost:** the units together do 1.6 times the old compile's work, and the binary's code is 12% larger. The
    outputs and the alignment's and profiling's instructions are the same (within 1%).
- **The mates of a pair vote on its taxon** (2026-10-15). `--mate_guidance consensus` is the default; `guide` is the
  guidance of 0.7.1-0.8.4, and `off` (`--no_mate_guidance`) is neither
  ([running.md](running.md#alignment-options)).
  - **The vote:** a pair whose mates did not align together as its best candidate is scored on every taxon either
    mate aligned to by the sum of both mates' scores. A mate without an alignment there is looked for, and where it
    is not found it scores as an alignment at the identity floor. The MAPQ comes from the margin over the next taxon.
  - **Why:** guidance completed the pair on a sure mate's taxon without weighing the other mate's own alignment, and
    doubled its MAPQ (audit L1).
  - **Effect, with the models trained on `guide`:** on the reduced r226 database's simulated test sets the calls are
    the same within 1-2 species (its own test set: F1 0.9597 to 0.9618), and on six paired-end mocks they are
    identical. Fragments counted for a wrong species fall by 0.3-0.4%, those for the right one by 0.02-0.05%.
  - **Cost:** +2.65% of the paired-end alignment loop's instructions against `guide` (200,000 pairs of a real sample).
  - **Databases:** a database should record the mode its models were trained with (`build_settings.tsv`); the
    databases built up to now used `guide` ([report](claude/2026-10-15-mate-consensus/README.md)).

## What 0.8.4 adds

0.8.4 collects the 34 commits after 0.8.3. The changes:
- **Presence models:**
  - the identity against each read's best strain allele;
  - features that do not change with the read length or a sparse index (`length_free`, `read_length`, and a record's
    `ZS`);
  - both kinds of models in the database (`--model_kind`), random forests their default again.
- **Read alignment:** the stride pick by default, with gates pe:4, se:4, pb:16, ont:8 and models at other gates. Opt-in:
  a sparse index by the genes' grid, a compact key map, `--flex_indels` and `--sparse_seeding`.
- **Read input:** zstd files and gzip files of many members read on several threads; each sample's log says how its
  input was read.
- **Database builds:**
  - training on paired-end reads of 75 to 250 bp, and on single-end Illumina and Ultima reads;
  - `--test-c-gates` also tests the stride pick and `--flex_indels`;
  - the database file holds its notes.
- **Mock benchmark:** scripts that download public mock communities and score protal on them (`scripts/mock_benchmark`).

The stride pick and its gates come from the r226 v27 build. The rest is measured on reduced or small databases only.

- **Identity against the best strain allele** (2026-10-13). Two new features, `allele_identity` and
  `allele_top_identity`.
  - **What they are:** `identity` and `top_identity` with each read measured against its copy's best strain allele
    (the representative's copy where no allele beats it). They form the group `strain_identity`, in the default set:
    93 features ([features.md](features.md#the-identities-against-the-best-strain-allele-strain_identity-2026-10-13)).
  - **Why:** a strain far from its representative but near a known genome of its species looked like a relative. HiFi
    reads of *E. coli* K-12, 1.6% from GTDB's type strain, had `top_identity` 0.9957 where the sample's other species
    had 1.0, and the HiFi model called a congener instead.
  - **Effect** on a reduced r226 database rebuilt with them
    ([report](claude/2026-10-13-consensus-allele-identity/README.md)):
    - both HiFi mocks call *E. coli* (P 0.12 → 0.87 and 0.97), and no longer the congener;
    - F1: ATCC MSA-1003 0.894 → 0.974, Zymo D6331 0.866 → 0.933;
    - the simulated test sets stay within noise.
  - **Databases:** the features need new models. An older database's models do not use them.
  - **Tried and left out:** a tie-break for long reads whose taxa beat each other in a cycle. It moved more reads to
    relatives than back.
- **No empty output folders** (2026-10-13). A map run made `alignments/`, `profiles/`, `strains/` and `misc/` as
  it read the map, and left the unused ones empty: `alignments/` and `profiles/` with a map without `SAM` and
  `PROFILE` columns (whose files go to the output folder itself), `strains/` with `--no_strains`. protal now makes a
  folder only when it writes into it ([running.md](running.md#where-the-outputs-go)).
- **The stride pick by default, gates 4/4/16/8, and models at other gates** (2026-10-13). The r226 v27 build
  tested the gates and both picks ([report](claude/2026-10-13-r226-v27-c-gates.md)): up to pe/se C 4, pb C 16 and
  ONT C 8 no measurable F1 was lost, for 36-44% less alignment CPU than at C 1 (pb 11% less than C 8), and the stride
  lost 2.5-17 times fewer present species than the hash, at equal or better F1.
  - `--build` records the stride pick for every read type unless `--sparse_first_pick` says otherwise. Runs follow
    the database's pick, stride where it records none; a database built on 2026-10-12/13 that records gates but no
    pick stays hash ([running.md](running.md#developer-options)).
  - A database may hold models trained at other gates (`model_pe_c1.xml`, ...; `--add_model --model_gate C`). A
    sample is profiled with the model trained at the gate nearest its own and warned only if that is not its gate.
  - `build_gtdb_database.py` builds at `--sparse-first pe:4,se:4,pb:16,ont:8` with `--sparse-first-pick stride` and
    trains models at `--model-gates pe:1/4/8,se:1/4/8,pb:8/16/32,ont:1/8/16`: both collections are profiled at each
    gate too (the streamed simulations streamed again), and the models go into the finished database. Each other
    gate costs the collections' alignment once more ([databases.md](databases.md#choosing-the-sparse-first-gates)).
  - `c_gates.tsv` counts the present species a gate leaves without a row as missed (`unseen`). Before, they were
    left out, so that the hash, which loses more, looked better. Each row is scored by the model of the gate nearest
    it (`model_gate`), and the pick is a column (`pick`). `--test-other-pick` (`--test-stride` before) tests the pick
    the database does not use, by default at its own gate; `--test-c-gates`' default is `pe:0/2/16,se:0/2/16,pb:0/4,
    ont:0/4/32`.
  - `measure_performance.sh --test_sparse` gives `--sparse_first_pick hash` with `--sparse_seeding`, which stands
    for an index of the hash pick's cores.

- **The gene grid and a sparse index** (2026-10-13; opt-in, not measured at GTDB scale). `--build` gives every index
  entry a grid level (0-5, three bits of its position): of each gene's cores, one per window of C consecutive cores
  is kept at C, the one a read of the species most likely carries intact (by its strain alleles, else by its columns'
  conservation within the family's genera) and the least shared with other genes. The pick `grid`
  (`--sparse_first_pick grid`) looks up every core of a read and counts only the grid's entries of its gate, about
  every C-th core of the stretch the read covers, and seeds a read kept from those alone; a short read's gene with one
  grid seed anchors on it (`--no_single_seed_anchors`: on two, as other runs). On the whole index (`--sparse_index off`,
  or samples not all grid-picked) the grid only gates, and a read kept is seeded from every core on two seeds a gene;
  `--sparse_index emulate` seeds as the sparse index on the whole index (the build's collections use it). A run
  whose every sample is grid-picked loads only those entries (`--sparse_index auto`): about 1/C of the values, so
  at r226 an estimated ~6.7 GB of values instead of 26.6 GB at C 4 (the key map stays 3.2 GB). On the mini database
  a paired-end run at C 4 wrote 0.05% fewer records than ungated (the stride on the whole index 0.23%, the hash
  pick's emulated sparse index 6.5%; without single-seed anchors 3.2%) for 1-5% more candidate alignments; HiFi lost at
  most 0.1% up to C 8 and 2% at C 32 ([report](claude/2026-10-13-gene-grid-index.md)).
  `build_gtdb_database.py --test-grid` (on by default with `--test-c-gates`) scores the test samples with the grid
  pick beside the database's, and `--sparse-first-pick grid` builds a grid-picked database. A database built before
  2026-10-13 has no grid levels: the grid pick stops on it (exit 2); the other picks are unchanged.

- **Presence features that do not change with the read length or a sparse index; training on several read lengths;
  random forests again** (2026-10-14). On a grid's sparse index the reads align as on the whole index, but a read
  seeded from 1/C of its cores counts 1/C of its unique k-mers, and a model trained on fully seeded runs called fewer
  strains (0.867 against 0.900 on the strain world); the same reads cut to 75 bp made a 150 bp model call six times the
  false positives (MAPQ, and relatives' short reads without a difference) ([report](claude/2026-10-14-sparse-index-calls/README.md)).
  Now:
  - A record's new `ZS` tag gives the seeds of its anchor that `ZU` and `ZT` were counted over (a single-seed anchor's
    seed counts too); a sample seeded from a grid level alone has that level's seed expectation in its SAM header
    (`@CO protal seed expectation: level=L unique=U two=T`, also logged).
  - The feature groups `length_free` (the reads' unique seeds over those their genes lead to expect, in place of the
    unique k-mer rates; the best reads' identity bounds by their aligned bases, in place of `top_identity` and
    `allele_top_identity`; the ancestry sites per kb) and `read_length` (the sample's mean and SD of aligned bases,
    and the taxon's against them) are in the default set ([features.md](features.md)). Their columns are in every
    training dump; a table of an older protal trains with the previous set.
  - A single-seed anchor's failure is no failed candidate (`ZF`) any more: it tripled `failed_candidate_rate` on the
    sparse index.
  - `build_gtdb_database.py` trains paired-end reads of 75, 100, 150 and 250 bp (`--read-setups`; 9 samples a point,
    the deep points rebalanced), single-end Illumina reads of 150 and 250 bp (`--se-read-setups`) and Ultima reads
    (`--se-ultima-setup`, `ultima:300:50:25:2`), and trains random forests as the databases' default models again
    (`--model forest,gbm`; `machine_learning_cmdline.py --model forest`); protal's `--model_kind auto` takes the
    database's default.
  - Simulated Ultima reads are cut at 320 bases (`ultima:...:MAX`), fall in quality over their last 100 bases, and
    their errors are indels with rare substitutions (Q35 to Q40), for scenarios and the design alike.

  **Effect** on the small real world and the strain world ([report](claude/2026-10-14-length-free-features/README.md)):
  - **The features hold:** they do not move with the read length or on the sparse index. A model trained on fully
    seeded runs calls on the sparse index at C 4 as it calls ungated (0.915 against 0.918 at the knob; the previous
    set 0.847). At equal false positives the previous set ranked as well there: the gain is calibration.
  - **Short reads still call more false positives:** 75 bp reads call 2 to 6 times those of 150 bp under most models,
    also those trained on 75 bp reads.
  - **New false positives:** the new set adds deep-sample false positives of taxa with a few fragments and no unique
    seed.

- **A compact key map** (2026-10-14; opt-in, `--compact_keys`). The key map gives every one of the 2^30 cores a
  16-bit offset (3.2 GB), of which r226 uses 9.2%. With the option it holds a presence bit per core with ranks (153
  MB) and 3 bytes a used core (~0.45 GB at r226), built as the index loads, so the direct key map is never held. On
  the reduced r226 database a run peaked at 0.8 GB instead of 3.4 GB with the same records; paired-end lookups cost
  the same within the noise, HiFi lookups ~7 ns more (+2.6% of the wall time) ([report](claude/2026-10-13-compact-key-map.md)).
  Not yet measured at r226.

- **How the input was read, in every run's log** (2026-10-13). After each sample's alignment a line says how each read
  file was read (format, members, threads; where it went on by members or fell back to one thread). It also says the
  share of the alignment its threads spent getting reads, a quarter or more being "the input held the alignment
  back", and the alignment and reading threads against the CPUs the run may use. `measure_performance.sh` writes it to
  `runs.tsv` (`input`, `reader_share`) and flags such runs; `--cold_reads` drops the read files from the page cache
  before each run.
  - Read files are asked to be read ahead (the next 64 MB); `PROTAL_READ_AHEAD=0` turns it off. On a local disk it
    changed nothing measurable; a network file system read cold is still to be tested.
  - A gzip file read in pieces holds 29-38% less memory, as fast (a shorter look-ahead: ~0.35 GB less for a
    paired-end run on 32 threads).
  - The reading threads beside a full `-t` cost the alignment nothing measurable (4 CPUs, `-t 4`, 1-4 reading
    threads per file).

  ([report](claude/2026-10-13-reading-input/README.md), [running.md](running.md#read-files))
- **gzip files of many members on several threads** (2026-10-13). 0.8.3 read a gzip file of many small members
  (members one after the other, not BGZF) in pieces at about a third of one thread's speed. Its block finder scanned
  the whole file for blocks to start a piece at and found none. At r226 the paired-end sample's alignment took 58 s
  instead of 24 s at 32 threads ([report](claude/2026-10-13-gzip-members-regression/README.md)). Such a file, its
  second member within 16 MiB of its start, is now read by members, as BGZF is: whole members inflated on 1 to 4
  threads per file. A header found by chance in a member's data, or a damaged member, sends the rest of the file to
  one thread ([running.md](running.md#read-files)). A file read in pieces goes on by members where small members
  follow one of its members (a large gzip file and a file of many members concatenated). A file read in pieces no
  longer keeps every page of the compressed file in memory: the pages behind the reader are released (the r226
  paired-end run held 7.3 GB more). The outputs are the same.
- **zstd read files on several threads** (2026-10-13). A zstd file's frames are independent, as BGZF blocks are:
  a zstd file of several frames is now decompressed on as many threads as a BGZF file is (1 to 4 per file, one per
  4 of `-t` and file). This covers the reads of GTDB builds, whose simulations write a frame per batch of reads,
  and `pzstd`'s files. A file of one frame (what `zstd` writes, also with `-T`) gains nothing; a frame whose header
  does not say its size (a stream's) and the rest of its file are read as a stream on one thread, as before
  ([running.md](running.md#read-files)). The outputs are the same; speed not measured yet.
- **`--test-c-gates` also tests `--flex_indels`** (2026-10-13). The ONT test samples at every gate and the Ultima
  ones at se's are also profiled with protal's `--flex_indels`, each beside the same samples at the database's gate
  (`c_gates.tsv`'s `flex_indels` and `baseline_*` columns; `--test-flex-indels`). The test set's Ultima simulations
  are then written to disk rather than streamed. This needs a protal with `--flex_indels`
  ([databases.md](databases.md#choosing-the-sparse-first-gates)).
- **`--test-c-gates` also tests the stride pick** (2026-10-13). Each read type's test samples are also profiled at
  each of their tested gates of 2 or more with the gate looking up every C-th core of a read rather than its hash
  sample (protal's `--sparse_first_pick stride`, per sample by the map's `SPARSE_FIRST_PICK`; `--test-stride`;
  `c_gates.tsv`'s `stride` column). With an older protal the build leaves these rows out, with a warning
  ([databases.md](databases.md#choosing-the-sparse-first-gates)).
- **The build's disk estimate counts the test samples at other gates** (2026-10-13). Each row a `--test-c-gates`
  point is profiled at writes its own SAMs and profiles: ~74 GB at r226 with the defaults, which the estimate left
  out. `--stream-above auto` now counts them for every point it writes.
- **Both kinds of models in the database, chosen per run** (2026-10-13). `build_gtdb_database.py` trains every model
  twice by default (`--model gbm,forest`): gradient-boosted trees, the database's default files as before, then random
  forests at the same gates, stored beside them (`model_pe_forest.xml`, `model_pe_c1_forest.xml`, ...; the trainers'
  memory as with one kind, the training time about doubled). Both are checked for parity with protal and scored side
  by side in `summary.txt`. `protal --model_kind forest` (or `gbm`, or per read type, `pe:forest`) profiles with that
  kind, each sample with the model of the gate nearest its own; `--profile_only` with it re-profiles a run's SAMs.
  `--add_model --model_kind` stores a model of a kind, checked to be one ([databases.md](databases.md#kinds-of-models)).
- **The database file holds its notes** (2026-10-13). `database.protal` also carries `build_metadata.tsv`,
  `gene2geneid.tsv` and `genome2tiid.tsv`, which no run reads, so the one file is the whole database: `--build` packs
  them from its folder and leaves them there, `--add_tables` stores them (with `--add_model` in one write), and the
  GTDB build's last step puts in the models and `build_metadata.tsv` together, in place
  ([databases.md](databases.md#the-notes-a-database-carries)). Every file a run reads was already in it.

## What 0.8.3 adds

0.8.3 collects the 30 commits after 0.8.2. The changes:
- **Read alignment:** the sparse-first gate gives up the reads that find nothing near-exact, recorded per read type
  by the database.
- **Read input:** gzip input inflated in pieces on several threads, and a sample's lanes read as one sample.
- **Database builds:** an F1 over the species the marker genes can see, and the gates trained on and tested.

None of it is measured in an r226 build yet; the first gated build is running.

- **The weighted site shift off by default again** (2026-10-12). The first r226 build trained with it
  (v24b) showed no measurable F1 gain, and its short-read species depths were 1-2% lower; on the small
  real world, with the same samples, the shift lowered 19-22% of the called species' depths
  ([report](claude/2026-10-12-r226-v24b/README.md), section 4). `--weighted_site_shift` turns it on, and
  `build_gtdb_database.py --weighted-site-shift` trains with it. r226 v24b was built with 0.8.2's default:
  profile it with `--weighted_site_shift` for the calls of its training.
- **A sample's lanes as one sample** (2026-10-12). A map's `FIRST` and `SECOND` cells list several files,
  comma-separated, and `--one_sample` makes the files of `-1`/`-2` one sample's: they are read one after the
  other as one stream, so the sample's depth and the features that depend on its size see all its reads
  (each lane a sample of its own before). A lane whose mate files hold different numbers of reads fails the
  sample. `protal_map_utils generate` makes Illumina's lane files of a sample one sample (`--split-lanes`: a
  sample per lane, as before) ([running.md](running.md#several-files-of-one-sample)).
- **F1 over the species the marker genes can see** (2026-10-12). `simulate_metagenomes --gene_positions`
  counts each genome's reads that overlap a marker gene; a GTDB build passes the release's positions, the
  training tables get `meta_marker_reads` and `<table>.unseen.tsv`, and the trainer and `summary.txt` add an
  F1 with only the species whose marker genes got reads present ([databases.md](databases.md#3-check-the-result)).
  At r226 v24b 60-76% of the false negatives had at most two reads of their own on the taxon
  ([report](claude/2026-10-12-r226-v24b-errors/README.md)). Untested at r226.
- **A soft fragment count, opt-in** (2026-10-12). `--fragment_count soft` (recorded by `--build` in
  `build_settings.tsv` and followed by the runs on the database) counts the model's `fragments` feature as the
  fragments the MAPQ and length filters keep plus each dropped one by its read-EM share, so that a divergent
  strain whose reads tie with a congener's reference is not counted as nothing;
  `build_gtdb_database.py --fragment-count soft` builds and trains with it. The default stays `mapq`. Untested
  at r226.
- **gzip files that are not BGZF on several threads** (2026-10-12, reworked 2026-10-13). A gzip read file of
  one member (as `gzip` writes it) is inflated by ISA-L in pieces on 1 to 8 threads per file (one per 4 of `-t`
  and file: 2 per file for a paired-end run on 16 threads, 4 on 32): pieces are decoded ahead twice with the
  32 KB before them as labels, which say where each copied byte comes from once that window is known, while
  the file's own thread decodes the others in order. On FASTQ about 2.0 and 2.4 GB/s on 2 and 4 threads
  against 1.7 on one (4 cores) ([report](claude/2026-10-12-isal-parallel-inflate/README.md)); BGZF input stays the
  fastest. It replaces rapidgzip (2026-10-12 only: 0.3 GB/s per thread on FASTQ,
  [report](claude/2026-10-12-rapidgzip-input/README.md)), and with it the copy of ISA-L in `lib/isa-l` and the
  nasm it needed: protal links the system's ISA-L again (libisal-dev, conda's isa-l), as before. The outputs
  are the same.
- **Sparse-first gates recorded by the database** (2026-10-12). `--sparse_first` takes a gate per read type
  (`pe:1,se:1,pb:8,ont:1`), `--build` records them in `build_settings.tsv`, and runs follow the database's
  gates unless told otherwise (`auto`, the default; a map's `SPARSE_FIRST` column gates one sample), saying
  which they use and warning when a sample is gated otherwise than the models were trained.
  `build_gtdb_database.py` builds with `--sparse-first pe:1,se:1,pb:8,ont:1` by default (at r226 paired-end
  C = 1 cost 0.01% of the records for 35% less alignment time, HiFi C = 8 0.54% for 62% less;
  [report](claude/2026-10-12-sparse-first/README.md)), so its models are trained on gated runs, and
  `--test-c-gates` profiles the test samples at other gates too and writes their F1 and alignment CPU to
  `model_logs/c_gates.tsv` ([databases.md](databases.md#choosing-the-sparse-first-gates)). Databases built
  before stay ungated. Untested at r226.
- **Smaller changes**, all with the same outputs:
  - `--stage1_probe` measures what a gate would keep.
  - At C = 1 a kept read is seeded from the gate's own lookups.
  - Every line of a lookup's key-map block is prefetched.
  - The profiling's caches are split in 64 locked shards (a 32-thread run's SAM parse 119 → ~25 thread-seconds).
  - `--scratch_dir` puts a run's intermediate files on a local disk; it replaces `--strain_spill`.
  - Build console: GTDB builds say each protal run of a collection and each evaluation as it starts.

## What 0.8.2 adds

0.8.2 collects the 14 commits after 0.8.1: start-up, the aligner and its speed. No r226 build or cluster run
of 0.8.2 has been measured yet; the numbers below come from r226's archaeal species, the small real world
(six GTDB r226 families) and protal's own WFA2 calls replayed.

### For everyone running protal

- **No 12-s wait before the alignment.** At r226, 0.8.1 took 51 s for a paired-end sample that 0.7.8
  took 33 s for; 12 s of it went to parsing `column_weights.tsv` on one thread before any read was
  aligned. The tables only profiling reads now load beside the alignment; the column weights and
  strain alleles are also binary members of `database.protal` (`column_weights.bin`,
  `strain_alleles.bin`), read without parsing; every table is read once, on several threads; and a
  table is read only where the run uses it. On all r226 archaeal species the column weights load in
  622 M instructions instead of 3,018 M and the strain alleles in 134 M instead of 587 M, with
  byte-identical outputs ([report](claude/2026-10-11-performance-audit/README.md)). A run with
  `--no_profile` now aligns with the column weights as a profiled run does (its SAM differed).
- **WFA2 faster, every alignment the same.** Four exact changes to WFA2-lib (loops that vectorize,
  fewer end checks, AVX2 copies chosen at run time), a reset of its wavefront memory that touches only
  what an alignment used, and the penalties halved (2/3/1) wherever X-drop cannot prune. On protal's
  own WFA2 calls that is 26% fewer instructions for paired-end reads, 46-48% fewer for PacBio and
  39-43% fewer for Nanopore, with the same outcome, score and CIGAR for every call; about −6% of
  short-read and −23 to −40% of long-read alignment are expected at r226
  ([report](claude/2026-10-11-wfa2-singletrack-speed/README.md)).
- **The Singletrack traceback** (`--wfa_traceback`; López-Villellas et al. 2025) keeps only WFA2's
  match wavefronts and recomputes a gap's steps while tracing back. It gives the same scores as the
  classic traceback but may choose otherwise among equally good alignments: 0.01% of paired-end,
  0.2-0.6% of PacBio and 47-49% of Nanopore records differed. With the changes above it needs
  1.6-2.4% fewer WFA2 instructions than classic and 8-20% less time. `--build` aligns with it by
  default and records it in `build_settings.tsv`; a run aligns as its database was built (classic for
  a database built before) and is warned when told otherwise
  ([report](claude/2026-10-11-singletrack/README.md)). WFA2-lib is now the fork that has Singletrack
  (tmokveld/WFA2-lib), with three fixes of protal's.
- **The weighted site shift on by default.** An unsure short read's site shift also discounts its
  differences at the columns the family's genera change freely (`column_weights.tsv`): neutral to
  +0.003 of F1 on the small real world, for ~1% of the alignment's instructions
  ([report](claude/2026-10-11-r226-v23/README.md), section 8). `--no_weighted_site_shift` turns it off.

### For anyone building a database

- `--build` aligns with Singletrack unless `--wfa_traceback classic`, and writes `build_settings.tsv`,
  `column_weights.bin` and `strain_alleles.bin`; `--compress_db` adds the binary tables to an existing
  database. `build_gtdb_database.py` trains with the weighted site shift (`--no-weighted-site-shift`
  leaves it out) and passes `--wfa-traceback` on to both builds.
- `measure_performance.sh --sylph SYLPH_DB` profiles the same reads with sylph for comparison;
  `--threads N` sets every run's threads.

### Compatibility

- 0.8.2 reads 0.8.1, 0.8.0 and 0.7 databases. Without the binary tables a run parses the text ones
  (`protal --compress_db` adds them); a database without `build_settings.tsv` was built with the
  classic traceback, and runs on it align classic.
- A database built with 0.8.1 or earlier (r226 v23 and v24 among them) was trained without the
  weighted site shift: profile it with `--no_weighted_site_shift` for the calls of its training.
- A database built with 0.8.2 trains its models on Singletrack alignments, whose Nanopore records
  differ from a classic database's in their choice among equally good alignments. No r226 build with
  Singletrack has been evaluated yet.
- `build_gtdb_database.py` refuses a protal or `simulate_metagenomes` binary of another version.
  Rebuild both after updating a checkout.

## What 0.8.1 adds

0.8.1 collects the 21 commits after 0.8.0. Three r226 builds have measured 0.8.0's features since,
with each model refitted without one feature group at a time. With species held out, the strain
alleles gained 0.0010-0.0017 of F1 and the ancestry sites 0.0010-0.0022 in v21
([r226 v21](claude/2026-10-09-r226-v21/README.md)); the polymorphic sites gained nothing
measurable there. In v22 the column weights gained nothing, two of their five features being
constant, and the deep alleles' index entries too little for their memory
([r226 v22](claude/2026-10-10-r226-v22/README.md)). v23 (at `8d8c89b` with `--index-alleles 0.01`,
the weights fixed, before the prebuilt ancestry sites, 16 alleles and the NCBI table) scored 0.0022-0.0037 above v22
with species held out (paired-end 0.9595, single-end 0.9650, PacBio 0.9791, Nanopore 0.9741).
Left out together, the three site groups (`ancestry`, `polymorphic`, `weights`) cost it 0.0012-0.0031;
each alone is replaceable by the others ([r226 v23](claude/2026-10-11-r226-v23/README.md)).

### For everyone running protal

- **Profiling at half the instructions, the same profiles.** `--build` computes every gene copy's ancestry
  sites once and stores them (`ancestry_sites.bin`). A run decodes them instead of computing them on
  a copy's first record, which was half of a sample's profiling instructions at r226 v23. A 1M-pair
  sample's profiling took 6.47 → 3.45 G instructions; its profiles are byte-identical. A database
  without the table computes the sites as before (`--no_ancestry_sites` at build time).
- **NCBI taxids.** The converter writes each species' NCBI species taxid, the one most of its GTDB
  genomes have (`ncbi_taxids.tsv`), and `--ncbi_taxids` adds the taxid and NCBI's name to every
  profile line. At r226 17.5% of the species get a named NCBI species and 82.5% one of NCBI's
  placeholders ("uncultured Lachnospiraceae bacterium"), mostly shared with other GTDB species.
  A database converted before 0.8.1 gets the table without a rebuild
  ([databases.md](databases.md#ncbi-taxids), [running.md](running.md#options-the-website-does-not-list),
  [report](claude/2026-10-11-ncbi-taxids.md)).
- **More strain alleles.** `--build` keeps up to 16 alleles per species and gene (4 in 0.8.0). On a
  strain-rich world of GTDB r226 marker genes 16 called 6 more test strains of 600 than 4, and 1
  fewer; at r226 they should cost about 0.9 GB of memory (an estimate)
  ([report](claude/2026-10-11-strain-false-negatives/README.md)). The deep alleles' k-mers can go
  into the index (`--index_alleles 0.01`), off by default: at r226 v22 they added 1.1 GB to the
  database and 2.8 GB to the run's memory for about 0.0025 of short-read F1, confounded with other
  changes of that build.
- **Ancestry sites the species shares.** The ancestry features count only the sites that the
  species' known genomes (its strain alleles) share, not the representative's own recent mutations,
  at which a real strain carries the congeners' base. On six real GTDB families a quarter of the real
  strains had agreed at under 88.5% of the sites ([real ancestry](claude/2026-10-10-real-ancestry/README.md),
  [small real world](claude/2026-10-10-real-small-world/README.md)). From a world grown on known trees
  ([ancestry true-positive test](claude/2026-10-09-ancestry-true-positive-test/README.md)): the
  agreement at the sites the species does not vary at is a feature of its own (`ancestry_fixed_agreement`),
  the alleles are chosen for the coverage they give a sample instead of farthest first, and from six
  congeners compared one with a base of its own no longer breaks a site's consensus.
- **The columns' conservation** (`weights`, five features, 91 in the default set): how conserved each
  column of a marker gene is within the genera of the species' family and among them, and the amino
  acid per codon (`column_weights.tsv`, `--build --column_weights`, 10 genera per family by
  default). The columns come from GTDB-Tk's own alignment and GTDB's tree where the converter wrote
  them (`--column_weights_alignment msa`, the default where present; 99.7-99.9% of GTDB's aligned
  residues placed on real families), else from the proteins aligned. Mismatches at conserved columns
  are read errors or reads from afar; a relative's mutations avoid them. The family's consensus
  also polarises the ancestry sites of genera with fewer than three species. `--weighted_site_shift`
  (opt-in, untested at r226) lets a difference at a weakly conserved column count half in the
  scores of the reads protal is unsure about ([features.md](features.md#the-columns-conservation-weights-2026-10-10)).

### For anyone building a database

- **Real strains over in-silico ones.** r226 v23 simulated real strains for 13% of its species and
  in-silico strains for 40%, which the models miss less often than real ones.
  `--insilico-strains auto` (the default) now makes only as many in-silico strains as keep the real
  ones at half of the simulated strains or more (`--min-real-strains`), and `download_gtdb.py` fetches
  up to 4 strains of every species that has them (`--species 40000 --per_species 4`; 36,515 species
  at r226). With the defaults, real strains should be ~54% of the simulated strains.
- **The finished database knows every strain.** Its alleles come from all of a species' genomes
  (`--final-allele-genome-share 1`); the training database's from half, the genomes no sample is
  simulated from (`--allele-genome-share 0.5`).
- The converter writes GTDB's alignments and trees for the column weights (`column_msa/`,
  `species_tree.nwk`; `download_gtdb.py` fetches them, 2.2 GB at r226, unless `--no_msa`) and
  `ncbi_taxids.tsv`. `--weighted-site-shift` and `--index-alleles` pass the run options on.
- The ancestry true-positive test (`simulate_ancestry_world.py`, `ancestry_oracle.py`,
  `ancestry_truth_test.py`) checks protal's ancestry, allele and polymorphic counting against a
  Python port on worlds grown on explicit trees.

### Compatibility

- 0.8.1 reads 0.8.0 and 0.7 databases. Without `column_weights.tsv` the `weights` features are −1;
  without `ancestry_sites.bin` a run computes the sites (the same features, slower); without
  `ncbi_taxids.tsv`, `--ncbi_taxids` stops before aligning
  (`gtdb_to_protal_db.py --priors_only` and `--add_tables` add it).
- The ancestry features count other sites than in 0.8.0: retrain a model trained on 0.8.0's
  tables. A training table of 0.8.0 lacks the `weights` columns: train it with the groups it has.
- `build_gtdb_database.py` refuses a protal or `simulate_metagenomes` binary of another version.
  Rebuild both after updating a checkout.

## What 0.8.0 adds, and why it matters

Since the r226 v15 build, the GTDB models' remaining errors are mostly of two kinds, both at the
same marker identity. A missed species is usually a strain of a species the database knows from
other genomes (60-80% of v18's false negatives). A false call is usually the nearest congener of a
species the database lacks (86-91% of v17's false positives). Nothing in a read's alignment to the
reference tells the two apart ([error-read signatures](claude/2026-10-07-error-read-signatures/README.md),
[r226 v17](claude/2026-10-08-r226-v17/README.md), [r226 v18](claude/2026-10-08-r226-v18/README.md)).
0.8.0 therefore gives protal what the species' other genomes and its congeners' gene copies say.
These features are in the default set, but no r226 build with them has been evaluated yet. Their
gain is still to be measured, and a group that does not pay will leave the set, as `foreign` did.
0.8.0 also adds a new profile line (the unknown share) and more speed.

### For everyone running protal

- **Profiles end with the unknown share.** A profile's last line, with `?` for name and taxonomy,
  gives the share of the sample's genomes that the called species do not explain. The species'
  abundances are now shares of all genomes, not of the called species alone. The share comes from
  each species' genome size (its GTDB genomes' mean, corrected by CheckM, in `species_priors.tsv`),
  its depth, and the bases read, which the aligner now writes into the SAM header. `<profile>.composition` adds the
  share of the reads explained, the average genome size, and how many species the rest would be.
  *Why it matters:* abundances among the called species alone overstate every species when much of
  a sample is not in the database (novel species, host, species too shallow to call). The `?` says
  how much that is, and makes abundances comparable between samples. *How accurate:* on a synthetic
  world with every species in the database, 0.965-0.976 of the reads were explained (truth 1). In
  the r226 v19 test samples, `?` was within ±0.02 for short reads from 50,000 fragments up. It was
  too low in shallower samples, where the called species take the reads of uncalled ones, and too
  low by 0.11-0.15 for PacBio at every depth (PacBio depths 5-35% high, cause open). Host reads count
  as unknown genomes. `--no_unknown_share` gives the old profile. A database converted before
  2026-10-08 has no genome sizes, so its profiles have no `?` line until `gtdb_to_protal_db.py
  --priors_only` and `protal --add_tables` give it them. Scripts that read profiles should expect the
  `?` line; `protal_profile_utils merge` keeps it last
  ([running.md](running.md#what-the-called-species-explain-the-unknown-share),
  [report](claude/2026-10-08-sample-composition/README.md),
  [r226 v19](claude/2026-10-09-r226-v19/README.md)).
- **Strains of known species.** With `--full_reference`, `--build` keeps up to four alleles of each
  species' marker genes from its other genomes (`strain_alleles.tsv`: edits of the representative's
  copy, each nearer it than the nearest congener's copy). protal uses them for the short reads it is
  unsure about, those with a candidate of another species within 3 mismatches of the best. There,
  each candidate is scored against its species' alleles: a difference a known strain carries counts
  as a match, and one where the species varies counts as half a difference (`--no_allele_scores`
  turns this off). A read of a known strain then stays on its species instead of going to a congener
  that happens to fit the reference better. A read that one species fits clearly better keeps its
  scores. Five features, in the groups `alleles` and `polymorphic`, say how a taxon's reads stand
  against the known strains and at the sites where the species varies. *Why it matters:* these are
  the misses that remain. On the r226 v19 error reads, the sites where a species is fixed separate
  missed real strains from false calls (AUC 0.59-0.67), which the alignment alone does not. The GTDB
  build takes the alleles from half of each species' genomes, chosen by a hash of the accession, and
  simulates strains from the other half only. Otherwise every simulated real strain would be one of
  its species' alleles, a gain no real sample gives
  ([features.md](features.md#the-strain-alleles-alleles-2026-10-08),
  [polymorphic sites](features.md#the-polymorphic-sites-polymorphic-2026-10-09)).
- **The congeners' gene copies.** `--build` aligns each marker gene's copies within their genus and
  stores, per copy, its nearest and median distance to the congeners' copies and its nearest
  congener (`congener_gaps.tsv`). A divergent short read (identity below 0.99) whose seeds fit more
  taxa equally well than `--align_top` took now tries up to 7 more of them in its best alignment's
  genus (`--adaptive_candidates`; the `ZC` tag lists those left untried). This finds a strain's own
  species when its seeds ranked below its congeners'. Four `gaps` features and one `untried` feature
  come from these ([features.md](features.md#the-gene-copies-gaps-and-foreign-reads-and-the-untried-candidates-gaps-foreign-untried-2026-10-07)).
- **Ancestry sites** (new in 0.7.9) are now the consensus of the congeners: where three or more were
  compared (the 16 nearest and the gene's nearest by alignment), a site needs nine in ten of them on
  one base other than the species'; with fewer, the nearest congener decides, as in 0.7.9. The comparison continues across indels,
  and the indels the congeners share are sites too (two more features). At r226 the group added
  0.003-0.012 of AUC in the identity band where strains and novel congeners overlap (v18), and
  +0.0005 of F1 in v17's ablation. The refined sites did no better there than 0.7.9's
  ([features.md](features.md#which-side-the-reads-take-ancestry-079)).
- **85 features in the default set** (73 in 0.7.9): `gaps`, `untried`, `alleles` and `polymorphic`
  joined it, and `foreign` left it (below).
- **Speed, with the same outputs:**
  - Nanopore reads: an exact bound on the indel distance between read and gene window (a
    bit-parallel LCS) refuses a candidate before WFA2 when no alignment within the score budget can
    exist. The k-mer screen refused almost no Nanopore candidate at the 0.85 identity floor, and
    Nanopore reads were 70% of the r226 v19 build's aligning. On 2,000 benchmark Nanopore reads the
    bound refuses 18.5% of the candidates and cuts the alignment's instructions by 45%.
    `--no_indel_bound` brings back the k-mer screen; PacBio and short reads keep the screen
    ([report](claude/2026-10-09-ont-wfa2-skipping.md)).
  - AVX-512 (Ice Lake, Zen 4) for the seed lookup's flex scan and the syncmer scan, chosen at run
    time (`PROTAL_SIMD=auto|avx512|avx2|scalar`). On a Zen 4 node at r226, seeding took 9.5% less
    time for paired-end reads and 10.2% less for HiFi, the whole run 3.7% and 3.9% less
    ([report](claude/2026-10-08-avx512-kernels/README.md)).
  - BGZF read files (from `bgzip`, for example) are inflated on 1-4 threads per file: 1.1 GB/s on
    one thread, 2.7-3.4 GB/s on three. After the faster seeding, a paired-end run at r226 waits for
    its gzip input again. Plain gzip is still read on one thread
    ([report](claude/2026-10-08-parallel-bgzf/README.md)).
- **Also:** BGZF output (`.sam.gz`, the simulator's `.fq.gz`) has the same bytes whichever thread
  writes it (ISA-L hashed a byte from the address of its state;
  [report](claude/2026-10-08-bgzf-thread-bytes/README.md)). `misc/cpu.tsv` has a row per stage of a
  run, with its wall and CPU time. Fixed: a `--map` run with a relative output folder (`-o out`, or
  `#OUTPUT_DIR out` in the map) wrote into `out/out/`.

### For anyone building a database

- **The foreign leak.** The r226 v17 models' jump in F1 (paired-end test 0.983) came from the
  `foreign` features. Their scan read the genomes the training samples are drawn from, so it told
  the models which species the simulation can draw. Refit without them, v17 scores 0.962, the honest
  level of v14, v15 and v18. The scan now tiles each database's full reference,
  every species alike, and runs only with `--foreign-rates`. The group is in no feature set
  ([features.md](features.md#the-foreign-features-leak)).
- **Species complexes held out whole.** Congeners within 0.01 of each other on their marker genes
  (`species_clouds.tsv`) are held out or kept together. The test of novel species no longer has
  near-identical twins on both sides of the split (2% of v17's novel-congener rows). The training
  table's `meta_novel_distance` and the trainer's errors by distance show how far each held-out
  species is from the nearest kept one.
- **In-silico strains like real ones.** They now carry as many marker substitutions at their nearest
  congener's sites, with its base, as real strains do (`insilico_strains.py --clouds
  --congener-share auto`). Before, they were too easy to tell from congeners: the species' base was
  missing at 3.5% of the congener sites, against 10% for real strains (r226 v19).
- **Outputs and logs.** `--outdir` holds `protal_db/`, `model_logs/`, `logs/`, `work/` and
  `console.log`, each file once: about 118 files at r226 instead of about 11,900
  ([outputs](databases.md#outputs)). The build checks each simulated sample's composition against
  the truth (`composition_accuracy.py`, per species too) and logs CPU use per command, every 15 s,
  and per protal stage. `--share-logs` also packs the simulators' logs and the in-silico strains'
  table.
- **Many cores.** The training database's files are derived on all cores, the setup steps run as a
  graph, the genomes are read once into a cache, the reports start right after training, and
  `--profile-ahead` is on by default. The r226 v19 build took 3:49 h at 84 threads; with these
  changes the [build-parallelism report](claude/2026-10-09-build-parallelism/README.md) estimates
  2:20-2:40 h at 120 cores. That has not been measured yet.

### Compatibility

- 0.8.0 reads 0.7 databases. A database without `strain_alleles.tsv` or `congener_gaps.tsv` works,
  but the features that need them are fixed values. A model trained with 0.8.0's default set should
  run with a database built by 0.8.0 (`--build` writes `congener_gaps.tsv`, and with
  `--full_reference` also `strain_alleles.tsv`; `--add_tables` stores either in an existing
  database).
- A training table of an older protal lacks the new columns: train it with the groups it has
  ([features.md](features.md)).
- `build_gtdb_database.py` refuses a protal or `simulate_metagenomes` binary of another version.
  Rebuild both after updating a checkout.

## Changes since 0.6.0a

### Reads and input

- Single-end (`se`), PacBio (`pb`) and Nanopore (`ont`) reads besides paired-end, each profiled with
  a model of its own (`--read_type`, or `READ_TYPE` in a map file). Long reads are aligned per gene,
  reads over 65 kb in chunks; since 0.7.4 through every link of their chain (31-38% fewer alignment
  instructions, the same F1). A single read file's type is told from its first reads (0.7.3).
- Pipes, gzip, BGZF and FASTA input. Input is inflated outside the reader lock (libdeflate and a
  vendored zlib-ng; since 0.7.8 ISA-L, 1.6-1.8x as fast on single-member gzip). Unusable inputs fail
  their sample instead of passing silently.

### The database

- One file, `database.protal` (seekable zstd, the index in columns), checked against its reference;
  107 MB instead of 3.4 GB for the 765-species tuning world, built 3x faster. The old raw files still
  work and can be converted both ways ([databases.md](databases.md#the-files-of-a-database)).
- `--add_model` installs a read type's model into the database (several at once since 0.7.3).
- Index format 2 with correct unique flags, built in parallel and the same for any thread count.
  No k-mer with an ambiguous base enters the index.
- Per-database tables that the profiling uses: the genes' within-species conservation
  (`gene_conservation.tsv`, 0.7.1), how the genes differ between congeners (`gene_congeners.tsv`,
  0.7.2), gene neighbours from every genome with lines for a species' own gene order
  (`gene_neighbours.tsv`, 0.7.1, in the database file since 0.7.3), gene copies near-identical to
  another genus's (`suspect_copies.tsv`, 0.7.5), species priors from GTDB's clusters
  (`species_priors.tsv`, 0.7.5; genome sizes since 0.8.0), each copy's distance to its congeners'
  copies (`congener_gaps.tsv`, 0.8.0) and alleles of the species' other genomes (`strain_alleles.tsv`,
  0.8.0; up to 16 per gene since 0.8.1), the columns' conservation in each family (`column_weights.tsv`,
  0.8.1), every copy's ancestry sites (`ancestry_sites.bin`, 0.8.1) and each species' NCBI species
  taxid (`ncbi_taxids.tsv`, 0.8.1).
- In memory the reference genes are held at two bits per base (0.7.1; 19.5 GB → 5.1 GB at GTDB
  r226 size) and the index's values packed at 42 bits (0.7.5; 35 → 27 GB at r226). The files are
  unchanged by either.
- `gene_table.bin` (0.7.6): `reference.map` and `unique_kmers.tsv` in binary, with the reference's
  fingerprint, as a member of `database.protal` beside the text tables. A run loads it on all threads
  without parsing (r226 size, six threads: 1.76 → 0.43 s). `--build` writes it; `--compress_db` adds
  it to an older single-file database; without it the text tables are read as before.
- Reduced databases (0.7.6): a subset of the marker genes (`--build_gene_subset`, the converter's
  `--genes`) now reaches every build phase. `scripts/rank_genes.py` picks the most distinctive genes
  with a share reserved for each domain, so that archaea keep genes
  ([databases.md](databases.md#building-a-database)).

### Profiles and strains

- Abundance from a species' own reads: `--depth_identity_margin` (default 0.08 since 0.7.1; 0.7.0's
  0.04 dropped 6-26% of the bases of strains 1-4.5% from the reference and biased abundances), a
  fragment's bases counted once where the mates overlap (0.7.3), foreign genes left out of depth
  and MSAs (`--keep_foreign_genes`, 0.7.3).
- Mates paired across neighbouring genes, and a mate rescued on the neighbour (0.7.1).
- Species calls: the presence model's features grew from the reads' counts and identities to the
  reads' other candidates (MAPQ, congener fits; 0.7.1), their divergence beyond what their base
  qualities expect and the conservation of the genes they hit (0.7.2), the sample's congeners
  (relatives and distance features, a read EM over the reads' alternatives; 0.7.4), the sample's
  depth and the evidence before the MAPQ filter (0.7.5). Species priors (`+priors`) are opt-in.
- Variants called per fragment with strand tests, deletions and unbiased depth; MSAs chosen by
  `--msa_knob`; strain rows from long reads with phasing (`--no_phasing` turns it off; 0.7.3).
  0.6.0a wrote 0-based partition files, which put every gene boundary one column off in IQ-TREE.
- Failures are reported and give a non-zero exit; outputs are written crash-safe.
- Outputs no longer depend on the order of the SAM's records (0.7.6), which multi-threaded alignment
  varies: sums of doubles are exact (fixed point), and ties break by taxid. Before, repeated runs of
  one build could differ in the last digit of a value or two
  ([order independence](claude/2026-10-05-order-independence/README.md)).
- The per-taxon files `misc/<taxon>.statistics.tsv` are written only with `--taxon_statistics`
  (0.7.6). They took 15 s of a 107 s r226 run on a network file system; the per-sample profile files
  hold the same numbers.
- In 0.8.0 (2026-10-08): the sample's composition. `species_priors.tsv` gives each species' genome
  size (GTDB's genomes, corrected by CheckM) and its marker genes' share of it; the aligner writes the
  reads it read into the SAM header; a profile ends with the share of the genomes its species do not
  explain (`?`), and its abundances are shares of all genomes (`--no_unknown_share`: of the called
  species, as before). `<profile>.composition` gives the share of the reads explained, the average
  genome size and its quantiles, and how many species the rest would be
  ([running.md](running.md#what-the-called-species-explain-the-unknown-share)). Older databases get the
  sizes with `gtdb_to_protal_db.py --priors_only` and `--add_tables`.
- In 0.8.1: `--ncbi_taxids` ends every profile line with the species' NCBI species taxid and NCBI's
  name of it ([databases.md](databases.md#ncbi-taxids)).

### Speed

- Short reads aligned from their anchors' exact matches, with flanks of up to three mismatches
  aligned without WFA2 (0.7.3 to 0.7.5); branch-free and AVX2 syncmer scans from 2-bit codes;
  k-mer lookups prefetched within the read and one read ahead (0.7.4).
- SAMs compressed in the alignment threads (`.sam.zst` by default), parsed as views and
  decompressed on threads of their own (0.7.2). One sample is profiled on all threads (0.7.2) and
  the profiling pass itself got cheaper (coverage as a difference array, the read EM on indices;
  0.7.4 and 0.7.5).
- One binary for every CPU: 13 hot functions are also compiled for x86-64-v3 (AVX2) and chosen at
  run time (0.7.3). Huge pages for the index where the kernel allows them.
- The database is loaded and the index built in parallel; since 0.7.6 the start-up's parts load
  side by side (the index beside the genome preload, the small tables each on a thread;
  `--sequential_load` restores the old order) and the gene tables are added per genome on all
  threads.
- At GTDB scale (0.7.6): a k-mer screen refuses a candidate before WFA2 when read and gene window
  share too few k-mers for any alignment within the score budget (exact; `--no_alignment_screen`
  turns it off). A long read's window ends at its chain's last link, so most long reads align through
  their chain. The failed candidates of the `ZF` tag are counted in vectors by taxid, not in a hash
  map per chunk; at r226 the paired-end SAM holds 45M unmapped records naming 140,000 taxa.
- Every run prints per sample the reads, the candidates tried, refused and aligned, and the records
  written. Every stage, the start-up and the teardown are timed (0.7.6).
- In 0.8.0: Nanopore candidates refused by an exact indel-distance bound before WFA2 instead of the
  k-mer screen (−45% of the alignment's instructions on benchmark reads; `--no_indel_bound`), the seed
  lookup's flex scan and the syncmer scan in AVX-512 where the CPU has it (`PROTAL_SIMD`; seeding about
  10% faster on Zen 4), BGZF input inflated on 1-4 threads per file, and `misc/cpu.tsv` with each
  stage's wall and CPU time ([above](#what-080-adds-and-why-it-matters)).
- In 0.8.1: the ancestry sites stored by `--build` (`ancestry_sites.bin`) and decoded by a run
  instead of computed on each copy's first record: a sample's profiling −47% of its instructions,
  the same profiles ([above](#what-081-adds)).

### Training a database

- `build_gtdb_database.py` builds and trains a database from a GTDB release: held-out clades and
  species, other strains, an independent test set, one model per read type, resumable reruns
  ([databases.md](databases.md#building-a-database), [databases.md](databases.md#the-presence-model)).
  0.6.0a shipped one model trained on older databases.
- In 0.8.0 (2026-10-08): the hold-out keeps species complexes whole (congeners within 0.01 on
  the references, `species_clouds.tsv` from `protal --write_species_neighbours` before anything is
  built; the training table's `meta_novel_distance` and the report's errors by that distance), and
  the foreign scan tiles each database's full reference (every species alike, every copy listed)
  instead of the genomes the samples are drawn from, which had told the r226 v17 models the
  simulation's species ([features.md](features.md#the-foreign-features-leak)); after the r226 v18
  build (no leak, 0.001 of AUC) the build scans only with `--foreign-rates`; and the ancestry sites
  are the consensus of the congeners (nine in ten of three or more compared, the nearest congener
  alone with fewer) rather than the nearest congener's differences, the comparison is chained
  across indels instead of ending at the first, and the indels the congeners share against the
  species are sites too, two more features (80 in the default set)
  ([features.md](features.md#which-side-the-reads-take-ancestry-079)). The build's `--outdir` holds
  `protal_db/`, `model_logs/`, `logs/`, `work/` and `console.log`, each file once: ~118 files at
  r226 instead of ~11,900 (the in-silico strains on the samples' disk, one error-read table per read
  type; [outputs](databases.md#outputs), [report](claude/2026-10-08-build-outputs.md)). Strain
  alleles: `--build` keeps up to 4 alleles of each species' other genomes per gene
  (`strain_alleles.tsv`, edits of the representative's copy, each nearer it than the nearest
  congener's), a short read's candidates are scored with their species' best allele
  (`--no_allele_scores`: not), and two `alleles` features join the default set (82); the GTDB
  build takes the alleles from half of the genomes by a hash of the accession and simulates strains
  from the other half only, so that no simulated strain is its own allele
  ([features.md](features.md#the-strain-alleles-alleles-2026-10-08)).
- In 0.8.0 (2026-10-09): the species' polymorphic sites (where one of its alleles differs). In the alignment, a
  read's difference that any of the species' alleles has counts as a match, and one at a polymorphic site with another
  base counts as half a difference. These shifts decide only the reads protal is unsure about (a candidate of another
  species within 3 mismatches). Three `polymorphic` features join the default set (85): of the polymorphic sites a
  taxon's reads cover, the shares with a known allele's base and with another, and the species' base at the ancestry
  sites fixed within the species against all of them
  ([features.md](features.md#the-polymorphic-sites-polymorphic-2026-10-09)). The in-silico strains put as many of their
  marker substitutions at their nearest congener's sites, with its base, as the real strains have there
  (`species_clouds.tsv` is now made before them; [report](claude/2026-10-09-r226-v19/README.md), section 7).
- In 0.8.1: the column weights (`weights`, 91 default features) from GTDB's own alignments and tree,
  the ancestry features on the sites the species' known genomes share, the fixed-site agreement as a
  feature, alleles chosen for coverage, up to 16 of them; the finished database's alleles from every
  genome, the training database's from half; `--insilico-strains auto` and more real strains
  downloaded, so that real strains are about half of the simulated strains; the deep alleles' index
  entries and the weighted site shift opt-in ([above](#what-081-adds)).
- 0.7.3 to 0.7.5 made the training data more realistic (PacBio HiFi reads with qualities by length
  instead of pbsim3's quality-0 reads; a mixed design of lognormal sigma 1.3 and 2.0 communities;
  30% of the species held out; deeper design points) and the build faster (one simulation queue,
  deep samples in chunks, the training database on `--scratch`, both collections profiled in one
  protal run).
- `--version` names the commit a binary was built from, which the build script checks.
- 0.7.6 changed the training defaults and design:
  - The species priors are opt-in (`+priors`); their gain at r226 was the cluster size alone, a rule
    the simulation cannot test.
  - Every species with one genome gets an in-silico strain to be simulated from: a codon-aware
    mutated copy of its representative, as far from it as the real strains are from theirs
    (`scripts/insilico_strains.py`, `--insilico-strains`).
  - Training depths of 2,000, 50,000 and 200,000 read pairs are added; r226 v10's test errors sat
    between the old points.
  - Databases of several GTDB releases, full and reduced, in two phases
    (`download_gtdb_releases.py`, `build_gtdb_releases.py`).
  - Simulated worlds can evolve each gene at the real r226 gene's speed
    (`simulate_gtdb_release.py --gene_rates r226`).
  - `trace_relatives.py` follows reads by their contig, so the build's gene-conservation trace works
    on GTDB's genomes.
  ([r226 v10 evaluation](claude/2026-10-04-r226-v10-evaluation/README.md),
  [mini database genes](claude/2026-10-05-mini-database-genes/README.md),
  [trace](claude/2026-10-05-trace-relatives.md)).

## Benchmarks

All versions were run on one simulated world: 900 species (12% archaea, 135 in no database),
marker genes in operon-like clusters, genes of different conservation
([v0.7.1 benchmark](claude/2026-10-01-v071-benchmark/README.md)). Each 0.7 version used the
databases and models its own pipeline trained on that world; 0.6.0a used its own databases and the
model it ships. "full" means every sampled species is in the database, "missing" that the training
database lacks 290 species (349 for 0.7.5, whose holdout design changed, so its missing rows are
not like-for-like). 0.7.4 never had a database of its own; 0.7.5 holds its changes. Machine: WSL2,
Intel Core Ultra 7 258V, 6 threads, 23 GB. The rows for 0.7.0 to 0.7.2 are the
[v0.7.3 benchmark](claude/2026-10-02-v073-benchmark/README.md)'s, those for 0.6.0a, 0.7.3 and
0.7.5 the [v0.7.5 benchmark](claude/2026-10-03-v075-benchmark/README.md)'s; the reruns of older
versions reproduced their calls to the fourth digit.

0.7.6 was not rerun on this world. Its short-read changes leave the alignments and profiles
identical to 0.7.5's (checked on 500k and 5M pairs), apart from the last digit of the few values
that depended on the SAM's record order. Long reads differ where a window's end moved (on the test
reads 51 → 56 of 60 now align through their chain; one HiFi read of the r226 sample). Its training
changes show in the GTDB models below.

### Detection and abundance, short reads

Means over 26 paired-end samples (1,000 to 5M pairs) and 24 single-end samples; F1 of the species
reported at knob 0.5, false positives per sample, Bray-Curtis dissimilarity of the abundances to
the truth (0 is identical), recall over all present archaea:

| reads | database | version | F1 | precision | recall | FP per sample | Bray-Curtis | archaea recall |
|---|---|---|---|---|---|---|---|---|
| pe | full | 0.6.0a | 0.634 | 0.993 | 0.556 | 1.00 | 0.345 | 0.242 |
| pe | full | 0.7.0 | 0.903 | 0.976 | 0.869 | 2.96 | 0.131 | 0.736 |
| pe | full | 0.7.1 | 0.908 | 0.980 | 0.873 | 2.12 | 0.087 | 0.747 |
| pe | full | 0.7.2 | 0.911 | 0.982 | 0.875 | 1.88 | 0.085 | 0.749 |
| pe | full | 0.7.3 | 0.913 | 0.982 | 0.876 | 1.35 | 0.084 | 0.742 |
| pe | full | **0.7.5** | **0.915** | **0.988** | 0.876 | **1.27** | 0.085 | 0.745 |
| pe | missing | 0.6.0a | 0.629 | 0.982 | 0.560 | 1.50 | 0.357 | 0.225 |
| pe | missing | 0.7.0 | 0.895 | 0.963 | 0.865 | 2.00 | 0.131 | 0.697 |
| pe | missing | 0.7.1 | 0.900 | 0.967 | 0.869 | 1.54 | 0.093 | 0.704 |
| pe | missing | 0.7.2 | **0.906** | 0.977 | 0.869 | 1.08 | 0.088 | 0.700 |
| pe | missing | 0.7.3 | 0.902 | 0.963 | 0.872 | 1.46 | 0.089 | 0.708 |
| pe | missing | 0.7.5 (349 species out) | 0.902 | 0.966 | 0.872 | 1.77 | 0.098 | 0.706 |
| se | full | 0.7.0 | 0.874 | 0.975 | 0.830 | 2.88 | 0.145 | 0.682 |
| se | full | 0.7.1 | 0.888 | 0.981 | 0.844 | 1.92 | 0.107 | 0.696 |
| se | full | 0.7.2 | 0.888 | 0.982 | 0.843 | 1.75 | 0.107 | 0.700 |
| se | full | 0.7.3 | **0.893** | 0.984 | 0.847 | 1.08 | 0.106 | 0.708 |
| se | full | **0.7.5** | **0.893** | **0.987** | 0.847 | 1.08 | 0.106 | 0.702 |
| se | missing | 0.7.0 | 0.866 | 0.962 | 0.826 | 2.04 | 0.146 | 0.648 |
| se | missing | 0.7.1 | 0.881 | 0.973 | 0.836 | 1.42 | 0.113 | 0.648 |
| se | missing | 0.7.2 | 0.884 | 0.977 | 0.839 | 1.33 | 0.109 | 0.648 |
| se | missing | 0.7.3 | **0.886** | 0.978 | 0.838 | 1.21 | 0.111 | 0.652 |
| se | missing | 0.7.5 (349 species out) | 0.881 | 0.962 | 0.846 | 1.71 | 0.113 | 0.647 |

0.6.0a profiles paired-end reads only. The gain from 0.6.0a to 0.7.0 is the model's: 0.7.0 with
0.6.0a's model detects as 0.6.0a does, so it is what a database trained for its own reference gives
([0.7 against 0.6.0a](claude/2026-09-30-v07-vs-v06/README.md)). The step from 0.7.0 to 0.7.1 in
Bray-Curtis (0.131 → 0.087) is the depth margin 0.08. Per sample, 0.7.1 over 0.7.0 is +0.0055 F1
(significant), 0.7.2 over 0.7.1 +0.0022 (significant), 0.7.3 over 0.7.2 +0.0020 and 0.7.5 over
0.7.3 +0.0019 (both within noise); 0.7.3 over 0.6.0a is +0.279 (+0.185, +0.370).

By depth, F1 / Bray-Curtis with the full database (from the
[v0.7.1 benchmark](claude/2026-10-01-v071-benchmark/README.md); the later versions change these
rows by at most a few thousandths):

| read pairs | 0.6.0a | 0.7.0 | 0.7.1 |
|---|---|---|---|
| 1,000 | 0.242 / 0.588 | 0.749 / 0.192 | 0.755 / 0.155 |
| 10,000 | 0.602 / 0.444 | 0.949 / 0.094 | 0.954 / 0.063 |
| 500,000 | 0.964 / 0.040 | 0.984 / 0.076 | 0.984 / 0.028 |
| 5,000,000 | 0.986 / 0.043 | 0.980 / 0.090 | 0.992 / 0.030 |

Archaea in samples of 10,000 pairs: 0.6.0a found 3%, the 0.7 versions 82-83%; at 500,000 pairs 63%
against 98-100%. Species at 0.01-0.1% abundance in 10,000 pairs: 8% against 84-85%.

### Long reads

F1 (false positives per sample) on 4 PacBio HiFi and 4 Nanopore samples of 3 and 90 Mb each, made
by 0.7.3's collector (HiFi reads with qualities by length, which is what real HiFi reads carry;
0.7.0 to 0.7.2 were trained on pbsim3 reads with every base quality 0):

| reads | database | 0.7.0 | 0.7.1 | 0.7.2 | 0.7.3 | 0.7.5 |
|---|---|---|---|---|---|---|
| pb (HiFi) | full | 0.818 (3.00) | 0.824 (1.38) | 0.817 (2.75) | 0.826 (0.62) | 0.826 (**0.12**) |
| pb (HiFi) | missing | 0.792 (2.50) | 0.798 (1.75) | 0.793 (2.62) | 0.800 (0.88) | **0.811** (**0.12**) |
| ont | full | 0.865 (1.25) | 0.870 (0.88) | 0.872 (0.38) | 0.872 (0.62) | 0.873 (0.38) |
| ont | missing | 0.843 (2.38) | 0.854 (1.75) | 0.845 (2.50) | 0.858 (1.12) | 0.860 (0.75) |

0.7.5 over 0.7.3 on PacBio with species missing is +0.0113 per sample (significant), with
0.75 fewer false positives and closer abundance. The long-read F1 is lower than the short-read one
because half the samples are 3 Mb, where the depth decides.

### Strains

Twelve species of 8 strains each, one strain per sample at 2-40x, each version against its own
database; SNPs of the strain MSAs after qcmsa against the true genotypes, trees by IQ-TREE against
the true topology:

| reads | version | MSA columns (raw / qcmsa) | SNP precision | SNP recall | SNP F1 | miscalled SNPs per 100 kb | true topology |
|---|---|---|---|---|---|---|---|
| pe | 0.6.0a | 114,125 / 108,841 | 1.000 | 0.768 | 0.869 | 0.0 | 4 of 12 |
| pe | 0.7.0 | 114,640 / 114,416 | 0.998 | 0.937 | 0.967 | 2.0 | 10 of 12 |
| pe | 0.7.1 to 0.7.5 | 114,695 / 114,695 | 0.998 | 0.939 | 0.967 | 1.9 | 11 of 12 |
| pb | 0.7.0 to 0.7.2, 0.7.3 and 0.7.5 `--no_phasing` | 115,220 / 114,112 | 0.951-0.953 | 0.867-0.868 | 0.907-0.908 | 45.8-47.9 | 6-7 of 12 |
| pb | 0.7.5 (phasing on) | 115,171 / 114,064 | 0.948 | 0.861 | 0.903 | 50.9 | 8 of 12 |
| ont | 0.7.0 to 0.7.5 | 114,681 / 114,110 | 0.994 | 0.889 | 0.939 | 5.7 | 8 of 12 |

0.6.0a left many true SNPs as N at low depth. On the strain audit's world (8 species, 1-50x), its SNP
recall at covered true SNPs was 0.57 at 1x and 0.65 at 2x against 0.993 and 0.996 for 0.7.0, and its
qcmsa removed more cells the deeper the sample (0.78 of the positions called at 50x against 0.96).
The PacBio samples hold one strain each, so 0.7.3's and 0.7.5's phasing, which split 8 of 96
sample rows, is wrong there and is scored without those rows.

### Speed

Mean wall / CPU seconds per run at `-t 6`, both databases, no qcmsa, from the v0.7.5 benchmark
(0.7.0 to 0.7.2 from the v0.7.3 benchmark, whose deep runs shared the machine with a build, so
their 5M-pair and long-read rows are inflated; 0.6.0a, 0.7.3 and 0.7.5 ran in turn on an otherwise
idle machine). Repeated runs vary by 10-20%.

| reads | depth | 0.6.0a | 0.7.0 | 0.7.1 | 0.7.2 | 0.7.3 | 0.7.5 |
|---|---|---|---|---|---|---|---|
| pe | 1,000 pairs | 3.1 / 3 | 1.1 / 4 | 1.1 / 4 | 1.1 / 4 | 0.9 / 3 | 1.0 / 3 |
| pe | 10,000 pairs | 3.9 / 4 | 1.1 / 4 | 1.1 / 4 | 1.0 / 4 | 1.3 / 4 | 1.5 / 4 |
| pe | 500,000 pairs | 13.2 / 36 | 7.6 / 32 | 6.5 / 25 | 5.1 / 25 | 6.1 / 28 | 5.7 / 25 |
| pe | 5M pairs | 59.8 / 276 | 76.0 / 339 | 81.0 / 311 | 58.4 / 307 | 30.5 / 167 | 29.7 / 158 |
| se | 500,000 reads | - | 4.0 / 16 | 3.7 / 14 | 2.9 / 12 | 3.0 / 13 | 3.3 / 13 |
| pb | 90 Mb | - | 10.1 / 54 | 11.1 / 60 | 10.6 / 59 | 6.0 / 32 | 4.7 / 25 |
| ont | 90 Mb | - | 10.4 / 54 | 10.4 / 54 | 10.0 / 55 | 9.1 / 49 | 7.2 / 39 |

Like-for-like comparisons on a quiet machine, from the reports that made them:

| comparison | data | 0.6.0a or older | newer | source |
|---|---|---|---|---|
| 0.6.0a → 0.7.0 | 5M pairs 2x150, 150 species, tuning world | 100 s (310 CPU s) | 49.1 s (187 CPU s) | [0.7 against 0.6.0a](claude/2026-09-30-v07-vs-v06/README.md) |
| 0.6.0a → 0.7.0 | 500k pairs 2x150 | 10.2 s | 6.0 s | same |
| 0.6.0a → 0.7.0 | strain run A, 42 samples with qcmsa | 116 s | 36.6 s | same |
| 0.6.0a → 0.7.0 | index load, 5M-pair run | 9.4 s | 1.7 s | same |
| 0.6.0a → 0.7.5 | 5M pairs, operon world | 59.8 s | 29.7 s | v0.7.5 benchmark |
| 0.6.0a → 0.7.5 | 500k pairs, operon world | 13.2 s | 5.7 s | v0.7.5 benchmark |
| 0.7.3 → 0.7.5 | 90 Mb Nanopore / PacBio HiFi | 9.1 / 6.0 s | 7.2 / 4.7 s | v0.7.5 benchmark |

At GTDB scale: 0.7.5 against 0.7.6 on the r226 database (143,614 species, 27 GB file) with real reads
(a gut paired-end sample, a PacBio HiFi barcode; the read counts were not recorded), AMD EPYC 9634 at
32 threads, warm runs, medians of three
([performance at GTDB scale](claude/2026-10-04-performance-gtdb-scale/README.md)). The 0.7.6 column
is the fourth cluster run (`8c37dd6`, wall time; the stage times from the third, `7693f28`), before
`gene_table.bin` and the concurrent start-up:

| | paired-end 0.7.5 | paired-end 0.7.6 | HiFi 0.7.5 | HiFi 0.7.6 |
|---|---|---|---|---|
| wall | 131.5 s | **58.5 s** (−55%) | 103.7 s | **24.6 s** (−76%) |
| aligning | 61.9 s | 38.3 s (the k-mer screen) | 83.7 s | 11.8 s (screen and chain end) |
| profiling | 36.0 s | 6.4 s (failed candidates by taxid) | 3.4 s | 0.5 s |
| per-taxon statistics files | 14.9 s | 0 (`--taxon_statistics`) | 1.6 s | 0 |
| gene tables | (untimed) | 3.4 s (6.0 s before the parallel load) | (untimed) | 3.4 s |
| untimed | 22-26 s | under 1 s | 9 s | under 1 s |

The remaining fixed cost per run is about 11 s (13.8 s in the third run): the gene tables, the
genome preload (3.2 s), the index (3.2 s warm, 131 s cold from NFS) and the small tables (1.4 s).
`gene_table.bin` and the concurrent start-up should take it lower; that run has not been made. Seeding (20.6 s per thread) is what is
left of the paired-end alignment.

Where the time went, by performance round (callgrind instruction counts do not depend on the
machine's load; wall times are one thread unless said):

- [Round 1](claude/2026-09-29-performance-profiling/README.md) (before 0.7.0): gzip was inflated
  inside the reader lock (1M pairs at 8 threads: 21 s gzipped, 8.5 s plain), the syncmer test was
  recomputed per window (40% of the loop on a realistic mix), the 3 GB key map was read through
  4 KB pages (huge pages: seeding 5.9 → 3.7 s), the model's parser threw 229k exceptions at start.
  Fixed in 0.7.0: threaded gzip reader, branch-free and AVX2 syncmers, anchored alignment, zstd SAM
  output, zlib-ng, parallel index passes, huge pages, cPMML patches.
- [Round 2](claude/2026-09-30-performance-round2/README.md) (just before 0.7.1, against a build of
  the day before 0.7.0): aligning 1M realistic pairs at one thread 18.0 → 9.5 s, at six 5.6 → 3.2 s;
  then three cheap fixes (no per-read clock
  reads, a table reverse complement, no per-anchor read copy) took another 19% of the wall time, the
  windowed gene decode and syncmers from 2-bit codes a few percent more.
- [Round 3](claude/2026-10-02-performance-round3/README.md) (0.7.3 → 0.7.4): flanks of one
  mismatch without WFA2 and a faster index decoder (−7.7% of a 100k-pair run's instructions);
  long reads through their chain (−31% to −38% of the long-read alignment's instructions, 15-21%
  less CPU on 90 Mb samples). [Prefetching](claude/2026-10-02-prefetch/README.md) the k-mer
  lookups: about 4-5% of alignment. The [one binary](claude/2026-10-02-one-binary/README.md):
  −3.5% instructions on short reads against plain x86-64, and the launcher is gone.
- [Review](claude/2026-10-03-performance-review/README.md) and
  [pass 2](claude/2026-10-03-performance-pass2/README.md) (0.7.4 → 0.7.5): the profiling stage had
  grown to 30-40% of a short-read run with the congener features (5M pairs at six threads 34.0 s,
  of which profiling 12.6 s, the read EM 7.5 s on one thread); the EM on indices and cheaper
  congener sketches took it back (profiling −65% at one thread), then coverage as a difference
  array and 2-3 mismatch flanks without WFA2 (profiling −16.7% of instructions, aligning −2.4%).
  At 0.7.5 a 500k-pair run at one thread is 13.0 s: aligning 9.4 s, profiling 2.0 s, loading 1.9 s.

The [multithreading audit](claude/2026-10-01-multithreading-audit/README.md) found the locks
cheap at 0.7.1 (the reader lock 1.3-1.5% of thread time at six threads, the SAM writer nothing);
what did not scale was one sample profiled on one thread, fixed in 0.7.2.

### Memory

Peak RSS of a run on the operon world (765-species database with gene neighbours, 6 threads):

| version | peak RSS | source |
|---|---|---|
| 0.6.0a | 3.21-3.39 GB | v0.7.1, v0.7.2 and v0.7.5 benchmarks |
| 0.7.0 | 3.48-3.56 GB | same |
| 0.7.1 | 3.44-3.61 GB | same |
| 0.7.2 | 3.44-3.49 GB | same |
| 0.7.3 | 3.45-3.48 GB | v0.7.5 benchmark |
| 0.7.5 | 3.42-3.45 GB | v0.7.5 benchmark |

At this size the index's fixed 3.2 GB key map is most of the memory, so the versions hardly differ
(0.7.0 takes about 0.15 GB more than 0.6.0a, the later versions a little less than 0.7.0). The changes
that matter are at GTDB scale, where the index values and the genes dominate. Breakdown of a full
GTDB r226 run (16.6M genes, 2.90 G index entries), estimated from measured per-value and per-gene
costs except where said ([memory profiling](claude/2026-09-30-memory-profiling/README.md),
[the 2-bit gene store](claude/2026-09-30-gene-store/README.md),
[memory audit](claude/2026-10-03-memory-audit/README.md)):

| part | 0.6.0a and 0.7.0 | 0.7.1 to 0.7.4 (2-bit genes, `Gene` 112 → 32 bytes) | 0.7.5 (index values packed at 42 bits) |
|---|---|---|---|
| index values | 35 GB (measured on the r226 v7 index: 35.0) | 35.0 | **27.3** (computed from the measured slot count) |
| key map | 3.2 | 3.2 | 3.2 |
| reference genes | 17 (1 byte per base) | **4.3** | 4.3 |
| gene tables | 1.9 | **~0.8** | ~0.8 |
| **run, about** | **~58-59 GB** | **~43 GB** | **~35 GB** |

The 59 GB is the figure the website's download page still gives for 0.6. Measured on r226 runs of
0.7.5 (the four cluster runs above, both read types): **38.0 GB** peak RSS, of which the index
27.3 + 3.2 GB as loaded and the genes, tables and run buffers 7.5 GB. 0.7.6 has not been measured
there; its changes add no table in memory. On the operon world the packed index saved 70 MB of
213 MB of values with byte-identical outputs. A database build at r226 holds the index and the reference
(peak 88 GiB with both the training and the final index built at once in the first r226 run; about
50 GB with `--one-build-at-a-time`). Reads, threads and depth add little (≤ 20 MB per thread;
profiling ~0.1 GB per million alignments, kept per sample only in `--map` runs).

### Databases: size and build time

The tuning world (765 species, `reference.fna` 95 MB), built from the same converted files at
6 threads ([0.7 against 0.6.0a](claude/2026-09-30-v07-vs-v06/README.md)):

| build | wall | CPU | peak RSS | on disk |
|---|---|---|---|---|
| 0.6.0a, raw files | 40.3 s (31.0 s at 1 thread) | 78.5 s | 3.47 GB | 3.82 GB |
| 0.7, raw files (`--no_compress`) | 13.1 s | 27.9 s | 3.65 GB | 3.80 GB |
| 0.7, `database.protal` at zstd level 3 | 4.9 s | 21.1 s | 4.35 GB | 129 MB + full reference |
| 0.7, `database.protal` at level 19 (the default) | 101.2 s | 180.5 s | 4.88 GB | 107 MB + full reference |

0.6.0a was slower on 6 threads than on one (its uniqueness check got slower with threads). The 0.7
build is bound by writing the index when it writes raw files, and by zstd at level 19.

GTDB r226 (143,614 species) with `build_gtdb_database.py`, on the HPC
([r226 build evaluation](claude/2026-10-02-r226-build-evaluation/README.md),
[build profiling](claude/2026-10-03-build-profiling-r226/README.md),
[r226 v9](claude/2026-10-03-r226-v9-evaluation/README.md)):

| run | code | threads | whole workflow | notes |
|---|---|---|---|---|
| v1 | `ff57266` (in 0.7.1) | 16 | 8:54 h | conversion 2:46, training database 0:40, final database (level 19) 1:20 in the background, training-data simulations 1:52 (paired-end) + 1:23 (long reads), protal on 480 samples 0:14, test set 1:01, training 6 min, `--add_model` x4 22 min; peak 88 GiB; `database.protal` 20.45 GB (its index 18.05 GB) |
| v5 | `27423c6` (0.7.4) | 64 | 1:55 h after the holdout | the collector waited ~40 min for the long-read simulations, which ran last and alone |
| v8 | `a3e397d` | 52 | 2:17 h from the end of the conversion | |
| v9 | `26b065c` (0.7.5's defaults) | 52 | 2:02 h from the end of the conversion | with 30% of the species held out, a 30M-pair point and the mixed design; index build 33.5 min, of which the suspect-copy scan 11.7 |

Between v1 and v5 the conversion no longer pipes 86 GB through one zstd process, the gene
neighbours take 9.7x less CPU, the simulations run on the node's scratch, and the collector profiles
all samples of a collection in one protal run. The 0.7.5 simulation queue (long reads largest first
and in chunks, beside the paired-end points) halved a six-core test collection's wall time
(6:05 → 3:11) with the same reads. 0.6.0a had no build-and-train pipeline, and its index build was
never run at r226 size; from the tuning-world builds it would take about three times 0.7.0's.

### The models at GTDB scale

The r226 runs trained a model per read type and scored it on an independent test set of the same
run. The test sets differ between runs (species held out 20% → 30%, abundances sigma 1.3 → 2.0,
deeper points), so the F1 columns are not strictly comparable; the trend is what they show. F1 of
the independent test set at knob 0.5 (false positives per sample where reported):

| run | code | pe | se | pb | ont | source |
|---|---|---|---|---|---|---|
| v1 | `ff57266` (0.7.1 era) | 0.948 (4.0) | 0.947 (3.9) | 0.886 (4.2) | 0.904 (5.9) | [r226 build evaluation](claude/2026-10-02-r226-build-evaluation/README.md) |
| v5 | `27423c6` (0.7.4) | 0.9615 | 0.9570 | 0.9714 | 0.9704 | [v5/v6](claude/2026-10-03-r226-v5-v6-training/README.md) |
| v7 | `a3e397d`, divergence features, knob curve | 0.9654 (2.50) | 0.9606 (2.91) | 0.9764 (0.93) | 0.9737 (1.00) | [v7/v8](claude/2026-10-03-r226-v7-v8-evaluation/README.md) |
| v8 | `a3e397d`, + the sample's depth | 0.9671 (1.64) | 0.9685 (1.46) | 0.9763 (0.64) | 0.9730 (0.78) | same |
| v9 | `26b065c` (0.7.5 defaults, `+priors`), harder test set | **0.9699** | 0.9668 | 0.9727 | 0.9741 | [v9](claude/2026-10-03-r226-v9-evaluation/README.md) |
| v10 | `281a4ba` (after 0.7.5: the default set without priors, as 0.7.6 trains), same samples as v9 | 0.9630 (2.15) | 0.9596 (2.26) | 0.9675 (0.90) | 0.9689 (1.10) | [v10](claude/2026-10-04-r226-v10-evaluation/README.md); [features.md](features.md) |

From v1 to v9 the paired-end false positives per test sample fell from 4.0 to about 1.6 at the same
or higher sensitivity; v1 was precision-limited, with false positives growing with depth (0.2 at
1,000 pairs to 8.6 at 500,000). The v9 ablation attributes +0.007 to +0.009 of paired-end and
single-end F1 to the species priors (mostly their cluster-size singleton rule, which still wants a
check on real data), which is why they are opt-in since `281a4ba` (0.7.6). v10 is the same samples
without them: 0.007 lower on the paired-end test set, all of it absent one-genome species called and
strains of multi-genome species missed. v11, the first build with 0.7.6's training design (in-silico
strains of one-genome species, training depths of 2,000, 50,000 and 200,000 pairs), will show
whether the priors' gain was the simulation's. The PacBio and Nanopore rows between v1 and v5 also
reflect the change from quality-0 to HiFi training reads and the knob curves.

## What these benchmarks do not cover

Accuracy on real samples (real samples were only timed, at r226); a 0.6.0a run at GTDB scale (its
memory and index build there are estimates); 0.7.6 on the benchmark world (by its checks the same
calls as 0.7.5 for short reads) and its memory at r226; strain mixtures and long-read strains of
HiFi reads; 0.7.4 on its own. 0.7.7 to 0.8.2 were not run on the benchmark world either: their
accuracy shows only in the r226 builds' test sets ([r226 v18](claude/2026-10-08-r226-v18/README.md),
[r226 v19](claude/2026-10-09-r226-v19/README.md), [r226 v21](claude/2026-10-09-r226-v21/README.md),
[r226 v22](claude/2026-10-10-r226-v22/README.md), [r226 v23](claude/2026-10-11-r226-v23/README.md),
and for 0.8.2, with Singletrack and the weighted site shift, [r226 v24b](claude/2026-10-12-r226-v24b/README.md)).
Every number above is from one laptop (6 threads) or one HPC node, on simulated worlds whose strains
are 0.4-4% from their references at the markers.

## Sources

The version benchmarks: [0.7 against 0.6.0a](claude/2026-09-30-v07-vs-v06/README.md),
[0.6.0a, 0.7.0 and 0.7.1](claude/2026-10-01-v071-benchmark/README.md),
[0.7.2](claude/2026-10-02-v072-benchmark/README.md), [0.7.3](claude/2026-10-02-v073-benchmark/README.md),
[0.7.5](claude/2026-10-03-v075-benchmark/README.md); for 0.7.6 the
[performance at GTDB scale](claude/2026-10-04-performance-gtdb-scale/README.md) and the
[r226 v10 evaluation](claude/2026-10-04-r226-v10-evaluation/README.md). Everything else is linked where it is used; the
full list is in [`docs/claude/README.md`](claude/README.md). The condensed changelog is the header
comment of `src/main.cpp`.
