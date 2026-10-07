# The database build: logic and speed

Date: 2026-10-07. Branch `audit-fixes` at `c9b46f9`. Read-only: nothing was changed, compiled or run (the local
cores are reserved for other work).

**Since implemented** (2026-10-07, the commit after this report's): L1, L2, L3, L5, L6, the log's core report (L16),
L8 as documentation (the table's distance is over all genes; a build no longer reads an old `unique_kmers.tsv`), and
S1-S4 with every output byte-identical (S3 also makes gene conservation deterministic for a full reference in frames).
Their speed was not measured: the laptop was busy, so the gain at r226 is for the next cluster build to show (its log
now counts the suspect scan's candidates merged in full and the uniqueness check's distinct k-mers and cells).

**Question.** Audit the database build for (a) logic and (b) speed.

**Scope.** Everything that turns a GTDB release into `database.protal`:
- the conversion (`scripts/mini_db/gtdb_to_protal_db.py`, `gene_neighbours.py`);
- `protal --build` (`src/Build.h`, `src/Hash/Seedmap.h`, `IndexCodec.h`, `KmerLookup.h`, the gene tables in
  `GeneConservation.h`, `GeneIncongruence.h`, `SpeciesNeighbours.h`, writing and packing in `Zstd.h` and
  `Database.h`);
- how `scripts/build_gtdb_database.py` drives its two builds: the training database without the held-out species,
  and the finished database.

Not in scope: the read simulator (audited the same day, [report](2026-10-07-simulate-metagenomes-audit/README.md)),
the collector's simulation and profiling, and the trainer.

**Sources.**
- The code at `c9b46f9`, read in four disjoint parts (index passes and flags; index writing and packing; gene-only
  tables; conversion and orchestration). Each finding below that is marked CONFIRMED was checked against the code
  path end to end. PLAUSIBLE means the mechanism is in the code but its size or trigger was not measured.
- The logs of the r226 builds v10-v14 (`local/vNN/training_db_index.log`, `index_and_package.log`, `convert.log`,
  `training_db.log`; SLURM jobs on 64-thread nodes, e.g. v14 = 24005369). Their extraction commands are at the end.
- The whole-build timings in the [v14](2026-10-06-r226-v14/README.md) and [v15](2026-10-07-r226-v15/README.md)
  reports.

## Summary

**Logic.** The core of the build is sound:
- The index is byte-identical for any `-t`.
- Every k-mer is checked against the full reference.
- Each write is read back and compared.
- The training database is exactly the finished one without the held-out species: the same gene ids, taxids and
  order, and the held-out species removed from both references and from the gene neighbours.

The problems are at the edges:
- Three resume paths of `build_gtdb_database.py` silently build the wrong thing (L1-L3).
- The gene-conservation factors still depend on thread timing (L4). The logs show it: the same input gave 33,583,317
  to 33,583,322 copies compared.
- A read error in the references ends a pass early and the build still exits 0 (L5). `unique_kmers.tsv` is written
  unchecked (L6).
- `species_neighbours.tsv` measures distances on all genes, not on the hittable genes a run uses (L8).

The 3% of k-mer values the log reports as not stored is a deliberate rule (cores of ≥2,048 values are left out
whole), not a loss biased by reference order.

**Speed.** An r226 `protal --build` takes 26-36 min on 64 threads. Most of it goes to four things, and changes that
leave the outputs byte-identical could bring a build to roughly 10-12 min (an estimate; nothing was timed):

| Phase | v14 training / finished | v10-v14 range | What would cut it |
|---|---|---|---|
| suspect-copy scan | 9:24 / 9:45 | 7:32-10:05 | 95% of its 15.5e9 sketch comparisons are thrown away at once: an exact prefilter, ~−7 min (S1) |
| index write + `database.protal` | 14:41 / 8:27 | 8:20-20:01 | the index is written, read back, copied and read back again: write it into `database.protal` directly, ~−4 to −6 min (S2) |
| uniqueness check | 3:06 / 3:29 | 2:46-3:29 | each core's cells scanned once per range instead of once per k-mer, and parallel decompression: ~−2:30 (S3) |
| gene conservation | 1:29 / 1:43 | 1:29-1:54 | bound by one decompressing thread: ~−1:15 (S3) |
| passes 1 + 2 | 1:19 / 0:53 | 0:41-3:10 | read again from disk after the preload: ~−1 min where it was evicted (S4) |
| unique k-mer statistics | 2:41 / 1:16 | 0:35-2:41 | not the comparisons; needs a timer split first (S5) |
| the rest (gene congeners, gene table, loading, checks) | 0:48 / 0:48 | | |
| **whole build** | **33:28 / 26:21** | **26:21-36:34** | |

For the GTDB pipeline at its defaults, none of this is on the critical path: the training simulations are (v15: 4:05
of 4:42). The build's speed matters for standalone builds, the reduced and multi-release builds, CPU hours, and how
soon the training database frees the scratch disk. On the pipeline's critical path, about 5-8 min can come off by
starting the simulations before the training database's files are derived, plus two small orchestration fixes
(P1, P4, P5).

## Where a build's time goes

The phase times of both builds of every r226 run since v10. `T` is the training database (zstd level 3, written to
the node's SSD), `F` the finished one (level 9, written to `/hpc-home`). Both run at the same time, with `-t 64` each,
next to the simulations.

| | v10 T / F | v11 T / F | v12 T / F | v13 T / F | v14 T / F |
|---|---|---|---|---|---|
| gene conservation | 1:33 / 1:39 | 1:36 / 1:43 | 1:36 / 1:48 | 1:39 / 1:54 | 1:29 / 1:43 |
| gene congeners | 0:16 / 0:18 | 0:17 / 0:17 | 0:22 / 0:23 | 0:21 / 0:23 | 0:21 / 0:23 |
| suspect copies | 7:48 / 7:53 | 7:32 / 7:40 | 9:23 / 9:47 | 9:32 / 10:05 | 9:24 / 9:45 |
| pass 1 | 0:09 / 0:12 | 0:10 / 0:14 | 0:10 / 1:54 | 0:39 / 2:16 | 0:40 / 0:10 |
| pass 2 | 0:32 / 0:36 | 0:33 / 0:35 | 0:34 / 0:54 | 0:37 / 0:54 | 0:39 / 0:43 |
| uniqueness check | 3:03 / 3:10 | 3:12 / 3:18 | 2:55 / 2:47 | 2:47 / 2:46 | 3:06 / 3:29 |
| unique k-mer statistics | 0:35 / 0:45 | 0:37 / 0:45 | 1:08 / 0:54 | 1:14 / 0:54 | 2:41 / 1:16 |
| write index | 7:48 / 4:24 | 9:09 / 4:18 | 15:10 / 4:25 | 4:39 / 3:08 | 12:22 / 3:31 |
| write `database.protal` | 5:28 / 10:26 | 4:44 / 8:00 | 4:51 / 5:39 | 9:21 / 5:12 | 2:19 / 4:56 |
| **Run build** | **27:22 / 29:34** | **28:11 / 27:07** | **36:34 / 29:03** | **31:24 / 28:27** | **33:28 / 26:21** |

What the table says:
- **The suspect-copy scan is the largest compute phase in every build.** It grew from ~7.7 to ~9.5 min between v11
  and v12 with the same number of comparisons (15.9e9 and 15.8e9 in the training database). v12 is the first build
  that runs it before the index (`2baef20`); why it got slower was not looked into.
- **The I/O phases vary 2-4x between identical builds.** The training database writes its index at level 3 more
  slowly than the finished one at level 9 in every run. Index plus bundle stays near 14 min in four of five runs:
  - v13: 4:39 + 9:21;
  - v14: 12:22 + 2:19.

  That looks like unflushed writes of one phase paid for in the next, on a disk the simulations write their reads
  to. The CPU work of the writer is ~20-60 s. In v14 the training index write went on ~4.6 min after the finished
  build had ended, so CPU contention alone does not explain it.
- **Pass 1 varies from 9 s to 2:16 for the same 15 GB.** The preload read `reference.fna` in under 2 s at the start.
  After the ~11 min of gene tables the file is evidently no longer cached, and it is read again from disk; in the
  finished build that disk is the network file system.
- **The smaller training database finishes last in four of five runs**, although it has 7.6% fewer sequences
  (13,424,837 against 14,527,738) and compresses at level 3 against 9.

## Logic findings

| # | Finding | Where | Severity | Status |
|---|---|---|---|---|
| L1 | A reduced-marker run interrupted, then a full run in the same OUTDIR: the "full" database is built from the subset's genes | `build_gtdb_database.py:1807, 1874-1879, 2058-2061` | medium | CONFIRMED |
| L2 | A rerun with another `--seed` (or any change to the in-silico key) while the finished database is kept: in-silico strains from a uniform ANI 97-99.5 instead of the real strains' divergence | `build_gtdb_database.py:1871-1873, 1909-1921` | medium | CONFIRMED |
| L3 | Two OUTDIRs sharing one `--scratch`: a rerun reuses the other build's training database (other species held out) | `build_gtdb_database.py:2074-2082` | medium-low | CONFIRMED |
| L4 | Gene-conservation factors depend on thread timing | `GeneConservation.h:522-530`, `Build.h:947-959` | medium | CONFIRMED, visible in the logs |
| L5 | A read error in `reference.fna` or `full_reference.fna.zst` ends a pass early; the build exits 0 | `FastaBatches.h:57-64`, `FastxReader.cpp:74-76`, `Build.h:1316-1319, 1615-1674, 947-959` | medium | CONFIRMED path |
| L6 | `unique_kmers.tsv` written in place, its stream never checked | `Build.h:1691-1693` | medium-low | CONFIRMED |
| L7 | The training database's ~41 GB of growth on scratch is not reserved against the simulations' `--keep-free` | `build_gtdb_database.py:2074, 2169-2171` | medium | PLAUSIBLE |
| L8 | `species_neighbours.tsv` distances are over all genes; a run's `relative_distance` over the hittable genes only | `SampleContext.h:112-125`, `GenomeLoader.h:631-648`, `Build.h:1095-1162` | medium-low | CONFIRMED |
| L9 | The training database's full reference holds the simulated strains' own copies | `gtdb_to_protal_db.py:327-348, 588-600` | design question | PLAUSIBLE effect |
| L10 | A failing background finished-database build stops the whole run, hours of simulations included | `build_gtdb_database.py:303-308, 2181` | low-medium | CONFIRMED (design choice) |
| L11 | No memory or core check for two concurrent builds plus the simulations | `build_gtdb_database.py:2181-2193`, `collect_training_data.py:139-142` | low-medium | CONFIRMED |
| L12 | A `--gene-ranking` of another release is accepted if its gene ids exist | `build_gtdb_database.py:2033-2036` | low-medium | CONFIRMED |
| L13 | Crash leftovers: an old separate index next to new tables, `.partial` files, no fsync before sources are removed, silent failure to remove a stale index | `Build.h:131-183, 170-176, 405-408`; `Zstd.h:1222-1231`; `Database.h:553-558` | low-medium | CONFIRMED paths |
| L14 | Weak input checks: the taxonomy is never validated in `--build`; a record-count mismatch is only a warning; the fingerprint covers `reference.fna`'s size, not its content | `RunProtal.h:168-215, 2950`; `ReferenceFingerprint.h:38` | low | CONFIRMED |
| L15 | Suspect-copy verdicts: congener search incomplete (bucket cap, one-sided probes), minima of noisy estimates; `--suspect_copy_distance` above 0.05 silently capped | `GeneIncongruence.h:209, 218-248`; `Options.h:177` | low | PLAUSIBLE / CONFIRMED |
| L16 | Small: IUPAC codes in `SameAsStored`, CRLF in `CheckGeneNeighbours`, latent caps and shifts, formal data races, misleading log lines | see below | low | CONFIRMED |

### L1. An interrupted reduced-marker run turns the next full run into a reduced one

`convert_key` (`build_gtdb_database.py:1807`) does not say whether a gene subset was asked for. A run with `--n-genes`
or `--genes`:
1. converts the release into `OUTDIR/.converted`;
2. marks the stage `convert`;
3. derives the subset's `reference.fna`, `reference.map` and `internal_taxonomy.dmp` into `protal_db` (2058-2061).

If it stops before the finished build has packed `protal_db` (an error, Ctrl-C, the SLURM limit), a later run
without a subset takes the branch at 1874: the stage is done and `protal_db` holds the three files, so it is "the
release converted by an earlier run".

What follows:
- The finished database is built from N genes.
- The training database is derived from it.
- `build_metadata.tsv` says `marker_genes all`.

Nothing warns.

**Fix:** record in the stage which folder holds the whole conversion. Or compare `protal_db/gene2geneid.tsv` with
`OUTDIR/gene2geneid.tsv` before resuming. Or forget `convert` when a subset folder is derived into `protal_db`.

### L2. A reseeded rerun makes its in-silico strains blind

When the finished database is kept (`final_done`), `converted` stays `None`. `positions_file()` then looks in:
- `protal_db`: `gene_positions.tsv` went into `database.protal` and was removed;
- `.converted`: removed after the derive.

So it finds none, and `insilico_strains.py` runs with `--ani 97-99.5` and without `--taxonomy`. Instead of each
species' strain being as far from its representative as the table's real strains are, gene by gene (the point of
0.7.6's in-silico strains), every one gets a uniform ANI.

Any change to `insilico_key` triggers it: `--seed` (the documented way to draw another holdout), `--insilico-strains`,
an edit of `insilico_strains.py`, or a lost `genomes_simulated.tsv`. The release is converted again anyway a few
steps later for the training database (2090), but after this step. The only sign is "(no gene_positions.tsv: ANI
97-99.5)" at the end of the step's line. `test_c_another_seed_then_sigterm` reruns with another seed but does not
look at the in-silico step.

**Fix:** do that reconversion before the in-silico step, or keep `gene_positions.tsv` (~250 MB) in OUTDIR. Or stop
with an error instead of falling back silently when `--insilico-ani` was not given.

### L3. A shared scratch folder mixes two builds' training databases

The training database lives in `SCRATCH/training_db`, but its stage key lives in `OUTDIR/.stages`. The only check
that the database is the one the key describes is `isfile(database.protal)` (2081-2082).

Scenario:
1. OUTDIR X builds its training database with holdout H1.
2. OUTDIR Y, with another seed or release, rebuilds `SCRATCH/training_db`.
3. X is rerun and reuses Y's database.

X's profiles are then consistent with Y's database, but X's `heldout_species.txt` names as novel species that the
database holds, and the other way round: the models learn wrong labels. `build_gtdb_releases.py` gives each database
its own scratch folder; runs started by hand are not protected.

**Fix:** keep the stage file inside `training_db`, or store the size and mtime of its `database.protal` in the key.

### L4. Gene-conservation factors depend on thread timing

`Estimator::Take` gives each (species, gene) the first `kMaxCopies` (16) copies that any thread reaches
(`fetch_add`, `GeneConservation.h:528`). Copies farther than 0.25 still use up a slot. The threads read 1 MB blocks
of the full reference at once, and a species' copies are scattered through each gene's region. So for every
species with more than 16 genomes (E. coli has thousands), which 16 are compared depends on timing and on `-t`.

The logs show it. The finished database's input was the same in v10-v14, yet it compared 33,583,321, 33,583,321,
33,583,321, 33,583,322 and 33,583,317 copies.

The factors feed:
- the conservation and divergence features;
- `--gene_conservation db`;
- `gene_congeners.tsv`;
- the in-silico strains (through the conservation factors).

A rebuild is therefore not reproducible bit for bit. An earlier note recorded this ("conservation factors
nondeterministic"); the 2026-10-05 order-independence fix covered ANI sums, not this.

**Fix (either):**
- Number each 1 MB block under the reader lock (block boundaries are deterministic), and keep each slot's 16
  smallest (block, record) ordinals.
- Or process each gene's frames in order on one thread, which goes with S3's parallel frames.

Related, low:
- **IUPAC codes.** `SameAsStored` (`GeneConservation.h:602-610`) maps every base other than C, G or T to A. The
  2-bit store maps Y, S and B to C, and K to G (`PackedSequence.h:11-14`). A representative gene with such a code is
  not recognised as its own copy, and is counted as another copy at d ≈ 0.0005.
- **The representative's own copy.** For species with more than 16 genomes it is often not among those taken, and
  then a real zero of an identical other genome is dropped instead.

### L5. A read error shortens a pass, and the build succeeds

Where:
- `FastaBatches::Read` (`FastaBatches.h:57-64`) takes a short read as the end of the file and never checks
  `m_is.bad()`.
- `BufferedFastxReader::LoadBlock` (`FastxReader.cpp:74-76`) returns false when the stream fails.
- A zstd error becomes `badbit` on the `std::istream`.
- The passes check only FASTA format errors (`Build.h:1316`). The uniqueness check (1615-1674) and gene conservation
  (947-959; `SeqReader::Success()` is never asked) check nothing.

So a truncated or corrupt `full_reference.fna.zst`, such as a disk-full write in the derive, prints one error line.
After that:
- The uniqueness check ends early: k-mers shared with the genomes after the cut stay "unique".
- Gene conservation gets fewer copies, or none, and says "no factors" instead of failing.
- `protal --build` exits 0, and the pipeline goes on.

A truncated `reference.fna` would be a consistent but incomplete index, and `CheckReferenceAgainstMap` only warns
when the record counts differ (L14).

**Fix:** after each pass, the uniqueness loop and the conservation loop, exit 8 if the stream is `bad()` or the
reader reports an error. Set `m_error` in `FastaBatches::Read` on `bad()`.

### L6. `unique_kmers.tsv` is written unchecked

`Build.h:1691-1693` opens the file in place, writes it and closes it. There is no `.partial` and no `fail()` check;
every other table writer exits 8 on a bad stream (968-974, 1056-1061, 1145-1150, 1223-1238).

On a full disk:
- A cut in the middle of a line is probably refused when `gene_table.bin` is made (9 columns expected).
- A cut at a line boundary passes. The genes after it get no unique-k-mer counts, so they are not hittable and runs
  ignore them. That table is packed and "verified" as it is.

**Fix:** write `.partial`, check the stream, rename.

### L7. The training database can fill the scratch disk under the simulations

The simulations keep `--keep-free` (30 GB) free on scratch and wait otherwise; profiling, which frees space, cannot
start before the training database exists. That database grows by ~41 GB in its last ~15 min (v14):
- the derived folder is ~19-20 GB;
- the index adds 18.3 GB;
- `database.protal` adds 21.8 GB before the separate files are removed.

The peak is ~61 GB. v14's simulation log already has 15 "waits for room" lines. (v12's "No space left on device"
came 2:42 into its simulations, long after the training database, before `--profile-blocks`; it is not this case.)

**Fix:** while the training database is built, give the simulations a reserve of `--keep-free` plus ~45 GB (or 2.5x
the derived folder), and drop it when the build ends. S2 also removes 18 GB of this peak.

### L8. `species_neighbours.tsv` is not the distance a run computes

`ReferenceSketch` takes the genome's hittable genes (`GetHittableGenes`). In a run those are the genes with
short_unique + long_unique > 0 (`GenomeLoader.h:1001`). The build writes the table before the index, when no gene
is known to be hittable yet, so it takes every gene (`GenomeLoader.h:631-636`).

Genes identical between two congeners have no unique k-mers, so a run leaves them out of both sketches, and the
build counts their zeros. For near-identical congeners, which are exactly the pairs `db_congeners_01` and `_02`
count, the build's median is smaller. The run's value can be about twice as large, or 1 (`kFarDistance`) if fewer
than 10 hittable genes are shared.

The models are trained and run on the same table, so training and use agree. But the table does not mean what
`docs/databases.md` and `Build.h:1087-1089` say ("the distance a run's `relative_distance` reads"), and features
that compare the two measures (e.g. `unexpected_congener_fit_share`) mix them.

**Fix:** write the table after the index (load the hittable set from the counts just made and reload the genes,
~10 s), or change the comments and the docs to "all genes".

### L9. Design question: the training database already knows the simulated strains

The derive removes held-out species by taxid only (`gtdb_to_protal_db.py:588-600`). The full reference keeps every
GTDB genome of the kept species, including the non-representative genomes the samples are simulated from. So a
k-mer of a simulated strain that another species' reference shares is already non-unique in the training database.
A real strain that GTDB lacks does not get that, so on cross-species false positives training is easier than real
samples. Only the in-silico strains are honest here.

The size of the effect is unknown.

**Option:** carry the accession in the full reference's headers (`>taxid_gene_ACC`; `ExtractHeaderInformation` stops
at the second underscore) and leave the simulated genomes out of the training database's full reference.

### L10-L12. Robustness of the orchestration

- **L10.** A failed background build of the finished database ends the run within ~5 s (`Job.ended`), killing hours
  of simulations. That database is needed only for `--add_model` at the end. Better: report the failure and rebuild
  it in the foreground after the trainer (the existing `--one-build-at-a-time` path), and fail only if that fails
  too.
- **L11.** Nothing reads `MemAvailable`, the cgroup's `memory.max` or the CPU affinity. During the builds there are
  64 build threads twice plus 64 simulation slots per collection, ~256 runnable threads on 64 cores. v14's two index
  writes overlapped at ~72 GB together. Better: below ~2 × 40 GB plus a margin, switch to `--one-build-at-a-time`
  with a note, and warn when `-t` exceeds the allowed CPUs.
- **L12.** Gene ids are ordinal per release (`gtdb_to_protal_db.py:688-691`). A ranking of r220 passed to an r226
  build picks other markers without an error. Better: also require that each row's marker name equals the release's
  name for that id.

Smaller resume-key points:
- `final_key["genes"]` hashes `gene_subset.txt` with its comment lines, which name the ranking's path.
- `file_identity` resolves a `protal` wrapper script, not the binary it calls.
- `release_identity` records only the marker folder's entry count and mtime.
- Duplicated markers within a genome: the first copy is taken separately for `reference.fna` and the full reference
  (moot at r226: none).

### L13. Crash leftovers

- **An old separate index during a rebuild.** A folder that still holds `index.prx(.zst)` (after `--unpack_db` or a
  `--no_bundle` build) and is rebuilt with the same reference has new `gene_conservation.tsv`,
  `suspect_copies.tsv`, `species_neighbours.tsv` and `unique_kmers.tsv` next to the old index for ~30 min, and for
  good if the build is killed. The fingerprint still matches, and separate files win over `database.protal`.
  **Fix:** remove or rename any old index and `gene_table.bin` in `--db` when the build starts, or write all tables
  as `.partial` and rename them after the index.
- **`.partial` files.** protal has no signal handler, so a kill leaves `index.prx.zst.partial` (18 GB) and
  `database.protal.partial` (22 GB). On scratch with `--keep-free` that is 40 GB nobody counts. **Fix:** remove the
  known names at the start of the build.
- **No fsync.** No `fsync` happens before the rename, and the sources are removed after a check that read the page
  cache: a node crash on a local file system can leave neither copy. Only `ReplaceTail` (`--add_model`) syncs.
  **Fix:** `fdatasync` the `.partial`, sync the directory, then remove the sources.
- **A stale raw index that cannot be removed** is ignored silently (`Build.h:176`). With `--no_bundle` it then
  shadows the new `index.prx.zst`.
- **`--compress_frame_mb 0` and `--no_compress`** write the index without reading it back, while the log says
  "verified".

### L14. Weak input checks in `--build`

- `internal_taxonomy.dmp` is read only by three lenient parsers (`Genera`, `LineagesFromTaxonomy`,
  `CheckGeneNeighbours`). Duplicate names, cycles, undefined parents or a reference taxid missing from it are packed
  and fail only at the first query (`LoadTaxonomy` runs when profiling, `RunProtal.h:2950`). **Fix:** load it with
  `IntTaxonomy` in `--build`, and check that every taxid of `reference.map` is a species in it.
- `CheckReferenceAgainstMap` only warns when the counts of records and map rows differ (`RunProtal.h:211-214`).
  **Fix:** make it an error.
- The reference fingerprint is the hash of `reference.map` plus the size of `reference.fna`. The index is built from
  `--reference`, the fingerprint and the bundle use the folder's `reference.fna`, so another `--reference` with the
  same headers and lengths goes unnoticed (the GTDB script passes the folder's own file).

### L15. Suspect copies

The rule is applied per gene to candidates found through 8 probe hashes, skipping buckets of more than 2,048 copies
(`GeneIncongruence.h:205-216`).
- **Incomplete congener search.** A copy whose near congener is reachable only through a skipped bucket, or only
  from the congener's side, gets `congener = 2`. It then counts as suspect as soon as any foreign copy is within 0.02.
- **Noise.** With 128 hashes the standard error of a distance is ~0.003 at 0.02, and `foreign` is a minimum over
  many such estimates, biased low.
- **Effect.** 0.04% of copies are flagged, so the effect on calls is small.
- **Fix.** Re-check the few thousand flagged copies exhaustively against their genus and with exact distances before
  writing them.
- **The option's range.** `--suspect_copy_distance` above 0.05 has no effect beyond 0.05, because pairs beyond
  `kReportDistance` are dropped first. It is not validated (`Options.h:177`).
- **Training against finished.** The two databases flag differently: 4,870 against 5,589 copies in v14, because
  held-out congeners are missing. Suffixed GTDB genera (Bacillus, Bacillus_A) count as different genera; a flag
  still needs the copy to be nearer the other genus by the 0.02 margin.

### L16. Small items

- **`CheckGeneNeighbours` and CRLF** (`Build.h:847`): `Table::Read` strips `\r` before testing for an empty line; the
  second pass does not, so a CRLF file with a blank line fails with "clade 0 is not in the taxonomy".
- **`GetExact`'s cap** (`KmerLookup.h:338`) empties the matches above its cap. The build's cap (2048) cannot be
  reached by a kept core today. With a lower cap, an over-full exact match would leave the k-mer unique. Do not cap
  in the build's lookup.
- **`IndexRangeBits`** (`Build.h:1324-1328`) goes negative for maps of exact k ≤ 4 (undefined shift). The build
  always uses k = 15.
- **Formal data races.** The uniqueness check and `CountUniqueKmers` change flag bits with atomic byte operations
  while other threads read the same bytes non-atomically (`Seedmap.h:479-496, 563-571`). This is undefined
  behaviour in C++ and ThreadSanitizer will report it. It is harmless on x86: only flag bits change, and no reader
  uses a flag another thread changes in that phase.
- **The log.** The 65,536-line histogram is 65,536 of the build log's 65,672 lines (7,099 rows are non-zero).
  Several figures are misleading:
  - "keys: 1.5 GB" counts cells (the key map is 3.0 GiB);
  - "values: 30.2 GB" is the file's layout, not the 24.7 GB in memory;
  - "Failed to meet demands? 22476" counts control blocks;
  - "Processed reads" adds pass 1's records before the subset filter, pass 2's after, and the full reference.

  Print the cores and values left out, a banded histogram, and per-pass counts.

### Not a bug: the 3% of values "not stored"

"Stored 96.99 of total values ... Failed to meet demands? 22476" looks like loss. It is the rule at
`Seedmap.h:1659`: a 15-base core with 2,048 values or more is left out of the index whole. The other two caps (16,384
cells per key, 65,536 per block) can never bind for a kept core.

The v14 histogram adds up exactly:
- 24,112 cores of ≥2,048 values hold 83,781,057 of the 2,774,572,586 k-mer occurrences (3.02%);
- the largest core has 36,641 values;
- the 2,690,791,529 kept values are the index's entries.

The rule is deterministic and does not depend on reference order. It keeps the uniqueness flags consistent: a
dropped core has no entries to clear or to hit.

It is not the query's `--max_key_ubiquity` 256, which caps the cells tied at the best score. Cores of 257-2,047
values (21% of the entries) are used by queries. Raising the 2,048 would add ~74% to the uniqueness check's
comparisons and ~0.76 GB.

## Speed findings

### In `protal --build` (per build, r226, 64 threads)

| # | Change | Expected gain | Outputs |
|---|---|---|---|
| S1 | Suspect-copy scan: an exact bound before each comparison, a merge that stops early, flat sketch arrays, no sort of the candidates; then one gene at a time with its copies on all threads | 9.5 → ~1.5-2 min; the gene phase's peak −12 GB | byte-identical |
| S2 | The index written straight into `database.protal` (not `index.prx.zst` first, then copied), verified by frame hashes | finished −4 to −6 min, training similar; −18 GB of scratch at the peak | `database.protal` byte-identical |
| S3 | Uniqueness check through `PartitionedPass` (each core's cells scanned once per range), and the full reference decompressed in parallel frames | uniqueness 3:06 → ~0:30-0:40; conservation 1:30 → ~0:20 | byte-identical (conservation deterministic with L4's fix) |
| S4 | Passes 1 and 2 from the preloaded genes instead of re-reading `reference.fna` | −0:30 to −2:00 where the file was evicted | byte-identical |
| S5 | Unique k-mer statistics: split the timer, then sort (value, row) pairs, format rows in parallel, first-touch the arrays on their threads | 0:36-2:41 → ~0:15-0:30 (uncertain) | byte-identical |
| S6 | The seekable writer without its batch barrier: a reader thread, compressors taking frames from a counter, one ordered writer | 10-60 s, more under contention | byte-identical |
| S7 | Smaller: check the flag before the atomic clear (6.45e9 locked operations for 0.70e9 non-unique entries), reuse the scan's sketches for species neighbours, parse the taxonomy once, dedupe identical full-reference records (count them first) | seconds to ~1 min | byte-identical |

**S1. The suspect-copy scan compares 15.5e9 pairs to keep ~2.4M.**
- **Candidates.** For each copy of a gene, `Scan` collects the copies sharing any of its 8 smallest sketch hashes
  (`GeneIncongruence.h:205-216`): ~1,160 candidates per copy.
- **Comparisons.** It computes the full 128-hash distance for each (`SketchDistance`, a branchy merge through a
  pointer to a separately allocated vector), ~2.3 core-µs per candidate.
- **What is kept.** Only distances ≤ 0.05 are kept (line 222). That is 1.19M cross-genus pairs, and at most ~2.4e9
  comparisons even if every congener pair were that close; probably more than 95% are discarded.
- **The bound.** A distance ≤ 0.05 needs at least 49 of the union's bottom 128 hashes shared (Jaccard ≥ 0.378). A
  1,024-bit signature per copy, built per gene (13 MB), gives an upper bound on the shared count by
  `popcount(a & b)`. The merge can also stop once the hashes still to come cannot reach 49.

Both are exact, so `suspect_copies.tsv` and `gene_incongruence.tsv` stay byte-identical and the change can be checked
by `cmp`.

**Load balance.** The scan runs in parallel over the 168 genes only. About 120 bacterial genes are heavy, so on 64
threads it runs in about two waves with a tail, and every sketch of the database is held at once (~13.5 GB peak). One
gene at a time, with its copies on all threads (each copy writes only its own slots), fixes both.

**S2. The index's bytes are written twice and read four times.** `SaveIndex` writes `index.prx.zst` (18 GB) and reads
it back. `db::Write` copies its frames into `database.protal.partial`, on one thread, and `db::Verify` reads both
files again (`Database.h:515-533, 441-464`). The file is deleted afterwards.

That is ~36 GB written and ~72 GB read for an 18 GB index. On `/hpc-home` at ~80-90 MB/s it accounts for ~3.5-4.5
of the finished build's 8.5 min of writing, and on the contended scratch disk for ~40% of the training build's I/O.

**Fix:**
- a `db::Source` whose frames come from the index writer (the frame count is known in advance, so the directory
  frame can still go first);
- verify it with `index_codec::Verify` on the member's slice of the seek table;
- free the index before the reference is compressed;
- verify compressed members by a hash of each frame made while compressing, instead of a second single-threaded
  decompression of the source (`Database.h:466-485`);
- keep `index.prx.zst` only for `--no_bundle`.

**S3. The uniqueness check scans ~192 cells per k-mer, and its input arrives on one thread.**
- **The lookups.** Each of 14.8e9 full-reference k-mers looks its core up and scans all of its cells. The
  size-weighted mean core size is 197, and each core is looked up ~5.3 times per value it holds. The predicted
  2.83e12 cells compared match the logged 2.84e12. Cores of more than 256 values take 72.6% of the comparisons.
- **The clears.** There are 6.45e9 clear events for 0.70e9 non-unique entries, so each is cleared ≥9 times.
- **Fix: grouped lookups.** Run the check through `PartitionedPass`, as the two index passes already are:
  1. items are (core, flex part, taxid);
  2. per key range and round, sort them and collapse them to (core, flex, min taxid, max taxid);
  3. scan each touched core once.

  The decision rule is unchanged and clears commute, so the flags come out the same. Because the full reference is
  grouped by gene, a core's lookups fall into 2-3 rounds: ~7e9 cells instead of 2.84e12.
- **The input.** The full reference is decompressed under the reader's one lock, on whichever thread holds it, and
  copied twice per block (`SeqReader.h:69-83`, `FastxReader.cpp:65-121`, `Zstd.h:220-320`). Gene conservation,
  which does little per record, runs at 0.75-0.97 GB/s in every build: one zstd thread. That is the floor of both
  phases.
- **Fix: parallel frames.** The converter writes one frame per marker file (173), but the derive's single `zstd -T8`
  writes one frame for the training database. Write record-aligned frames with protal's seek table (protal reads
  that format already) and decompress them on 8-16 workers (128 MB windows each because of `--long=27`).

**S4. Passes 1 and 2 read `reference.fna` again.** The genes are preloaded (2-bit) when the build starts. After the
gene tables the file has been evicted from the cache in some builds: pass 1 took 9 s to 2:16 for the same 15 GB.
Feeding the passes from the preloaded genes in `reference.map` order removes the I/O. Records with non-ACGT bytes
must be re-read by offset, since the 2-bit store keeps IUPAC codes as bases. `CheckReferenceAgainstMap`'s third read
of the file (≤17 s, outside the timers) can fold into the preload.

**S5. Unique k-mer statistics.** The same comparison count took 1:14 in v13 and 2:41 in v14 (101,622,076,238 both
times), and 0:36 in v10/v11. The comparisons are ~2-5 s of it. What remains is likely:
- three locked `fetch_add`s per entry on random rows;
- a single-threaded sort of 13.4M rows by a random key;
- 13.4M string keys in a `tsl::sparse_map`;
- iostream formatting.

Split the timer first.

### In the GTDB pipeline (`build_gtdb_database.py`)

| # | Change | Gain on the whole build | Status |
|---|---|---|---|
| P1 | Start the simulations before deriving the training database's files: `--simulate_only` needs only the genome table and `heldout_species.txt`, both ready at line ~1972 | 3-5 min (the derive's duration, estimated from the 2026-10-03 report) | order CONFIRMED |
| P2 | Start the finished database's build after the training database's, or at nice 10 with half the threads | training database ~11-16 min sooner; less CPU taken from the simulations for ~26 min; less of L7 and L11 | PLAUSIBLE |
| P3 | Disk contention on scratch: `ionice -c2 -n7` for the simulations, or the training database on another disk; confirm with timers in `SaveIndex` / `db::Write` (compute, write, verify) and iostat first | 5-10 min of the training build (off the critical path at the defaults) | pattern CONFIRMED, cause PLAUSIBLE |
| P4 | The genome-table step globs every genome root six times recursively (`build_gtdb_database.py:395-405`): one `os.walk` with a suffix filter | 1-2 min (the step took ~3 min in v11) | PLAUSIBLE |
| P5 | Run `gene_neighbours.py` while the converter writes the full reference: it reads only `reference.fna`, `reference.map`, the taxonomy and `genome2tiid.tsv`, written before | ~1:13 | CONFIRMED |

Not worth doing:
- **Sharing tables between the two databases.** Uniqueness, suspect copies and species neighbours must differ
  (held-out species), and filtering the finished database's tables would leak them into training. Only gene
  conservation and congeners (~2 min) could be shared.
- **Speeding up the conversion further.** The 48 s of spooling, 31 s of sorting and 51 s of writing to the network
  file system would give about 1 min.

On P3: fix 4 of the [2026-10-03 report](2026-10-03-build-profiling-r226/README.md) (the training database on scratch)
removed 6 min of network writes in v5. Since the simulations run beside the build on the same disk, that database's
write phases are 2-4x slower than the finished one's. It still beats the network file system, but the disk is shared.

## Checked and found correct

- **Passes.** Passes 1 and 2 extract the same k-mers per record: the same subset filter, the same handler and
  `DropAmbiguousKmers`. The serial and partitioned paths are equivalent. The index and `unique_kmers.tsv` are
  byte-identical for any `-t`: counting commutes, pointers are laid out per part and then offset, values are placed
  per range in batch order, and flag clears are atomic, commutative and idempotent.
- **Uniqueness.** The rule is per taxon, the same in the single-entry and multi-entry branches: same taxon, skip;
  another taxon with an equal whole k-mer, clear all. Two entries of one taxon with the same whole k-mer are both
  cleared. `IndexedKmer`'s canonical k-mer matches the handler's; ties cannot occur (15-base cores, odd length). The
  build's "exact" (FlexCell equality) is the query's best score of 16. Single-entry cores get no two-flag, and the
  query treats them alike.
- **Widths.** Taxids, gene ids and lengths are limited to 20 bits up front (`RunProtal.h:173-190`). `PackedLayout`
  fits every stored position. `PackValue` exits on overflow rather than truncating. `SlotBits32` reserves enough
  bits for every core. Block offsets, main keys and batch offsets fit their integer types.
- **Writing the index.** Packed values are unpacked exactly into the file's 8-byte layout. `Verify` compares every
  key-map cell, value cell, the header (fingerprint, features) and the layout. Copied frames are compared byte for
  byte, compressed members by their whole decompressed content. Frames carry zstd checksums. The writer's output is
  byte-identical for any thread count. The models stay last for `--add_model` in place, and the index is freed
  before packing.
- **Crash safety, apart from L13.** Every step writes `.partial` and renames, and the sources are removed only after
  a verified write. The index goes first, so `database.protal` takes over at once.
- **Gene-only tables.** `Genera()` reads the rank column and stops at the root, and its genus definition is the
  run's (`AncestorsOfRank("genus")`). Species-neighbour batching compares a genus larger than the batch whole.
  `Set` keeps the 16 nearest by (distance, taxid). `SketchedTaxonDistance` equals the run's pairwise distance, and
  the sketch prefix property holds. `Scan` is deterministic for any thread count. `CheckGeneNeighbours`' gene-id
  vector has no off-by-one, and `CheckGenePositions` checks what it says. `gene_congeners.tsv` is deterministic.
- **The training database.** It is the finished one minus the held-out species: the same gene ids, taxids and order,
  and `reference.map` offsets consistent. The held-out species leave the full reference too (71,424,841 of
  79,520,648 records in v14), so uniqueness, conservation, suspect copies, congeners and species neighbours never
  see them. `gene_neighbours.tsv` and `gene_positions.tsv` are counted again without their genomes. The taxonomy,
  priors and `genome2tiid.tsv` keep them, as intended.
- **Conversion and resume.** The conversion is deterministic for any `-t`. In a single OUTDIR with its own scratch,
  the resume keys cover the converter and gene-neighbour scripts, the release, the genome table, the held-out
  species, the protal binary and both zstd levels. Stages are marked only after success, and interrupted
  conversions or builds are redone after `clear_build_outputs`.

## Tests

`BuildIndexTest` (`tests/e2e/test_protal_e2e.py`) checks that the index and `unique_kmers.tsv` are identical at
`-t 4` and with `--serial_index_passes`. Nothing checks the other tables:
- `gene_conservation.tsv` is not deterministic (L4);
- `suspect_copies.tsv` and `species_neighbours.tsv` are deterministic but untested.

A build at `-t 1` and `-t 4` with `cmp` of every table would catch L4 and guard S1-S3.

The resume paths of L1-L3 are untested. `test_c_another_seed_then_sigterm` reruns with another seed but does not
look at the in-silico step.

## Documentation

`docs/databases.md`, "Step by step", needs three changes:
- It lists the index and the uniqueness check as step 1; the build makes the gene tables first.
- It says the species-neighbour distance is "as a run's `relative_distance`" (L8).
- "The index is the same for any `-t`" is true of the index but not of `gene_conservation.tsv` (L4).

The website is not affected: none of this changes how protal is run.

## What to do first

1. **Before the next r226 build:** L1 and L2 (a few lines of the script each, both silent), L5 and L6 (exit on read
   and write errors), L3 (the stage file inside `training_db`).
2. **L4,** deterministic conservation, with S3's parallel frames: the build becomes reproducible.
3. **S1,** the suspect-copy prefilter: the largest single gain, byte-identical, easy to verify.
4. **P1 and P5:** the only items that shorten the whole pipeline at its defaults.
5. **Timers** splitting the index write, the bundle and the statistics into compute, write and verify. The next
   r226 run then tells whether P3 is the disk; S2 follows.

## Reproducing the numbers

All from the logs in `local/` (no runs):

```bash
L=local
# phase times of both builds
for v in v10 v11 v12 v13 v14; do for f in training_db_index.log index_and_package.log; do
  echo "$v $f"; grep -E " took " $L/$v/$f | grep -E "^(Gene|Suspect|Pass|Value|Uniqueness|Unique|Write|Run build)"
done; done
# copies compared in gene conservation (L4), sketch pairs (S1), values stored
grep -h -oE "\([0-9]+ copies compared\)|[0-9]+ sketch pairs compared|Stored [0-9.]+ of total values" $L/v1?/index_and_package.log
# the core-size histogram: cores of >= 2048 values and their values
awk 'NR>=90 && NR<=65631 && NF==2 && $1 ~ /^[0-9]+$/ { if ($1 >= 2048) { c += $2; v += $1 * $2 } }
     END { print c, v }' $L/v14/training_db_index.log
```
