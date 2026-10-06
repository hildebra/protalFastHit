# Why `--add_model` took ~20 minutes, and replacing the models in place

2026-10-06, branch `audit-fixes` at `15006c2` plus the working-tree change described here (`src/Utilities/Database.h`,
`src/Utilities/Zstd.h`, `src/Build.h`, `src/Options.h`, tests, `docs/databases.md`). Measured in WSL (Ubuntu 24.04,
Core Ultra 7 258V, ext4 on the laptop's NVMe, 23 GB RAM, other sessions' load 2-4).

## The question

Storing the trained models in a GTDB `database.protal` (`protal --add_model pe.xml,se.xml,pb.xml,ont.xml --read_type
pe,se,pb,ont`, the last step of `scripts/build_gtdb_database.py`) took about 20 minutes per build on the HPC. Could the
models get pre-reserved slots, or could more cores be used?

## Why it took that long

`--add_model` (HEAD: `Build.h` `AddModel`, `Database.h` `db::Write`):

1. checks every model (`LoadModel`: cPMML parse, contract, knobs): about a second for four r226-size models;
2. writes a new `database.protal.partial` with every member: the ~21 GB of index, reference and tables copied frame by
   frame (one thread, `pread` + `fwrite`), the models compressed at level 19;
3. checks it (`detail::Verify`): reads every copied frame again from both files and compares them (parallel), and
   decompresses the models;
4. renames it over the old file.

To store ~30 MB of models (1.6 MB compressed) it moves the whole file four times: 21.5 GB read, 21.5 GB written, 43 GB
read for the check, **86 GB of I/O**. CPU is not the limit: HEAD's run here took 170 s wall-clock for 26 s of user CPU
time (216 s system time, all of it the copies). At the user's ~20 minutes on the HPC network file system, the file
system delivered ~70 MB/s on average. The 2026-10-01 build profile had ~5.5 min per rewrite on the same file system
(`docs/version_changes.md`, v1: four rewrites in 22 min), so the speed varies with the file system's load. The HPC
logs of the recent builds are not in this repository, so the 20 minutes could not be split further. In the pipeline
log, the step line `adding the ... models` can include waiting for the final database's background build;
`added in H:MM:SS` times `--add_model` alone.

More cores would not change this. The copy moves bytes without compressing them. The check is already parallel. A
parallel copy (`pwrite` at offsets planned in advance) could keep more requests in flight on a parallel file system,
but still moves 86 GB.

## What was changed: the models replaced in place

Pre-reserved fixed-size slots were considered. A model larger than its slot would need the full rewrite again, and the
padding wastes space. Instead, the models are the **last members** of the file, so nothing follows them: they can grow
or shrink, and replacing them rewrites only the end of the file. The one other thing that changes is the directory
(frame 0: each member's first frame and frame count), at the start of the file. It is now stored **raw**: a zstd frame
of uncompressed blocks, whose size depends only on the length of its content. With the same member names, the new
directory is as long as the old one and is overwritten in place.

- `db::detail::RawFrame` writes the directory as a raw zstd frame. It is still an ordinary zstd frame: older protal
  versions read such a file (checked below), and the format version stays 1.
- `ModelsLast` puts the models after the other members wherever a `database.protal` is written: `--build` and
  `--compress_db` (`BundleDatabase`), `--compress_db` adding a gene table (`AddGeneTable`), and `--add_model`'s full
  rewrite.
- `db::InPlaceFrom` decides whether the update can be done in place. Conditions: the same members in the same order;
  those before the first replaced member unchanged (the bundle's own frames); at most 1 GiB from there on (held in
  memory); and a raw directory frame.
- `db::ReplaceTail` does the in-place update:
  1. It reads the members from the first replaced one on: frames to copy are read before anything is overwritten,
     and the others are compressed, all of their frames in parallel (`-t`). A single model is one frame, so it uses
     one thread.
  2. It writes the old directory frame and old tail to `database.protal.journal` and fsyncs them.
  3. It writes the new tail and seek table at the old tail's offset, truncates the file there, overwrites the
     directory frame, and fsyncs.
  4. It checks the result: the members before the first replaced one are where they were, the new frames read back
     as written, and the compressed members read back as their sources.
  5. On any failure it writes the journal's bytes back, so the database is unchanged. On success it removes the
     journal.
- **A run stopped between steps 3 and 5** leaves a file that does not open, with the journal beside it.
  `Options::ResolveDatabase` handles this:
  - With `--add_model`, it writes the journal's bytes back first (`db::RestoreFromJournal`). That happens only if
    the journal belongs to this file: 4 KiB of the untouched bytes before the tail must match.
  - Any other command fails with a message that names the journal and says to run `--add_model` again.
  - `db::Write` removes a stale journal next to a file it replaces.
- `--add_model` (`AddModel`) says which way it went (`Replace the last N of M members ... in place`, or `Rewrite ...`)
  and how long storing took.

A database packed before this change (compressed directory, models before `gene_table.bin`), or one that gains a model
for a read type it had none for (a new member), is rewritten once by `db::Write`, with the models put last. Later
`--add_model` runs on it are in place. The GTDB build's final database is packed by the same binary and already holds
a model or placeholder for every read type (unless `--no-placeholder-models`), so its `--add_model` runs in place.

Readers: an in-place update changes the file under a protal run that is just starting (one that has loaded its
database is unaffected; it reads the members at start-up). `docs/databases.md` says not to start runs on the database
while `--add_model` replaces its models.

## Measurements

[`scripts/bench.sh`](scripts/bench.sh) built a 21.5 GB (20.0 GiB) `database.protal` with
[`scripts/make_bigdb.cpp`](scripts/make_bigdb.cpp): the 2 MB mini database (`data/mini_db`, single file) plus a
20 GiB member of incompressible data in 64 MB frames, the models last. It then ran `--add_model` with the r226 v5
models: `trained_model.xml` 12.6 MB, `_se` 12.6 MB, `_pb` 3.2 MB, `_ont` 3.3 MB, all random forests. The current
gradient-boosted models (250 trees of 63 leaves) should be of a similar size. Each step ran on the file the previous
step left; `-t 6` unless noted. `head` is `15006c2` built from `git archive`
([`scripts/build_head.sh`](scripts/build_head.sh)); `new` is the working tree. Timings and filtered logs are in
[`results/logs/`](results/logs/).

| Run | Binary | What happened | Wall | User | System | Max RSS |
|---|---|---|---:|---:|---:|---:|
| `old_rewrite` | head | full rewrite (adds se, pb, ont) | 169.7 s | 26.3 s | 215.8 s | 0.92 GB |
| `new_first` | new | full rewrite: directory compressed by head, so it cannot be done in place | 148.8 s | 21.8 s | 437.7 s | 0.92 GB |
| `new_in_place` | new | 4 models in place | 9.7 s | 22.0 s | 1.0 s | 0.42 GB |
| `new_in_place_again` | new | 4 models in place | 13.6 s | 36.0 s | 0.9 s | 0.41 GB |
| `new_in_place_1t` | new, `-t 1` | 4 models in place | 20.0 s | 19.8 s | 0.2 s | 0.21 GB |
| `new_one_model` | new | se model in place (se, PB, ONT rewritten) | 8.2 s | 8.0 s | 0.1 s | 0.17 GB |
| `new_one_model_1t` | new, `-t 1` | the same | 8.0 s | 7.8 s | 0.1 s | 0.17 GB |

- In place, `Store the models in database.protal` took 8.8 s of the 9.7 s run. It is level-19 compression of the
  largest model: a 12.6 MB model fits in one frame (frames are 64 MB) and takes ~8 s on one thread, so the four
  models in parallel take as long as the largest one. The file I/O is the 1.6 MB journal, the 1.6 MB new tail, and a 1 KB
  directory frame.
- The local full rewrite (150-170 s) is much faster than on the HPC (~20 min), because the page cache and a local
  NVMe disk serve it. In place, the HPC run should take about as long as here: the ~10 s is CPU, not I/O.
- **Compatibility:** after an in-place update of the mini database, head's `--unpack_db` gives the same files as the
  new binary's (`diff -r`), with the models as given. Head's `--add_model` rewrites the file without complaint.
- **Recovery** ([`scripts/recover.sh`](scripts/recover.sh), [`scripts/make_journal.cpp`](scripts/make_journal.cpp)):
  1. A journal is written for a mini database as it was before an in-place update, and the updated file is cut by
     20 bytes.
  2. A run on it fails: `... the seek table at its end is missing ... (an interrupted protal --add_model left
     db.protal.journal: run --add_model again, which first writes the database's old content back from it)`.
  3. `--add_model` prints `Wrote the old content of db.protal back from db.protal.journal`, then replaces the model
     in place. No journal is left, and the models are those of the state before the stopped run, plus the new one.
  4. A journal whose kept bytes do not match the file is refused (`... is not of other.protal`).
- **Tests:** four new unit tests in `tests/test_Database.cpp`:
  - `TheDirectoryFrameIsRaw` (`ZSTD_decompress` reads `RawFrame` at 1 B to 300 KB);
  - `ModelsAreReplacedInPlace` (frame counts growing and shrinking; bytes before the tail unchanged; the same model
    gives the same file);
  - `InPlaceNeedsTheSameMembersAndARawDirectory` (including a forged file with a compressed directory, as head
    wrote, rewritten once by `Write`);
  - `AnInterruptedReplacementIsWrittenBack`.

  One new end-to-end test: `ReadTypeModelTest.test_models_are_replaced_in_place` (in place, byte-identical with the
  same model and level, the same profile). Results:
  - `Database.*` and `Zstd*`: 35 passed;
  - the whole unit suite without `PackedIndex.*` (13 GB): 367 passed, 3 skipped benches;
  - `ReadTypeModelTest` and `CompressedDatabaseTest` on the mini database: 17 passed.

## What remains

- An existing `database.protal` (any database packed before this change, such as the r226 builds so far) still pays
  for one full rewrite at its first `--add_model`. After that, updates are in place.
- Storing four models now costs about one level-19 compression of the largest model (~8 s). Compressing the models at
  a low level instead would take well under a second (the zstd CLI: 0.07 s at level 3 for the 12.6 MB model, 9.1 s
  at 19) and add about 1-2 MB to a 21 GB file. Not done: 10 s is not worth a second compression setting.
- The full rewrite (first update of an old file, or a new member) is unchanged. Making it cheaper would mean a
  parallel copy, checking the copied frames by checksum instead of re-reading the source (half the check's reads),
  or `copy_file_range`, which depends on the file system. Not done.
