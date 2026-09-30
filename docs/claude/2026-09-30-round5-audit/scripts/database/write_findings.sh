#!/usr/bin/env bash
# Persist findings.md into the WSL work dir (caller requested), and free large raw-index copies.
set -u
mkdir -p ~/audit6/database
cat > ~/audit6/database/findings.md <<'EOF'
# Audit round 5 — database container and compression

Reviewer area: src/Utilities/Database.h, src/Utilities/Zstd.h, src/Hash/IndexCodec.h,
Seedmap load/save, src/Build.h (bundle/compress/unpack/add_model), db location/validation in
src/Options.h, docs/database-files.md. Repo: worktree strain-fixes @ 39a8585. All experiments
used at most 2 threads. Mid-audit the shared WSL VM crashed (Wsl/Service/E_UNEXPECTED, machine
load ~13 from other sessions) and rebooted; ~/audit6 persisted and every experiment completed.
wsl --shutdown was never used.

## Findings

M1 (Medium) Zstd.h ~829 ReadFrame / ~730-735 ParallelReadFrames — seek-table decompressed_size
(<=4GB) sizes the output buffer BEFORE validation. protal's frames store the content size so the
ZSTD_getFrameContentSize guard catches a forged size; but a frame with content size omitted
(legal zstd, e.g. `cat f | zstd`, pzstd) skips the guard and the code resizes to the claimed size.
A tiny crafted file can hold many content-less frames each claiming ~4GB; only compressed sizes
are cross-checked (compressed==table_start), not decompressed. Per worker => x threads.
Evidence: harness.cpp under `ulimit -v 1200000`: Case A (size stored) -> "holds 6370 bytes, the
seek table says 4294967280; out buffer resized to 0" (caught, no alloc). Case B
(ZSTD_c_contentSizeFlag=0) -> "header content size = UNKNOWN ... ReadFrame THREW: std::bad_alloc".

M2 (Medium) Database.h ~291-294 Locate; Options.h ~1604-1628 — a partial/aborted --unpack_db (or
any stray index.prx/.zst) next to a valid database.protal makes Locate pick separate-file mode and
silently ignore the complete bundle; a profiling run then emits a bare "... does not exist"
cascade that never mentions the bundle (the unused_bundle remedy hint is only added under --build).
--unpack_db writes the index member first, so a crash after that rename yields exactly this state.
Evidence: s13 case C: folder {index.prx.zst, database.protal} only -> exit 30, "Sequence file
does not exist ... reference.fna", ".map", "internal_taxonomy.dmp", model, unique_kmers — no
mention of the present database.protal. (Case D, all files + bundle: loads separate, prints
"database: separate files ... (not .../database.protal)".)

L1 (Low) Build.h AddModel ~444-462 — folder-mode --add_model writes only the separate model file
and leaves a co-located database.protal untouched/stale (no RemoveStaleBundle, no warning). If the
separate files are later removed, the stale bundle (missing the new model) is used silently.
Evidence: s08 B: bundle md5 unchanged before/after; folder gains model_se.xml; bundle still 0
model_se members.

L2 (Low) Options.h ~1633 — compression-option range checks are guarded by
((m_build&&m_compress)||m_compress_db), so --add_model skips them. Benign (zstd clamps level; bad
frame size errors later) but inconsistent: --build rejects --compress_level 99, --add_model accepts
it (exit 0). Evidence: s08 E (level 99 -> stored) vs D (frame_mb 0 -> clean "needs frames").

L3 (Low) Zstd.h CompressFramesTo ~1144-1160 / WriteSeekable ~1048 — peak compress memory ~=
threads x frame_size (one in+one out buffer per worker). --compress_frame_mb 4095 -t 16 ~ 100+ GB.
Within documented ranges; GTDB-scale footgun. Code review.

L4 (Low) Zstd.h ~793 — file_size(path, ec) result used without checking ec in the raw ParallelRead
path; on failure total=SIZE_MAX drives the chunk loop (reads then fail). Bounded per worker to
64MB; benign latent trap. Code review.

L5 (Low) Zstd.h CompressFile verify ~1279-1298 — the verify pass of CompressReference /
--compress_db --no_bundle decompresses+compares the whole reference.fna single-threaded after the
parallel compress (serial tail at GTDB scale). Intentional safety, not a bug. Code review.

L6 (Low) Options.h ModelDbFile ~762-765 — --model NAME resolves an existing CWD path before the DB
member of that name; a CWD file named like a member shadows the DB model. Documented precedence.

## Verified to hold

- Round-trip deterministic + byte-identical both ways: two --compress_db runs byte-identical;
  --decompress_db reproduces index.prx/reference.fna/map/taxonomy/unique_kmers/model exactly
  (s07). Confirms the docs' byte-identity claim.
- Atomicity: all writers use <path>.partial + verify + rename; leftover *.partial ignored by
  Locate, folder still loads+profiles (s07 #4); failure => "the database is unchanged", no partial
  (s08 D); success => no partial (s08 A).
- Path traversal blocked: forged names "../escapee.map", "/tmp/escapee.map" rejected by IsFileName
  in Bundle::Open; nothing created outside unpack dir. Duplicate/non-tiling/non-covering members,
  bad version, and directory>16MB (rejected pre-decompression) all -> exit 30, clear messages
  (s09, forge.py).
- Content checksum catches a 1-byte flip in a member frame -> exit 8 "Restored data doesn't match
  checksum". Seek-table truncation -> exit 30 "seek table at its end is missing"; mid-frame
  corruption -> exit 8 (s13, existing tests).
- Wrong --db (random file, plain reference.fna.zst) -> exit 30 "neither a folder nor a single-file
  protal database" (s09).
- DecodeChunk ASan+UBSan fuzz (fuzz_decode.cpp): ~4.2M mutated/truncated inputs into exactly-sized
  buffers -> decoded_ok=3589241 returned_error=623295, no ASan/UBSan abort, no crash (s12b).
- 33/33 existing unit tests pass: Database.*, Zstd.*, ZstdSeekable.*, IndexCodec.* (s03).
- --add_model bundle mode adds member, rewrites via .partial + verify, prints coverage, no leftover
  partial; rewrites the whole file (transient ~2x bundle size on disk) (s08 A).

## Docs vs behaviour
database-files.md accurate for tested workflows. Gaps: (1) folder-mode --add_model does not
update/remove a co-located database.protal (L1); (2) incomplete separate files next to a bundle
shadow it with unhelpful errors (M2).

## Scripts (this folder)
s00_env, s01_mem, s02_survey, s03_unit, dbinfo.py, s06_setup, mkreads.py, s07_roundtrip,
s08_addmodel, forge.py, s09_forge, harness.cpp+s11_harness, fuzz_decode.cpp+s12_fuzz+s12b_run,
s13_protal_corrupt.
EOF
echo "wrote ~/audit6/database/findings.md ($(wc -l < ~/audit6/database/findings.md) lines)"
# free the big raw-index copies (3.2 GB each), keep the small DBs and scripts
du -sh ~/audit6/database 2>/dev/null
rm -rf ~/audit6/database/raw ~/audit6/database/rtA ~/audit6/database/rtB ~/audit6/database/rtC ~/audit6/database/out_* 2>/dev/null
du -sh ~/audit6/database 2>/dev/null
