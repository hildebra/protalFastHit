# Reviewer report: I/O layer after the performance merge (round 5, strain-fixes @ 39a8585)

Saved from the reviewer's hand-back (it could not write findings.md). Mini DB unpacked at WSL ~/audit6/io/db; ASan/UBSan
Debug copy ~/audit6/io/src/build-asan; outputs ~/audit6/io/t1..t10; scripts ../scripts/io/ (lib.sh, mkreads.py,
setup.sh, build_asan.sh, t1_formats.sh ... t10_strains.sh). WSL crashed once ("Catastrophic failure"), restarted with 24 GB;
before that the OOM killer killed a protal while ~23 protal processes ran (most not the reviewer's).

| id | sev | file:line | finding | evidence |
|---|---|---|---|---|
| F1 | Medium | src/SequenceUtils/SeqReader.h:205-208 | paired-end only: neither file parses as FASTQ/FASTA -> "sequence reader - unrecognized file format" but exit 0 with a header-only SAM and empty profile; LoadBatch sets the error but Next() returns on "both blocks invalid" before checking Error(); hits .fq.zst, .bz2, .xz, BAM. Single-end fails (rc 1). Pre-existing | t3b.sh: -1 zst_R1.fq.zst -2 zst_R2.fq.zst -> rc 0, 0 records; same for .bz2 and non-FASTQ text; se_* rc 1 |
| F2 | Medium (regression, 401694b) | src/IO/ThreadedGzStream.h:51, src/IO/Bgzf.h:139-144 | FIFO and process-substitution input broken: open() sniffs for BGZF with a separate ifstream, eating the first 18 bytes of a pipe; paired then hits F1 (rc 0, empty), single-end rc 1. Before 401694b open() only called gzopen | t9b_fifo_se.sh: -1 <(zcat a1.gz) -2 <(zcat a2.gz) -> rc 0, 0 records; mkfifo plain and gzip PE rc 0 empty; SE rc 1 |
| F3 | Medium | src/RunProtal.h:320,334,358; src/IO/FastxReader.cpp:124 | an unreadable input (permission denied) -> rc 0, empty SAM and profile, no message: ThreadedGzIstream only sets badbit, LoadBatch takes !ifs as end, RunProtal never checks is_open(); paired and single, plain, gzip, BGZF | t9_special.sh chmod 000: nr, nrse, nrgz, nrb all rc 0, 0 records, empty log |
| F4 | Medium | src/IO/ThreadedGzStream.h:180-192 (zlib-ng gz_look) | after a complete gzip member, anything that is not a gzip header is ignored as trailing garbage, no error: a damaged later member header drops all reads after it silently (gzip -t warns, exit 2) | t3b.sh: mbad_R1.fq.gz SE rc 0, 1962 records (full 3924); both mates damaged PE rc 0, 3924 (full 7848); one mate -> rc 1 via the pair count |
| F5 | Low | src/IO/FastxReader.cpp:16-19 (callers 240, 300, 306, 315, 324, 328) | StripString calls str.back() on an empty string: FASTA with a blank or \r-only line -> ASan heap-buffer-overflow (1 byte before); also --build's reference reading; release reads correctly but UB. Pre-existing | t6_asan.sh: fablank, facrlf, fablse each 1 ASan error in StripString <- ReadNextSequence; release 400/400 records |
| F6 | Low | reads never uppercased (only genes: GenomeLoader.h:209,507) | lowercase reads seed (k-mer code accepts acgt) but never align (raw character compare): rc 0, 0 records, empty profile; anchored and whole-read. Pre-existing | t2 lower: rc 0, 0 records (baseline 7848); t3b lower_whole 0 |
| F7 | Low | src/IO/FastxReader.h:40 (IsValidPair unused) | mismatched R1/R2 names never checked: a shuffled R2 -> rc 0, every pair discordant, no warning. Pre-existing | t2 names: rc 0, 7848, QNAME "1/1"; shuf: flags 65/129/113/177, no warning |
| F8 | Low | src/IO/FastxReader.cpp:221-227 | a plain single-end FASTQ whose last record is cut short loses it silently | t6 cut_R1.fq SE: rc 0, 99 names vs 100; PE rc 1 |
| F9 | Low | src/IO/FastxReader.cpp:223-229 | SEQ/QUAL lengths and the + line not checked: a short quality gives an invalid SAM record; profiling skips it | t6 qshort: rc 0; profiler "skipped 1 record(s): ... QUAL and SEQ lengths differ" |
| F10 | Low | src/IO/ThreadedGzStream.h:236-240; src/RunProtal.h:402-404 | BGZF followed by an ordinary gzip member (valid gzip) rejected as "truncated or corrupt"; read_error_message() never printed for any gzip failure | t2 bgzfplusgz, t3: rc 1, generic message; gzip -t passes; reverse order reads fine |
| F11 | Low | src/IO/SamFile.h:315-319, src/RunProtal.h:294 | two runs on the same prefix overwrite each other's fixed-name <sam>.partial (O_TRUNC, no lock): a corrupt SAM can land under its final name; later reruns skip it and fail at profiling | t5_crash.sh (kill -9 hit the timeout wrapper): second run "cannot read ....partial.records.partial"; final big.sam.gz fails --profile_only "no BGZF block at byte 4944725" |
| F12 | Low (code) | src/RunProtal.h:107-112, src/IO/SamFile.h:281 | SAM not fsynced before rename: after power loss the final name can point at a short file (.zst/.gz caught by end markers; a plain .sam cut at a line boundary not) | code |
| F13 | Low | src/IO/FastxReader.cpp:133; src/SequenceUtils/GenomeLoader.h:602,839 | empty inputs pass without warning (header-only SAM, empty profile); a 0-byte unique_kmers.tsv -> "No taxon passes the model", rc 0; a 0-byte reference.map accepted, the error then blames unique_kmers.tsv | t2 empty/emptygz, t3 se_empty/se_emptygz: rc 0; t8 u_uempty rc 0; t8 empty map rc 8 with the unique_kmers message |
| F14 | Low (outside I/O) | src/Hash/IndexCodec.h:203 | UBSan: memcpy with a null source (empty vector .data(), size 0) on every index load; with halt_on_error=1 stops every sanitizer run | t6 |

Suggested fixes: F1 check Error() after LoadBlockOMP(); F2 skip the BGZF sniff for non-regular files or sniff through
the one descriptor; F3 fail the sample if !is.rdbuf()->is_open(); F4 inflate members with zng_inflate, require a header or
EOF after each; F5 guard !str.empty(); F6 uppercase reads; F10 print read_error_message().

Verified to hold: sorted SAM records identical for -t 1/-t 2 x sam/gz/zst (7848, md5 f56e37a7, 329 header lines), -t 1
unsorted order identical across formats; @SQ exactly the 327 genes named; profiles identical; two samples with MSAs: every
output identical for -t 1 sam, -t 2 zst, -t 2 gz and --profile_only of the zst SAMs; gzip -t and zstd -t accept outputs;
reruns reuse each format (mtime unchanged), a different --sam_format reuses the plain .sam; --profile_only same profile per
format and with --full_sam_header; inputs matching baseline: CRLF, no trailing newline, gzip, multi-member gzip, BGZF,
gzip+BGZF, two BGZF concatenated, trailing blank line, FASTA single/multi-line, 3 MB read line (plain, gzip, BGZF);
failures detected (rc 1, no SAM): gzip truncated, flipped byte, zeroed CRC, BGZF truncated mid-block or missing EOF block,
BGZF CRC, R1/R2 count mismatch (1, 32, 33), one empty mate, blank line mid-file, a directory, non-FASTQ SE; SAM reading:
.sam.zst cut/missing seek table/frame boundary/dropped frame/flipped byte rc 1 with reason, .sam.gz cut/missing EOF/flip,
plain cut mid-line rc 1; 0-byte SAMs "the file is empty"; plain text named .zst/.gz, single-member gzip, zstd CLI output
read correctly; wrong LN or foreign gene in @SQ rejected; a SAM with no aligned genes works in all formats. Undetectable by
design: plain SAM or multi-member gzip cut exactly at a line/member boundary. Crash safety: kill -9 during alignment leaves
only <sam>.partial.records.partial, rerun realigns (313920 records, all formats); EFBIG via ulimit -f in each format gives a
reason, rc 1, nothing left behind; unwritable -o rc 2. Gene arena: identical records and profiles from raw reference.fna,
database.protal, single-frame zstd, seekable 997-byte frames (genes span frames), --preload_genomes_off, 1 and 2 threads;
lifetime OK (arena owned by GenomeLoader in ProtalDB, never copied/moved after load). Parallel parsers: 1.19 MB
reference.map / 1.34 MB unique_kmers.tsv, -t 1/2, errors with correct physical line numbers incl. CRLF, blank lines, last
line without newline, for several error kinds. Sanitizers: PE runs in all formats and --profile_only clean except F14; 46 I/O
unit tests pass under ASan+UBSan with leak detection.
