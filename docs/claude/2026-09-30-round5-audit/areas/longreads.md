# Reviewer report: long reads, PacBio HiFi and ONT (round 5, strain-fixes @ 39a8585)

Saved from the reviewer's hand-back. Prebuilt ~/strain-build/bin/protal; an ASan/UBSan -O1 -g build of a copy with
PROTAL_LR_DEBUG prints (patch_debug.py) at ~/audit6/longreads/build-asan/protal. DBs: mini DB copy (only model_pe.xml, so
pb/ont runs used --model_pb/--model_ont mini_db/model_pe.xml: profiles test the plumbing, not a model); a close-relatives
mini DB (~1% apart) ~/audit6/longreads/mini_db_close. Reads from gen_reads.py with a per-base truth map: HiFi 0.1% errors;
ONT 1.6% (30% subs, 25% ins, 45% del, homopolymer losses) and a 5% set; random qualities. eval_sam.py checks records against
truth. -t 1/2. Scripts ../scripts/longreads/ (run1.sh ... run11.sh, lrlib.py, gen_reads.py, eval_sam.py, eval_msa.py,
msa_gaps.py, check_world.py, patch_debug.py, build_asan.sh, chunk_cap.cpp, run_chunkcap.sh), outputs out_e*.txt beside them.

| id | sev | file:line | finding | evidence |
|---|---|---|---|---|
| LR1 | Medium (regression from the performance merge, ff71972) | Options.h:43 (DEFAULT_X_DROP 1000); Alignment/WFA2Wrapper2.h:78-87,106-108; RunProtal.h:231; Core/LongReads.h:289,133-135 | ff71972 ("Pass --x_drop to WFA2") was tested on 150 bp pairs only; with the default --x_drop 1000, WFA2 fails or returns a CIGAR with mismatches inside M blocks for the own species' alignment of a gene an ONT read ends in; that candidate is lost, the relative's weaker alignment is the only hit and gets MAPQ 139-149 (counted unique); with --x_drop 0 the right species wins; the help's "prunes nothing in practice" does not hold for long-read windows | run3b.sh four reads: 1000 -> 2_35 MAPQ 139, 2_38 145, 2_103 149, 3_82 140, no own-species secondary; 0 -> 3_35 26, 3_38 31, 3_103 28, 2_82 23 with the relative as 0x100; debug: "LRDBG align-fail AlignAnchor 3_38" only with 1000; ont_152 "Invalid alignment (298S...55M3S): mismatch G-T inside an M block". run8.sh/run11.sh: e1 ONT 300 reads: 3 wrong-species records, 2 invalid, 1189/1220 partial genes -> x_drop 0: 0, 0, 1193/1220; close-relatives ONT 1 wrong -> 0; debug failed AlignAnchor 28 -> 15; HiFi and 5%-error ONT identical either way; no time difference |
| LR2 | Medium | Options.h:1845-1853 | the length check reads only the first 100 reads: long reads without --read_type (single-end) pass if those 100 are short, then aligned as short reads (one record per read, read soft-clipped, one gene per read; a 140 kb read no record), exit 0, no warning | run4.sh: 100 x 150 bp then 5 x 20 kb, 2 x 70 kb, 1 x 140 kb -> rc 0, e.g. l70_1 16 3_64 ... 21723S1017M...47052S; run4b.sh 5/40 and 2/60 whole genes found; a long read within the first 100 -> rc 30 |
| LR3 | Low | Options.h:1533-1541, RunProtal.h:273-277 | a rerun that reuses a sample's SAM takes the read type from --read_type alone; the SAM's @CO read type is ignored (a PacBio SAM profiled with the ONT model and AF floor, no warning); --profile_only does warn | run7.sh |
| LR4 | Low | Options.h:541,551 | for an all-ONT run the options summary prints the global values, not the ONT ones used (0.85, 0.2); ONT's 0.2 AF floor printed nowhere | e1 ONT log: "max score ani: 0.900000", "snp min af: 0.150000", run line "(-a 0.85)" |
| LR5 | Low | LongReads.h:62 (overlap capped at max_chunk/2), :284 (uint16_t cast) | "every gene lies wholly in its owning chunk" fails silently if the DB's longest gene exceeds ~31.9 kb; window positions wrap for genes over ~59 kb; no guard (theoretical for GTDB markers) | chunk_cap.cpp: longest 31 kb -> 0; 35 kb -> 26/671 (100 kb read), 104/1702 (200 kb); 50 kb -> 181/516, 724/1547 |
| LR6 | Low | LongReads.h:133-135,154-192 | a segment with a single hit gets MAPQv2(best,0), high whatever its length, never settled by the read's other genes: a short end fragment whose own-species candidate was never seeded goes to a relative with high MAPQ | run11.sh m64001_c/379/ccs 3_91 MAPQ 84 61M24991H (read from taxon 2); run8c.sh ont5_267 3_29 MAPQ 94, 141 bp; 1 of 2,758 and 1 of 3,285 records |
| LR7 | Low (code) | RunProtal.h:324, LongReads.h:521 | output handler's CIGAR identity floor hard-coded 0.8 for pb and ont; CigarANI counts a 1 bp indel once, GetProxyANI (-a) twice, so 0.8 never binds when -a >= 0.8 and silently becomes the floor when -a < 0.8; short-read handlers the same | AlignmentUtils.h:175-180,412-424 |
| LR8 | Low (code) | LongReads.h:524-549 | if a segment's best hit fails ExtractSNPs, the next hit is written with the MAPQ computed against the dropped one; ZR:i:1 never read; -m 0 writes every hit (as short reads) | run7.sh: -m 0 and -m 3 3,610 records (1,530 secondary), -m 1 2,080 |
| LR9 | Low (tests) | tests/test_LongReads.cpp; tests/e2e/test_protal_e2e.py:1704-1770 | no unit test of LongReadAligner::operator()/Align or a read over 65 kb; e2e long reads have constant qualities (a QUAL orientation bug on 0x10 would pass); POS/CIGAR never checked against truth | tests |
| LR10 | Low | Hash/IndexCodec.h:203 | UBSan null pointer argument 2 in DecodeChunk on index load (= I/O F14) | every ASan run |

Verified to hold: ~33,000 primary/supplementary records, both strands, random qualities: hard clips + SEQ = read length, H
only at ends; SEQ = read bases (revcomp for 0x10), QUAL reversed for 0x10; no pair/mate/unmapped flags; one primary per read,
rest 0x800; secondaries MAPQ 0; POS/CIGAR put >= 98% (HiFi) / 95% (ONT) of aligned bases at their true position (exceptions:
1 bp ONT indel placement). Chunking (run2.sh; 65,000, 65,001, 65,535, 65,536, 70 k, 130 k, 200 k reads): HiFi 19,363/19,363
whole genes found, none twice, exact coordinates; ONT 5,268/5,268 whole, 2,909/2,910 half-on-read; only wrong record an LR1
case. Taxon consensus (close relatives): HiFi settled 643/702 ambiguous segments, ONT 825/825; every ZR:i:1 on the read's own
species; other wrong-species records: genes the species lacks (marker loss), LR1, LR6. Profiles vs truth within ~0.02:
HiFi .474/.350/.175 vs .493/.348/.159; ONT .494/.259/.247 vs .489/.264/.247; ONT 5% .463/.322/.215 vs .445/.317/.238; close
HiFi .500/.317/.183 vs .497/.299/.205; close ONT .549/.277/.174 vs .529/.282/.189; LowIdentityShare <= 0.05. FASTA input:
Q30 (pb) / Q18 (ont), SAM fields 1-6 and profiles identical to FASTQ. Mixed map (2 pe + 2 pb + 2 ont, READ_TYPE column): each
sample aligned with its own -a (0.9/0.9/0.85) and profiled with its own model; with a placeholder model_PB.xml the pb samples
report nothing, "Accepted 0" in statistics, no MSA row, while pe and ont pass (statistics and strains honour each sample's
cached score). Strain MSAs from long reads (20-24x): HiFi 938/942 and 818/818 true SNPs, 0 and 2 wrong bases; ONT 924/942,
806/818, 4 and 1 wrong; pe 941/942, 817/818, 0; ONT rows 0-2 IUPAC and 20-24 N. --profile_only warns on read-type mismatch;
long read within the first 100 stops the run (rc 30). ASan/UBSan clean on long-read paths (7.4 Mb ONT, 200 kb reads, 20-140 kb
through the short-read path); 19 long-read and option unit tests pass under ASan.
