# Single-end and long reads: evaluation, and single-end support

- **Evaluation**: 2026-09-27, branch `audit-fixes` at `c3c89ae`, by reading the code (nothing run).
- **Single-end support**: 2026-09-29, branch `audit-fixes` at `995c4f1` plus the working-tree
  changes described below (and other uncommitted work in the same tree: training scripts, docs). Built and tested in WSL Ubuntu 24.04 (gcc 13), 8 threads, on the mini
  database of `scripts/mini_db/build_mini_db.sh` (3 species x 3 genomes, synthetic GTDB r226).

## 1. Evaluation (2026-09-27)

The question: how much of protal must change to take single reads of 100-300 bp, and later PacBio
HiFi (~15 kb, ~Q30) and ONT (>50 kb, ~Q17) reads.

**Single-end, 100-300 bp: small to medium.** The profiler already holds mates as optional
(`AlignmentPair`), the SAM reader takes lone records, depth is aligned bases per gene length, and
the SNP/MSA code does not depend on read length. Missing were the command line and map file
(paired-end only), the run loop (a single-end branch that called `exit(8)`, over stale code), and
a working single-end SAM writer: the old `ProtalOutputHandler` flagged every record unmapped (so the
profiler would have skipped all), inverted the reverse-strand flag, wrote the raw negative score as
MAPQ, and did not check CIGARs. The model was trained on pairs (`total_hits` counts mates, MAPQ
comes from both mates' scores, mate rescue adds hits), so single-end reads need their own model.

*Correction:* the evaluation also said that the single-end reader's `LoadBlock` resynchronises on
the next line starting with `@` and so breaks on quality lines that start with `@` (Q31). Reading it
again, it looks at the line after that one and handles the case. It was still replaced (section 2),
because it takes 1 kB per lock, about four reads.

**PacBio HiFi: large.** Seeding and alignment accuracy are not the problem (~98.5% of 15-mer seed
cores survive Q30). protal assumes that a read maps to one locus. With ~120 markers per genome
(bac120/ar53), a dozen or more of them in the ribosomal-protein operons within 10-15 kb, a long read
often spans several indexed genes. Today only the top `--align_top` (3) anchors of a read are aligned
(`AlignmentStrategy.h`), MAPQ ranks a read's gene hits against each other (correct hits on several
genes give each other MAPQ ~0), one primary record is written, and the profiler groups a read's
records as candidates of one locus. Needed: chains per query segment, alignment and MAPQ per
segment, primary + supplementary (0x800) + secondary records per segment, and profiler grouping by
read and segment. Also: align (and write, hard-clipped) only the read's window on the gene instead
of the whole read (every record repeats the full SEQ/QUAL); scale the seeding cap (128 seeds, taken
by ubiquity not position) and the 9 bp end margins with read length; widen the `uint8_t` unique-seed
counts of anchors and alignments (a whole-gene anchor exceeds 255; they feed the ZU tag and the
lu/lsu features); train a model on HiFi simulations (pbsim3).

**ONT: HiFi plus more.** Read positions and lengths are 16 bit in `LookupResult`, `ChainLink`,
`ChainAlignmentAnchor::total_length` and the seed extension limits, `AlignmentInfo` counters are
`uint16_t` and `CigarInfo` counters `int16_t`: reads over 65,535 bp (common at 50 kb N50) wrap.
Seeds are chained only within +-6 of the first seed's diagonal, and the alignment window is
projected from that seed, which ~2% indel-heavy errors break over a 1-2 kb gene: needs DP chaining
(minimap2-like) projected from both chain ends. Scoring (4, 6, 2) with the 0.9 identity cut charges
a 1 bp indel twice a mismatch, which puts ONT reads of a 95-97% ANI strain at the cut. SNP calling
needs an allele-frequency floor (default 0), homopolymer-indel suppression and recalibrated quality
gates; the depth filter's identity margin assumes a species' reads are close in identity. The index
and database need no change for any of the three (~74% of 15-mer cores survive 2% error, ~2 in 9
positions sampled).

## 2. Single-end support (2026-09-29)

As the user asked: general single-end support in protal, with the single-end model assumed to be
`model_se.xml` in the database; training that model is a separate, later step.

### What changed

| Where | Change |
|---|---|
| `Options.h` | `-1` without `-2` is single-end; a map's `SECOND` may be `-` (that sample is single-end) or absent (all are). Per-sample read type (`ResolveReadTypes`): no second file, or with `--profile_only` the SAM's first usable record (no 0x1 = single-end; `HoldsPairedReads`). Default prefix of a single-end file: its name without FASTQ/FASTA and compression extensions. `--model_se`; `ModelDbFile(single_end)`: `--model_se`, else `--model`, else `model_se.xml` (single-end) / `model.xml` (paired-end). Up-front check that each read type in the run has its model. `PairedMode()` removed |
| `SeqReader.h` | `SeqReaderSE`: batches of 32 records per lock, as `SeqReaderPE` |
| `Classify.h` | `RunSingleEnd` replaces the stale `Run`: the paired loop without mate rescue or pairing. The runtime and histogram files in `misc/` are written by one function for both |
| `AlignmentOutputHandler.h` | `ProtalSingleOutputHandler<DEBUG>` replaces `ProtalOutputHandler`: candidates ranked by bitscore, each alignment once (the same taxon, gene, position and strand from two anchors would give MAPQ 0), MAPQ as for a mate aligned alone (`MAPQv2`, best vs second), FLAG 0/16 (+256 after the first), TLEN 0, CIGAR check per candidate, a read's records written as one buffer unit, `-m` honoured |
| `SamHandler.h` | an unpaired record is returned in the first slot (was: second); `PairedRecords()`; `HoldsPairedReads` |
| `RunProtal.h` | one loop over samples, each aligned as single- or paired-end; both models loaded up front if the samples need them; `ProfileWrapper` scores each sample with its model. Later stages (statistics, strains) use the scores cached per taxon, so a run may mix read types |
| `Build.h` | `--build` / `--compress_db` pack a `model_se.xml` found in the database folder; `--unpack_db` writes every member already |
| docs | `running.md` (single-end section, `--model_se`), `database-files.md`, `model-training.md`, `testing.md` |

Unchanged: seeding, chaining, alignment, the profiler's read processing, SNP calling, strains.

### Tests

Unit tests (`just test`), new: `SamReader.SingleEndReadsAreFirstReads`,
`SamReader.TellsPairedFromSingleEndReads`, `SingleOutputHandler.*` (orientation and flags, ranking
and secondaries, `-m`, duplicate candidates, an inconsistent candidate), `SampleMap.SingleEndSamplesHaveNoSecondFile`,
`ReadFileStem.*`, `Options.EachReadTypeHasItsModel`, `Options.ReadTypesOfSamples`.

End-to-end (`tests/e2e/test_protal_e2e.py`, class `SingleEndTest`): the first mates of the
simulated samples as single-end reads, on a copy of the database with its `model.xml` as
`model_se.xml` (so the profiles test the plumbing, not a model): unpaired records with SEQ in
reference orientation, profiles and joint MSAs, `--profile_only` choosing the single-end model from
the SAM, the missing-model check before aligning, a map mixing both read types, the default prefix,
and `model_se.xml` inside `database.protal`.

Commands (run from scripts in the session scratch space; the checkout is rsynced to the WSL
filesystem, as on /mnt/c the case-insensitive NTFS lets cPMML include protal's `Options.h`):

```bash
rsync -a --delete --exclude /data --exclude '/build*/' --exclude '/cmake-build-*/' --exclude /.git <checkout>/ ~/protal-se-build/src/
cmake -S ~/protal-se-build/src -B ~/protal-se-build/build -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON
cmake --build ~/protal-se-build/build -j 8 && ctest --test-dir ~/protal-se-build/build
PROTAL=~/protal-se-build/build/protal bash scripts/mini_db/build_mini_db.sh ~/protal-se-build/mini_db
PROTAL_TEST_DB=~/protal-se-build/mini_db/protal_db PROTAL=~/protal-se-build/build/protal \
    SIMULATE=~/protal-se-build/build/simulate_metagenomes python3 -m unittest -v tests/e2e/test_protal_e2e.py
```

Results: unit tests 111 of 111 passed (12 new). End-to-end 65 tests: the 57 existing ones passed
unchanged (paired-end behaviour is as before); of the 8 new `SingleEndTest` ones, one failed on a
wrong output path in the test itself (a map without a SAM column writes `<OUTPUT_DIR>/<PREFIX>.sam`),
and after fixing that the class passed 8 of 8. Final run on `995c4f1` plus these changes, with a
freshly built mini database: unit tests 111 of 111, end-to-end 65 of 65.

Paired-end vs single-end on the same reads: `scripts/mini_db/simulate_reads.py --community
examples/mini_db/community.tsv --pairs 30000 --seed 7` (reads from whole genomes, 0.2%
substitutions), aligned as pairs and as their R1 files alone, `-t 8`, the database's `model.xml`
also as `model_se.xml`:

| | primary records | MAPQ >= 4 | Mockella alpha / beta / Fakibacter gamma |
|---|---|---|---|
| truth | | | 0.500 / 0.300 / 0.200 |
| paired-end | 29,303 (14,614 read1, 14,689 read2) | 99.9% | 0.499 / 0.302 / 0.199 |
| single-end (R1) | 14,614 | 99.7% | 0.503 / 0.297 / 0.200 |

The same R1 reads align in both modes, and the profiles agree. The mini database's genomes are the
reference itself, so this shows that single-end reads go through protal correctly; it does not
show what mate rescue is worth on divergent strains, nor how the paired-end model scores
single-end samples of real complexity.

### Open

- `model_se.xml`: to be trained on single-end simulations (the training pipeline simulates pairs).
  Until then a database has none and single-end samples need `--model_se`.
- The website's usage and map-file pages describe paired-end reads only.
- Single-end reads lose mate rescue: expect fewer hits of divergent strains, and more MAPQ-0 reads
  between close relatives, than for the same reads as pairs. Not measured yet.
